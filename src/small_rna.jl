module small_rna

#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#
# Small RNA pipeline
#
# This script can be used to analyze Small RNA V3 and V4 libraries
#
#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#<><><><><><>#
using CairoMakie
using CSV
using DataFrames
using ProgressMeter
using MultivariateStats
using StatsBase
using StatsPlots
using Statistics
using MLBase
using GLM
using Measures
using GZip
using UMAP
using Clustering
using Distances
using PDFmerger
using XAM
using ArgParse

export MissingReferenceError
export Config
export get_genus_names
export create_target_organism_fasta_file
export capture_target_files
export progress_bar_update
export get_read_q_score!
export parse_fastqs
export trim_adapters
export generate_mirna_counts
export mirna_discovery_calculation
export plot_mirna_counts
export calculate_read_length_distribution
export plot_fragment_lengths
export align_with_bowtie2
export calculate_metrics
export plot_metrics
export find_common_mirnas
export write_common_mirna_file
export plot_clustering
export cq_vs_mirna_count_correlation
export plot_cq_linear_regression
export remove_files
export remove_intermediate_files

"""
    struct MissingReferenceError
        message::String
    end

A structure that represents a custom exception for missing reference.
"""
struct MissingReferenceError <: Exception
    message::String
end

"""
    struct Config
        library::Union{String, Nothing}
        fasta::Union{String, Nothing}
        organism::Union{String, Nothing}
        threads::Int
    end

A structure that stores the command-line arguments.

# Functionality

- Contains a check to ensure that a miRNA fasta file is passed as an argument. It is was
not, an error will be returned. This is to reduce the chances of an alignment errors downstream.
"""
struct Config
    library::Union{String, Nothing}
    fasta::Union{String, Nothing}
    organism::Union{String, Nothing}
    threads::Int

    # Constructor with error checks, Union stores possibility of no argument being passed
    function Config(library::Union{String, Nothing}, fasta::Union{String, Nothing}
                    , organism::Union{String, Nothing}, threads::Int
        )

        # Check if user specified a real miRNA file
        if !isfile(fasta)
            throw(MissingReferenceError("No miRNA reference found"))
        end
        
        new(library, fasta, organism, threads)
    end
end

"""
Create a file with a list of unique species names given a miRNA mature or hairpin fasta.
"""
function get_genus_names(fasta::String)
    genuses = Set{String}()

    open(fasta) do miRNA_reference
        for line in eachline(miRNA_reference)
            if startswith(line, ">")
                genus_abbreviation = line[2:4]
                push!(genuses, genus_abbreviation)
            end
        end
    end

    return genuses
end

"""
Function takes in mirBase mature miRNA fasta file with data from all available
organisms and pulls out only the miRNA data pertaining to the target organism.
"""
function create_target_organism_fasta_file(organism_name::String, fasta_file::String)

    # Skip making species miRNA reference if it exists in data/
    bowtie2_reference_files = capture_target_files(".bt2", "data")
    if isempty(bowtie2_reference_files)

        # Create file with all genuses with miRNA annotations
        genuses = get_genus_names(fasta_file)

        # Set the target or organism's genus and species for labeling the output file
        target_genus_abbreviation = ""
        output_fasta_file_name = organism_name * ".fa"

        # Loop through file containing the genus of all the organisms with miRNA data
        # If the target organism is in that list, the genus abbreviation is stored for use below
        for genus in genuses
            if occursin(lowercase(organism_name), lowercase(genus))
                target_genus_abbreviation = genus
                break
            end
        end

        # Open the IO streams and sets the output file names for the fasta and bowtie2 index
        input_fasta_file = open(fasta_file)
        output_fasta_file = open(string("data/", organism_name, ".fa"), "w")

        #=
        Loop through miRNA fasta looking for the header and sequence information
        for the target genus. If there's a match, that information is added to a new fasta file.
        Must convert RNA sequences in miRNA fasta to DNA.
        =#
        for line in eachline(input_fasta_file)
            input_genus_abbreviation = match(r">\K\w+", line)
            if !isnothing(input_genus_abbreviation)
                if target_genus_abbreviation == input_genus_abbreviation.match
                    write(output_fasta_file, line, "\n")
                    write(output_fasta_file, string(replace(readline(input_fasta_file), "U" => "T"), "\n"))
                end
            end
        end

        close(input_fasta_file)
        close(output_fasta_file)

        #=
        Use the newly created single organism fasta to build a bowtie2 index; will be used
        for alignment with bowtie2 to determine miRNA count information.
        =#
        wait(run(pipeline(
            `bowtie2-build 
            data/$output_fasta_file_name
            data/$organism_name`
            , devnull)
            , wait = false
            )
        )
        
    end
end

"""
    capture_target_files(files_to_capture::AbstractString, directory::AbstractString=".")

List all files in a directory. 

Check to see if each file contains the target file's string. 

# Example
```julia
julia> capture_target_files(".txt")
3-element Vector{String}:
 "file1.txt"
 "file2.txt"
 "file3.txt"
```
"""
function capture_target_files(files_to_capture::AbstractString, directory::AbstractString=".")
    return filter(file -> occursin(files_to_capture, file), readdir(directory))
end

"""
Update progress bar on the command line.
"""
function progress_bar_update(number_of_records::Int64
                            , interval::Number
                            , description::String
                            )
    progress_bar_update = Progress(number_of_records
    , dt=interval
    , barglyphs=BarGlyphs("[=> ]")
    , barlen=100
    , desc=description
    )

    return progress_bar_update
end

"""
    get_read_q_score!(line::String, q_score_list::Vector{Number})

Calculate the quality score for a given quality string.

This function accepts ASCII-encoded quality score and produces the average q score 
using a scalar value that estimates the 'total error probability' for that record.
It means that if you want to calculate average quality score, you don't just sum
up all the phred scores and find the average, but rather consider the 
probability distribution of the scores.

The Phred score for a base `Q` is calculated as `Q = ASCII value of quality score - 33`.
The error probability `P` for that base is then calculated as `P = 10^(-Q/10)`.
Then these probabilities are summed up for all the bases to derive total error.

# Arguments

`q_score` - The quality scores of the current line.
`q_score_list` - The vector containing all calculated quality scores.

# Returns

The updated q score list the line's average quality score.

# Reference

Illumina's explanation on quality scores and their corresponding error rates:
<https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html>

# Example
```julia
julia> get_read_q_score!("FGFEEE<FC<FGGGGGGGGGGGFFGFG8<@8CFFGF8EFGGCGFGGGGGGFG", [36.2, 35.9])
3-element Vector{Float64}:
 36.2
 35.9
 35.7
```
"""
function get_read_q_score!(line::String, q_score_list::Vector{Float64})
    probability_sum = 0

    # Find average read quality by converting quality encodings to Unicode numbers
    for char in codeunits(line)
        phred = char - 33
        error_probability = 10^(-phred / 10.0)
        probability_sum += error_probability
    end

    prob_q_score = probability_sum / length(line)
    q_score = -10.0 * log10(prob_q_score)

    # Add Q score to list
    push!(q_score_list, q_score) 

    return q_score_list
end

"""
    parse_fastqs(fastqs::Vector{String}, sample_names::Vector{SubString{String}}, library::String)

Determine the library type, calculate read counts, percent dimers, and average Q scores 
for a given set of fastq files.

# Parameters
- `fastqs`: A vector of filenames for the fastq files.
- `sample_names`: A vector of sample names corresponding to the fastq files.
- `library`: The type of library ("v3" or "v4").

# Returns
Three dictionaries containing read counts, percent dimers, and average Q scores 
indexed by sample name.

# Example
```julia
julia> parse_fastqs(["sample1.fastq.gz", "sample2.fastq.gz", "sample3.fastq.gz"], 
                    ["sample1", "sample2", "sample3"], "v4")
# Returns:
# - A dictionary of read counts.
# - A dictionary of percent dimers.
# - A dictionary of average Q scores.
```
"""
function parse_fastqs(fastqs::Vector{String}
                        , sample_names::Vector{SubString{String}}
                        , library::String
                        )

    number_of_records = length(fastqs)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .5
                                            , "Parsing fasta file..."
                                            )

    # Skip first 8 bp for v3 libraries, which would be 4N + 4N
    library_dimer_starts = Dict("v4" => 1, "v3" => 9)
    dimer_sequence = "TGGAATTCTCGGGTGCCAAGG"

    read_count_dict = Dict{String, Int64}()
    dimer_count_dict = Dict{String, Float64}()
    q_score_dict = Dict{String, Float64}()
    
    for (fastq_file, sample_name) in zip(fastqs, sample_names)
        q_score_list = Vector{Float64}()
        dimer_count = 0
    
        fastq = GZip.open(fastq_file)
        sample_read_count = 0
    
        # Should use FASTX here, but there seems to be incompatibility with the other
        # packages when creating apps using PackageCompiler.jl
        for lines in Iterators.partition(eachline(fastq), 4)
            sequence_line = lines[2]
            quality_line = lines[4]
    
            dimer_start = library_dimer_starts[library]
            if startswith(sequence_line[dimer_start:end], dimer_sequence)
                dimer_count += 1
            end
    
            # Calculate average quality score
            get_read_q_score!(quality_line, q_score_list)
            sample_read_count += 1
        end
        close(fastq)
    
        read_count_dict[sample_name] = sample_read_count
        q_score_dict[sample_name] = round(mean(q_score_list), digits=2)
        dimer_count_dict[sample_name] = 100 * dimer_count / sample_read_count

        next!(update_progress_bar)
    end
    
    return read_count_dict, dimer_count_dict, q_score_dict
end    

"""
    trim_adapters(fastqs::Vector{String}, library::String, sample_names::Vector{SubString{String}})

Trim the 3' adapter from each read in the provided fastq files based on the given library type.

# Library Adapters:
- V3 Read 1:
   5' Adapter: GATCGTCGGACTGTAGAACTCTGAACNNNN 
   3' Adapter: NNNNTGGAATTCTCGGGTGCCAAGG

- V4 Read 1:
   5' Adapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGA
   3' Adapter: TGGAATTCTCGGGTGCCAAGG 

3' bases with a quality score below 20 are trimmed. Reads shorter than 16 bases or 
those that weren't trimmed are discarded.

# Parameters
- `fastqs`: A vector of filenames for the fastq files.
- `library`: The type of library ("v3" or "v4").
- `sample_names`: A vector of sample names corresponding to the fastq files.

# Returns
A list of filenames for the trimmed fastq files.

# Example
```julia
julia> trim_adapters(["sub_sample1.fastq", "sub_sample2.fastq", "sub_sample3.fastq"], 
                    "v4", 
                    ["sample1", "sample2", "sample3"])
# Returns:
["sample1.cut.fastq", "sample2.cut.fastq", "sample3.cut.fastq"]
```
"""
function trim_adapters(fastqs::Vector{String}, library::String, sample_names::Vector{SubString{String}})
    number_of_records = length(fastqs)
    update_progress_bar = progress_bar_update(number_of_records, .5, "Trimming adapters...")

    adapter_mapping = Dict(
        "v4" => ["--adapter", "TGGAATTCTCGGGTGCCAAGG"],
        "v3" => ["--cut", "4", "--adapter", string("N", "{4}", "TGGAATTCTCGGGTGCCAAGG")]
    )

    min_length = 16

    for (fastq_file, sample_name) in zip(fastqs, sample_names)
        # Building up the command manually to account for different library types
        base_command = ["cutadapt", "--cores=0", "--discard-untrimmed", "--quality-cutoff", "20"]
        output_commands = ["--output", "$sample_name.cut.fastq", "--minimum-length", "$min_length"
                            , "--too-short-output", "$sample_name.short.fastq", fastq_file
        ]
        
        cutadapt_command = Cmd(vcat(base_command, adapter_mapping[library], output_commands))

        wait(run(pipeline(
            cutadapt_command, 
            stdout="$sample_name.cutadapt_information"
            ), wait=false
            )
        )

        next!(update_progress_bar)
    end

    trimmed_fastq_files = capture_target_files(".cut.fastq")

    return trimmed_fastq_files
end

"""
    generate_mirna_counts(input_sam_file::AbstractString)

This function creates a DataFrame containing all aligned miRNA and the number of times 
they appear.

# Returns

A DataFrame with the following columns:

* `name`: The name of the miRNA.
* `count`: The number of times the miRNA appears in the SAM file.
"""
function generate_mirna_counts(input_sam_file::AbstractString)
    miRNA_counts = Dict{String, Int64}()

    # Open the input sam file
    sam_file = open(SAM.Reader, input_sam_file)
    record = SAM.Record()

    # Go through each sequence in the sam file
    while !eof(sam_file)
        # Empty sam record to reduce memory usage
        empty!(record)
        read!(sam_file, record)

        # Unmapped reads have no ID
        if SAM.ismapped(record)
            miRNA_name = SAM.refname(record)

            # Check if the miRNA has already been added, update count if needed
            if haskey(miRNA_counts, miRNA_name)
                miRNA_counts[miRNA_name] += 1
            else
                miRNA_counts[miRNA_name] = 1
            end
        end
    end

    # Create a DataFrame with the miRNA names and counts
    miRNA_counts_dataframe = DataFrame(name = collect(keys(miRNA_counts))
                                    , count = collect(values(miRNA_counts))
                                    )
    # Sort the DataFrame by miRNA count in descending order
    sort!(miRNA_counts_dataframe, :count, rev = true)

    return miRNA_counts_dataframe
end

"""
    mirna_discovery_calculation(trimmed_fastqs::Vector{String}, sample_names::Vector{SubString{String}}, bowtie2_reference_name::String)

Align trimmed fastq files to a given human miRNA bowtie2 reference and calculate the number of unique miRNAs present in each sample.

The function returns:
1. A `Vector` of DataFrames where each DataFrame contains the counts of miRNAs for a specific sample.
2. A `Vector` of filenames of the SAM files generated during the alignment.

# Parameters
- `trimmed_fastqs`: A vector of filenames for the trimmed fastq files.
- `sample_names`: A vector of sample names corresponding to the fastq files.
- `bowtie2_reference_name`: The name of the bowtie2 reference against which the alignment should be done.
- `cores`: The number of threads to use for alignment.

# Returns
1. A `Vector` of DataFrames. Each DataFrame has two columns: `name` representing the miRNA name and `count` representing the number of times the miRNA was observed in the sample.
2. A `Vector` of filenames for the SAM files generated during the alignment.

# Example
```julia
julia> mirna_discovery_calculation(["sample1.cut.fastq", "sample2.cut.fastq", "sample3.cut.fastq"], 
                                   ["sample1", "sample2", "sample3"], 
                                   "human_miRNA_reference")
# Returns:
# ([469x2 DataFrame
#    Row │ name             count
#    │ String           Int64
#    ─────┼────────────────────────
#      1 │ hsa-miR-1307-3p    8548
#      2 │ hsa-miR-425-5p     3221
#      3 │ hsa-miR-488-5p     2951
#   ...], 
#  ["sample1.miRNA.sam", "sample2.miRNA.sam", "sample3.miRNA.sam"])
```
"""
function mirna_discovery_calculation(trimmed_fastqs::Vector{String}
                                    , sample_names::Vector{SubString{String}}
                                    , reference::AbstractString
                                    , cores::Int
                                    )
    number_of_records = length(trimmed_fastqs)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .5
                                            , "Calculating the number of miRNA present..."
                                            )

    vector_of_miRNA_counts_dfs = Vector{DataFrame}()

    for (fastq_file, sample_name) in zip(trimmed_fastqs, sample_names)
        #Align to miRBase reference
        wait(run(pipeline(
                `bowtie2
                --norc
                --threads $cores
                -x data/$reference
                -U $fastq_file
                -S $sample_name.miRNA.sam`
                , devnull)
                , wait = false
                ))

        counts_df = generate_mirna_counts(string(sample_name, ".miRNA.sam"))
        push!(vector_of_miRNA_counts_dfs, counts_df)
                
        next!(update_progress_bar)
    end

    sam_files::Vector{String} = capture_target_files(".miRNA.sam")

    return vector_of_miRNA_counts_dfs, sam_files
end

# Helper function to compute thresholds
function compute_thresholds(counts_column)
    threshold_of_one = sum(counts_column .>= 1)
    threshold_of_three = sum(counts_column .>= 3)
    threshold_of_five = sum(counts_column .>= 5)
    threshold_of_ten = sum(counts_column .>= 10)
    return [threshold_of_one, threshold_of_three, threshold_of_five, threshold_of_ten]
end

"""
    plot_mirna_counts(mirna_counts_dfs::Vector{DataFrame}, sample_names::Vector{SubString{String}})

Plot counts of the number of unique miRNA each sample aligned to. The function returns the names of the .csv files generated.

# Example
```julia
julia> plot_mirna_counts(["sample1.hairpin.hist","sample2.hairpin.hist","sample3.hairpin.hist"], 
                         ["sample1", "sample2", "sample3"])
["sample1_miRNA_counts.csv", "sample2_miRNA_counts.csv", "sample3_miRNA_counts.csv"]
```
"""
function plot_mirna_counts(mirna_counts_dfs::Vector{DataFrame}
                        , sample_names::Vector{SubString{String}}
                        )
    number_of_records = length(mirna_counts_dfs)
    update_progress_bar = progress_bar_update(number_of_records, .25, "Plotting miRNA counts from each sample...")

    # Create an empty vector to hold miRNA_counts file names
    miRNA_counts_files = Vector{String}()

    for (miRNA_df, sample_name) in zip(mirna_counts_dfs, sample_names)

        # Clean miRNA names
        miRNA_df[!, :name] = replace.(miRNA_df[!, :name], "(unk)" => "")

        # Calculate RPM and add as new column
        total_mapped_reads = sum(miRNA_df[!, :count])
        miRNA_df[!, :RPM] = round.(miRNA_df[!, :count] / total_mapped_reads * 10^6, digits=2)

        # Compute miRNA thresholds
        thresholds = compute_thresholds(miRNA_df[!, :count])
        
        # Barplot colors
        colors = [:grey88, :skyblue2, :peachpuff, :lightsalmon]

        # Set x and y-axis values; y-axis is a vector containing each sample's counted miRNAs
        x, y = 1:4, [thresholds[1], thresholds[2], thresholds[3], thresholds[4]]

        fig = Figure(resolution = (1800, 1200))
        increments = round(thresholds[1] / 10, sigdigits=2)
        ax = Axis(fig[1, 1]
                , xticks = (1:4
                , ["Threshold 1", "Threshold 3", "Threshold 5", "Threshold 10"])
                , yticks = 0:increments:thresholds[1]
                , ylabel = "Number of miRNA"
                , title = string(sample_name, " miRNA Counts")
        )
        CairoMakie.ylims!(ax, 0, thresholds[1] * 1.05)

        # Make sample barplot
        barplot!(ax
                , x
                , y
                , strokewidth = 1
                , color = colors
                , strokecolor = :black
                , bar_labels = :y
        )

        save(string(sample_name, "_miRNA_counts.pdf"), fig)
        push!(miRNA_counts_files, string(sample_name, "_miRNA_counts.csv"))
        CSV.write(miRNA_counts_files[end], miRNA_df)
        
        next!(update_progress_bar)
    end

    pdf_miRNA_files::Vector{String} = capture_target_files("_miRNA_counts.pdf")

    # Merge output PDFs into one file
    merge_pdfs([pdf_miRNA_files...], "Small_RNA_miRNA_plots.pdf")
    
    # Remove individual sample plots
    remove_files(pdf_miRNA_files)

    miRNA_counts_files::Vector{String} = capture_target_files("_miRNA_counts.csv")

    return miRNA_counts_files
end

"""
Use aligned SAM files to create a read length distribution file.
"""
function calculate_read_length_distribution(sam_files::Vector{String}
                                            , sample_names::Vector{SubString{String}}
                                            , cores::Int
                                            )
    number_of_records = length(sam_files)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .5
                                            , "Calculating Read Length Distribution..."
                                            )

    for (sam_file, sample_name) in zip(sam_files, sample_names)
	    #=
        If the SAM file is empty or almost empty, it will cause errors 
        downstream so the script exits.
        =#
        read_length_distibution = 
        try
            read(pipeline(
                `samtools 
                stats 
                -@ $cores
                $sam_file`
                , `grep ^RL`
                , `cut -f 2-`)
                , String
            )
        catch
            #=
            If the SAM file is empty or almost empty, it will cause errors 
            downstream so the script exits.
            =#
            println("Error, SAM file contains no read length distribution.")
            println("Exiting program...")
            exit()
        end

        next!(update_progress_bar)

        # Write read length data to output file
        output_length_file::IOStream = open("read_lengths_$sample_name.tsv", "w")
        write(output_length_file, string("Length", "\t", "Reads", "\n"))
        write(output_length_file, read_length_distibution)

        close(output_length_file)
    end

    length_files::Vector{String} = capture_target_files("read_lengths_")

    return length_files
end

"""
Create barplot of each sample's fragment lengths based on the number of reads found at each length
"""
function plot_fragment_lengths(length_files::Vector{String}
                                , sample_names::Vector{SubString{String}}
                                )
    number_of_records = length(length_files)
    update_progress_bar = progress_bar_update(number_of_records
                                            , 1
                                            , "Creating Read Length Plot..."
                                            )

    for (file, sample_name) in zip(length_files, sample_names)
        length_file::DataFrame = CSV.read(file, DataFrame)

        # Isolate fragment lengths and the number of reads at each length
        fragment_lengths::Vector{Int64} = length_file[!, :Length]
        reads_per_length::Vector{Int64} = length_file[!, :Reads]
        next!(update_progress_bar)

        plot = Figure(resolution = (1800, 1200))
        next!(update_progress_bar)

        Axis(plot[1, 1]
        , xticks = 15:5:maximum(fragment_lengths)
        , xlabel = "Fragment Lengths"
        , ylabel = "Number of Reads"
        , title = sample_name
        )
        next!(update_progress_bar)

        # Make sample barplot
        barplot!(fragment_lengths
        , reads_per_length
        , fillto = -1
        , color = fragment_lengths
        , strokecolor = :black
        , strokewidth = 1
        )

        save(string(sample_name, "_fragment_lengths.pdf"), plot)

        next!(update_progress_bar)
    end

    pdf_length_files::Vector{String} = capture_target_files("_fragment_lengths.pdf")

    merge_pdfs([pdf_length_files...], "Small_RNA_fragment_length_plots.pdf")

    remove_files(pdf_length_files)
end

"""
    align_with_bowtie2(fastq_file::String, read_count_dict::Dict{String, Int64})

Align fastq with specified reference RNA types.

Divide reads aligned with total reads to calculate percent alignment. Store 
alignment information in a dictionary.

# Example
```julia
julia> align_with_bowtie2("sample1.cut.fastq"
                        , Dict("sample1" => 6390314, "sample2" => 5000000, "sample3" => 7052928))
Dict{String, Float64} with 5 entries:
  "miRNA" => 65.0
  "tRNA" => 4.2
  "piRNA" => 1.2
  "snoRNA" => 1.9
  "rRNA" => 3.3
```
"""
function align_with_bowtie2(fastq_file::String
                            , read_count_dict::Dict{String, Int64}
                            , organism_name::String
                            , cores::Int
                            )
    sample_name = first(split(fastq_file, "."))

    # RNA types to be aligned to
    reference_RNA_types = 
    if organism_name != "hsa"
        ["miRNA"]
    else
        ["miRNA", "tRNA", "piRNA", "snoRNA", "rRNA"]
    end

    # Dictonary to hold each RNA types' alignment data
    metrics_dict = Dict{String, Number}()

    for reference_RNA in reference_RNA_types
        # Can skip aligning sample with miRNA reference since that was already done.
        sam_filename = ""
        if reference_RNA != "miRNA"
            sam_filename = string(sample_name, ".aligned.metrics.sam")
            # Align with bowtie2
            wait(run(pipeline(
                `bowtie2
                --threads $cores
                --norc
                -x data/$reference_RNA
                -U $fastq_file
                -S $sam_filename`
                , devnull)
                , wait = false
                )
            )
        else
            sam_filename = string(sample_name, ".miRNA.sam")
        end

        # Count number of read alignments
        aligned_reads_string = read(`
            samtools 
            view 
            -c 
            -F 4 
            $sam_filename`
            , String
        )
        aligned_reads = parse(Float64, aligned_reads_string)

        # Transfer unmapped reads to SAM file
        wait(run(pipeline(`
            samtools 
            view 
            -f 4 
            -o $sample_name.aligned.metrics.sam
            $sam_filename`
            , devnull)
            , wait = false
            )
        )

        # Convert SAM file back into a fastq
        wait(run(pipeline(`
            samtools 
            fastq 
            -@ 48 
            -n 
            $sample_name.aligned.metrics.sam`
            , stdout = fastq_file)
            , wait = false
            )
        )

        percent_alignment = round(100 * aligned_reads / read_count_dict[sample_name], digits = 3)
        metrics_dict[reference_RNA] = percent_alignment
    end

    return metrics_dict
end

"""
Calculate the bowtie2 alignment metrics for all samples. Return filename for plotting.
"""
function calculate_metrics(trimmed_fastqs::Vector{String}
                            , read_count_dict::Dict{String, Int64}
                            , sample_names::Vector{SubString{String}}
                            , dimer_count_dict::Dict{String, Float64}
                            , q_score_dict::Dict{String, Float64}
                            , organism_name::String
                            , cores::Int
                            )
    number_of_records = length(trimmed_fastqs)
    update_progress_bar = progress_bar_update(number_of_records, .5, "Calculating Sample Metrics...")
    output_metrics_file = open("Alignment_Metrics.csv", "w")

    # Mapping for organism specific metrics
    metrics_mapping = Dict(
        "default" => ["Read_Count", "Aligned_Reads", "Average_Q_score", "Dimer", "<16bp"
                    , "miRNA", "Unaligned_Reads"
        ],
        "hsa" => ["Read_Count", "Aligned_Reads", "Average_Q_score", "Dimer", "<16bp"
                , "miRNA", "tRNA", "piRNA", "snoRNA", "rRNA", "Unaligned_Reads"
        ]
    )

    # Determine which metrics to write based on organism_name
    metrics_to_write = get(metrics_mapping, organism_name, metrics_mapping["default"])

    # Write header
    write(output_metrics_file, join(["Sample"; metrics_to_write], ",") * "\n")

    for (fastq_file, sample_name) in zip(trimmed_fastqs, sample_names)
        # Calculate q_score score
        average_q_score = q_score_dict[sample_name]

        # Gather the number of canonical dimer reads
        percent_dimer = dimer_count_dict[sample_name]

        # Get number of short fragments from cutadapt output
        short_fragments_string = match(r"Reads that.+\(\K.+(?=%)"
                                        , readchomp(`cat $sample_name.cutadapt_information`)
                                        ).match
        short_fragments = abs(parse(Float64, short_fragments_string) - percent_dimer)

        # Function to align each fastq to each RNA type and calculate the number of reads mapped
        metrics_dict::Dict{String, Number} = align_with_bowtie2(fastq_file, read_count_dict, organism_name, cores)

        aligned_reads = round(values(metrics_dict) |> sum, digits = 2)
        unaligned_reads = round(100 - aligned_reads - percent_dimer - short_fragments, digits = 2)

        # Write data to output file based on determined metrics
        data_to_write = [
            sample_name, 
            read_count_dict[sample_name], 
            aligned_reads, 
            average_q_score, 
            round(percent_dimer, digits=2), 
            round(short_fragments, digits=2),
            metrics_dict["miRNA"],
            get(metrics_dict, "tRNA", ""),
            get(metrics_dict, "piRNA", ""),
            get(metrics_dict, "snoRNA", ""),
            get(metrics_dict, "rRNA", ""),
            unaligned_reads
        ]

        # Filter out empty strings from the data to write
        data_to_write = filter(!isempty, data_to_write)

        write(output_metrics_file, join(data_to_write, ",") * "\n")
        next!(update_progress_bar)
    end

    close(output_metrics_file)
    return "Alignment_Metrics.csv"
end

"""
Create violin plot of each sample's RNA aligment metrics
"""
function violin_plot_metrics(metrics_file::String)

    # Import sample metrics file
    metrics_file::DataFrame = CSV.read(metrics_file, DataFrame)

    # Create new dataframe without sample names and add sample and column names to arrays
    select!(metrics_file, Not([:Sample, Symbol("Read_Count")]))
    column_names::Vector{String} = names(metrics_file)
    num_of_cols = last(size(metrics_file))

    number_of_records = length(column_names)
    update_progress_bar = progress_bar_update(number_of_records, .25
    , "Creating Violin Plots..."
    )
    
    # Gather each metric's column for plotting
    metrics_columns(col) = [number for number in metrics_file[!, col]]

    # Make empty figure and set x and y axes
    # x-axis is vector of integers from 1 to the number of metrics
    # y-axis is vector of vectors containing the data from each metric
    fig = Figure(resolution = (1800, 1200))
    x, y = 1:num_of_cols, [metrics_columns(column) for column in 1:num_of_cols]
    
    # Create violin plot with data from all samples analyzed
    for column in 1:num_of_cols
        ax = Axis(fig[1, column]
        , yticks = 0:5:100
        , ylabel = "Percent"
        , title = column_names[column]
        )
        CairoMakie.ylims!(ax, 0, 100)

        # Make violin plot with combined sample data
        CairoMakie.violin!(fig[1, column]
        , repeat([x[column]]
        , first(size(metrics_file)))
        , y[column]
        , show_median=true
        )
        next!(update_progress_bar)

        # Save file once all columns have been added
        if column == num_of_cols
            save("Violin_plot_metrics.png", fig)
        end

        next!(update_progress_bar)
    end
end

"""
Create barplot of each sample's RNA aligment metrics.
"""
function plot_metrics(metrics_file::String
                        , library::String
                        , sample_names::Vector{SubString{String}}
                        )

    metrics_file::DataFrame = CSV.read(metrics_file, DataFrame)

    select!(metrics_file, Not([:Sample, Symbol("Read_Count")]))
    column_names::Vector{String} = names(metrics_file)
    num_of_cols = last(size(metrics_file))

    colors = [:snow3, :honeydew2, :lightpink, :bisque2, :deepskyblue, :aquamarine3
            , :peachpuff, :paleturquoise1, :orange1, :lightsalmon2
    ]
    
    sample_rows(row) = [values for values in metrics_file[row, :]]

    number_of_records = length(sample_names)
    update_progress_bar = progress_bar_update(number_of_records
                                            , .25
                                            , "Creating Metrics Plots..."
                                            )
    
    for sample in eachindex(sample_names)
        sample_name = sample_names[sample]

        # Set axis values; y-axis is a vector containing each sample's calculated metrics
        x, y = 1:num_of_cols, sample_rows(sample)
        fig = Figure(resolution = (1800, 1200))

        ax = Axis(fig[1, 1]
            , xticks = (1:num_of_cols
            , column_names)
            , yticks = 0:10:100
            , ylabel = "Percent"
            , title = string(sample_name, " Metrics - ", library)
        )

        # Set y-axis limits
        CairoMakie.ylims!(ax, 0, 105)

        barplot!(ax
            , x
            , y
            , strokewidth = 1
            , strokecolor = :black
            , color = colors
            , bar_labels = :y
        )

        next!(update_progress_bar)

        save(string(sample_name, "_metrics.pdf"), fig)
    end

    pdf_metrics_files::Vector{String} = capture_target_files("_metrics.pdf")

    merge_pdfs([pdf_metrics_files...], "Small_RNA_metrics_plots.pdf")
    
    remove_files(pdf_metrics_files)
end

"""
    find_common_mirnas(mirna_counts_files::Vector{String}
                            , read_count_dict::Dict{String, Int64}
                            , sample_names::Vector{SubString{String}}
                            )

Return the miRNAs present in all samples, together with their counts and reads per million reads (RPM).

# Arguments
- `mirna_counts_files`: A vector of files containing sample miRNA counts.
- `read_count_dict`: A dictionary where keys are sample names and values are the total number of reads 
for the corresponding sample.
- `sample_names`: A vector of sample names.

# Returns
- A Tuple of three items:
    1. A list of miRNAs that are common to all samples.
    2. A list of dictionaries where each dictionary contains the miRNA counts for one sample.
    3. A list of dictionaries where each dictionary contains the RPM values for the miRNAs in one sample.

# Example
```julia
julia> find_common_mirnas(["sample1_miRNA_counts.csv", "sample2_miRNA_counts.csv"]...,
                        , Dict("sample1" => 6390314, "sample2" => 5000000, "sample3" => 7052928)
                        , ["sample1", "sample2", "sample3"])
```
"""
function find_common_mirnas(mirna_counts_files::Vector{String}
                            , read_count_dict::Dict{String, Int64}
                            , sample_names::Vector{SubString{String}}
                            )

    # Vector of dictionaries containing the miRNA counts from each sample
    miRNA_info = Vector{Dict{String, Int64}}()
    RPM_info = Vector{Dict{String, Float64}}()
    miRNA_names_dict = Dict{String, Int64}()

    for (file, sample_name) in zip(mirna_counts_files, sample_names)

        # Read all lines from the file at once
        miRNA_file = readlines(file)

        sample_miRNA_counts_dictionary = Dict{String, Int64}()
        sample_miRNA_RPM_dictionary = Dict{String, Float64}()

        # Process each line in the file (skipping the first line, which is a header)
        for line in miRNA_file[2:end]
            split_line = split(line, ",")
            miRNA_count = parse(Int64, split_line[2])
            miRNA_name = first(split_line)

            # Calculate RPM for the current miRNA
            RPM = round(miRNA_count / (read_count_dict[sample_name] / 10^6), digits = 4)

            # Store miRNA count and RPM in the dictionaries for the current sample
            sample_miRNA_counts_dictionary[miRNA_name] = miRNA_count
            sample_miRNA_RPM_dictionary[miRNA_name] = RPM

            # Increment the count for the current miRNA name in the global dictionary
            miRNA_names_dict[miRNA_name] = get(miRNA_names_dict, miRNA_name, 0) + 1
        end

        # Add the dictionaries for the current sample to the storage vectors
        push!(miRNA_info, sample_miRNA_counts_dictionary)
        push!(RPM_info, sample_miRNA_RPM_dictionary)
    end

    # Find miRNA present in all samples, i.e., have counts equal to the number of samples analyzed
    miRNAs_in_common = filter(miRNA -> last(miRNA) === length(sample_names), miRNA_names_dict) |> keys

    return miRNAs_in_common, miRNA_info, RPM_info
end

"""
Write miRNAs all samples have in common to output file.
"""
function write_common_mirna_file(miRNA_Names::Base.KeySet{String, Dict{String, Int64}}
                                , mirna_info::Vector{Dict{String, Int64}}
                                , rpm_info::Vector{Dict{String, Float64}}
                                , sample_names::Vector{SubString{String}}
                                )

    for (index, sample) in enumerate(sample_names)
        common_miRNA_file::IOStream = open(string(sample, "_common_miRNAs.tsv"), "w")
        common_miRNA_file_RPM::IOStream = open(string(sample, "_common_miRNAs_RPM.tsv"), "w")

        if index == 1
            write(common_miRNA_file, string("miRNA", "\t", sample, "\n"))
            write(common_miRNA_file_RPM, string("miRNA", "\t", sample, "\n"))
            for miRNA in miRNA_Names
                write(common_miRNA_file, string(miRNA, "\t", mirna_info[index][miRNA], "\n"))
                write(common_miRNA_file_RPM, string(miRNA, "\t", rpm_info[index][miRNA], "\n"))
            end
        else
            write(common_miRNA_file, string(sample, "\n"))
            write(common_miRNA_file_RPM, string(sample, "\n"))
            for miRNA in miRNA_Names
                write(common_miRNA_file, string(mirna_info[index][miRNA], "\n"))
                write(common_miRNA_file_RPM, string(rpm_info[index][miRNA], "\n"))
            end
        end

        close(common_miRNA_file)
        close(common_miRNA_file_RPM)
    end

    common_files::Vector{String} = capture_target_files("_common_miRNAs.tsv")
    common_RPM_files::Vector{String} = capture_target_files("_common_miRNAs_RPM.tsv")

    # Combine all the separate common miRNA files into one file
    run(pipeline(`paste $common_files`, stdout = "Common_miRNAs.tsv"))
    run(pipeline(`paste $common_RPM_files`, stdout = "Common_miRNAs_RPM.tsv"))

    # Run principal component analysis and UMAP on common miRNA RPM
    plot_clustering("Common_miRNAs_RPM.tsv")

    remove_files(common_files)
    remove_files(common_RPM_files)
end

"""
    plot_clustering(Common_miRNA_File::String)

Generates Principal Component Analysis (PCA) and Uniform Manifold Approximation and 
Projection (UMAP) plots for libraries given a common miRNA RPM file.

This function performs the following steps:
1. Reads the common miRNA RPM file and creates a DataFrame.
2. Extracts sample names and miRNA names.
3. Ensures that there are at least two samples and more than one miRNA present in the dataset.
4. Transforms the data using PCA and prepares it for plotting.
5. If there are at least two principal components, plots the PCA and saves it as "common_miRNA_PCA.png".
6. If there are at least three principal components, performs UMAP dimensionality reduction and clustering using K-medoids.
7. Plots the UMAP and saves it as "common_miRNA_UMAP.png".
8. Writes the PCA and UMAP information to separate CSV files for cluster tracking.

# Arguments
- `Common_miRNA_File::String`: Path to the common miRNA RPM file.

# Outputs
- Create and saves PCA and UMAP plots.
- Write PCA and UMAP information.
"""
function plot_clustering(Common_miRNA_File::String)
    common_miRNA_counts = DataFrame(CSV.File(Common_miRNA_File))
    sample_names = DataFrames.names(common_miRNA_counts, Not([:miRNA]))
    miRNA_names = common_miRNA_counts[!, :miRNA]

    # There must be at least two samples to perform the principal component analysis
    if last(size(common_miRNA_counts)) > 2 && first(size(common_miRNA_counts)) > 1
        common_miRNA_counts_matrix = Matrix{Float64}(select(common_miRNA_counts, Not([:miRNA])))

        # PCA covariance matrix for PCA plot
        pca_matrix = fit(PCA, common_miRNA_counts_matrix; maxoutdim=20)
        transformed_counts = MultivariateStats.transform(pca_matrix, common_miRNA_counts_matrix)

        # Transposed PCA covariance matrix for UMAP plot
        umap_pca_matrix = fit(PCA, common_miRNA_counts_matrix'; maxoutdim=20)
        transposed_transformed_counts = MultivariateStats.transform(umap_pca_matrix, common_miRNA_counts_matrix')

        # Plot PCA if there are at least two principal components
        if first(size(transformed_counts)) >= 2
            sample_name_df = DataFrame("samples" => sample_names)

            # Write PCA information to file for cluster tracking
            pca_values = hcat(sample_name_df, DataFrame(transformed_counts[1:2, :]', ["PC1", "PC2"]))
            CSV.write("PCA_information.csv", pca_values)

            # Plot PCA
            StatsPlots.scatter(size=(1200, 800), dpi=300, titlefont=(16, "Computer Modern"),
                    xlabel="PC1", ylabel="PC2", title="Common miRNA: PCA",
                    transformed_counts[1, :], transformed_counts[2, :],
                    left_margin=23mm, right_margin=8mm, bottom_margin=8mm,
                    leg=false
            )
            savefig("common_miRNA_PCA.png")
        end

        # Plot UMAP if there are at least three principal components
        if first(size(transposed_transformed_counts)) >= 3
            miRNA_df = DataFrame("miRNA" => miRNA_names)

            # Create low-dimensional embedding for UMAP
            near_neighbors = first(size(common_miRNA_counts_matrix)) - 1
            embedding = umap(transposed_transformed_counts, 2;
                             n_neighbors=min(near_neighbors, 10)
                             , min_dist=0.01
            )

            # Create distance matrix for K-means clustering
            distance_matrix = pairwise(Euclidean(), embedding, embedding)
            kmeans_clusters = kmeans(distance_matrix, first(size(transposed_transformed_counts)))

            # Make dataframe with miRNA IDs, UMAP values, and cluster IDs
            DataFrames.hcat!(miRNA_df
                            , DataFrame(embedding', ["UMAP1", "UMAP2"])
                            , DataFrame("Cluster_ID" => kmeans_clusters.assignments)
            )
            # Sort DataFrame by cluster IDs
            sort!(miRNA_df, :Cluster_ID)

            # Write UMAP information to file for cluster tracking
            CSV.write("UMAP_information.csv", miRNA_df)

            # Plot UMAP of common miRNAs; color UMAP clusters based on K-means clustering
            StatsPlots.scatter(size=(1200, 800), embedding[1, :], embedding[2, :],
                    title="Common miRNA: UMAP", left_margin=13mm,
                    bottom_margin=10mm, dpi=300,
                    marker_z=kmeans_clusters.assignments,
                    color=:lighttest, xlabel="UMAP1", ylabel="UMAP2",
                    titlefont=(16, "Computer Modern"), leg=false,
                    markersize=9, markerstrokewidth = 0.1
            )
            savefig("common_miRNA_UMAP.png")
        end
    end
end

"""
    cq_vs_mirna_count_correlation(qPCR_Data_File::String,
                                  Common_miRNA_RPM_File::String,
                                  Common_miRNA_File::String,
                                  sample_names::Vector{SubString{String}}) -> DataFrame

Calculate the Pearson's correlation coefficient between qPCR data and common miRNA counts.

Functionality:
- Filters the qPCR data to exclude measurements without a quantification cycle (Cq) value.
- Intersects the filtered qPCR data with the miRNAs from the common miRNA files.
- Computes the average Cq value for each miRNA across replicates.
- Determines the Pearson's correlation coefficient between the qPCR data and the common miRNA counts.
- If the Pearson's coefficient is < -0.6 for any sample, a linear regression is plotted for that sample.

## Parameters:
- `qPCR_Data_File`: Path to the file containing qPCR data.
- `Common_miRNA_RPM_File`: Path to the file containing RPM data for common miRNAs.
- `Common_miRNA_File`: Path to the file listing common miRNAs.
- `sample_names`: List of sample names to process.

## Returns:
- A DataFrame containing miRNA names, their average Cq values, and RPM counts for each sample.

## Example:
```julia
julia> cq_vs_mirna_count_correlation("qpcr_raw_data.csv", "Common_RPM_miRNAs.csv",
                                     "Common_miRNAs.csv", ["sample1", "sample2", "sample3"])                   
# Outputs:
# A DataFrame with columns for miRNA, Cq, and RPM values for each sample.
# A CSV file named "Sample_correlation_coefficients.csv" containing sample names and their correlation coefficients.
```
"""
function cq_vs_mirna_count_correlation(qPCR_Data_File::String,
                                    Common_miRNA_RPM_File::String,
                                    Common_miRNA_File::String,
                                    sample_names::Vector{SubString{String}}
                                    )
    update_progress_bar = progress_bar_update(4
                                            , .25
                                            , "Calculating qPCR data vs miRNA count correlation data..."
                                            )

    # Import qPCR data file and both common miRNA files
    qPCR_data_table_full = DataFrame(CSV.File(qPCR_Data_File))
    common_miRNA_counts = DataFrame(CSV.File(Common_miRNA_File))
    RPM_dataframe = DataFrame(CSV.File(Common_miRNA_RPM_File))

    # Ensure there are at least two miRNAs for calculating correlation data
    if first(size(common_miRNA_counts)) > 1
        next!(update_progress_bar)

        # Filter out measurements with no Cq value
        qPCR_data_table = filter(row -> row.Cq != "NA", qPCR_data_table_full)

        # Get the intersection between all common miRNA and miRNA list after filtering
        common_miRNA = intersect(common_miRNA_counts[!, :miRNA], qPCR_data_table[!, :miRname])

        # Filter the RPM dataframe to include only common miRNAs
        filtered_RPM_dataframe = filter(row -> row.miRNA in common_miRNA, RPM_dataframe)
        next!(update_progress_bar)

        # Calculate the average Cq values for all common miRNAs
        average_cqs_values = [
            mean(parse(Float64, cq) for cq in qPCR_data_table[qPCR_data_table.miRname .== miRNA, :Cq])
            for miRNA in common_miRNA
        ]

        next!(update_progress_bar)

        # Create a new dataframe with the names of all common miRNAs and their Cq values
        cq_miRNA_dataframe = DataFrame(miRNA=common_miRNA, Cq=average_cqs_values)

        # Combine the Cq values and miRNA names with the RPM data
        combined_cq_RPM_dataframe = hcat(cq_miRNA_dataframe, filtered_RPM_dataframe[:, 2:end])

        # Calculate Pearson's correlation coefficients between the qPCR values and the log10 of miRNA counts
        correlation_coeffs = [
            cor(combined_cq_RPM_dataframe[:, :Cq], log10.(combined_cq_RPM_dataframe[:, sample]))
            for sample in sample_names
        ]
        
        # Write correlation values to output file
        correlation_dataframe = DataFrame(hcat(sample_names, correlation_coeffs), [:sample, :corr_coeff])
        CSV.write("Sample_correlation_coefficients.csv", correlation_dataframe)

        # Find samples with correlation values < -.6
        negative_correlation_samples = sample_names[correlation_coeffs .< -0.6]

        next!(update_progress_bar)

        plot_cq_linear_regression(negative_correlation_samples, combined_cq_RPM_dataframe)
    end
end

"""
Make linear regression plot for samples correlated with the qPCR data.
"""
function plot_cq_linear_regression(Negative_Correlation_Samples::Vector{SubString{String}}
                                    , Combined_Cq_RPM_Dataframe::DataFrame
                                    )
    #=
    Loop through samples with a Pearson correlation coefficient of at least -.5 and create 
    linear regression plot Since Cq and miRNA count are inversely correlated, we want the 
    samples with a negative correlation coefficient.
    =#
    for sample in Negative_Correlation_Samples
        x_values::Vector{Float64} = map(log10, Combined_Cq_RPM_Dataframe[!, sample])
        y_values::Vector{Float64} = Combined_Cq_RPM_Dataframe[!, :2]
    
        # Make new dataframe with Cq values and log10 converted miRNA counts
        sample_cq_dataframe::DataFrame = hcat(DataFrame(name = x_values), DataFrame(Cq = y_values))
    
        # Create linear regression model
        linear_model = lm(@formula(Cq ~ name), sample_cq_dataframe)
    
        # Use linear model to predict the Cq values expected based on the miRNA counts
        predictions = predict(linear_model, sample_cq_dataframe, interval = :confidence, level = 0.95)
    
        # Create scatter plot with Cq vs the log of miRNA read counts
        linear_plot = @df sample_cq_dataframe Plots.scatter(:name, :Cq, leg = false, markersize = 7)
    
        # Sort data on log10 miRNA counts; reverse order since the values are negatively correlated
        sorted_predictions::DataFrame = predictions[sortperm(predictions[!, :prediction], rev = true), : ]
        sorted_x_values::Vector{Float64} = sort(x_values)
    
        # Extends plot y-axis by 1/7 each side
        dY::Float64 = first(diff([extrema(y_values)...])) / 7  
        y1, y2 = (-dY, dY) .+ extrema(y_values)
    
        # Goodness of fit
        M::Matrix{Float64} = [ones(length(x_values)) x_values]
        a, b = M \ y_values
        Ŝ = round(sqrt(sum((y_values - b*x_values .- a).^2) / (length(x_values) - 2)), digits = 3)

        # Make plot Cq vs log10 miRNA read counts with R-squared, Goodness of fit, and Pearson's coefficient in title
        Plots.plot!(size = (1200, 800), title = string(sample, " Linear Regression", "\n", " R-Squared = "
            , round(r2(linear_model), digits = 3), "\n", "Goodness of Fit = ", Ŝ, "\n", "Pearson \U03c1 = "
            , round(cor(x_values, y_values), digits = 3)), xlabel = "log10 miRNA Read Count"
            , ylabel = "Average Cq", dpi = 300, titlefont = (14, "Computer Modern")
            , ylims = (y1, y2), guidefont = (12, :black), tickfont = (11, :black), left_margin = 13mm, bottom_margin = 10mm
            , linear_plot, sorted_x_values, sorted_predictions.prediction, linewidth = 4
            , ribbon = (sorted_predictions.prediction .- sorted_predictions.lower
            , sorted_predictions.upper .- sorted_predictions.prediction
            )
        )
    
        savefig(string(sample, "_miRNA_correlation.png"))
    end
end

"""
Spot remove unecessary intermediate files.
"""
function remove_files(files_to_delete::Vector{String})
    for file in files_to_delete
        rm(file)
    end
end

function remove_files(file_to_delete::String)
    rm(file_to_delete)
end

"""
Remove all intermediate files.
"""
function remove_intermediate_files()
    files_to_delete = Set(vcat(
    capture_target_files(".bam")
    ,capture_target_files(".cut.fastq")
    ,capture_target_files("read_lengths_")
    ,capture_target_files(".miRNA.")
    ,capture_target_files(".cutadapt_information")
    ,capture_target_files(".sam")
    ))

    for file in files_to_delete
        rm(file)
    end
end

function parse_commandline()
    arguments = ArgParseSettings(prog="Small RNA Analysis"
                                , description = "Basic analysis of NEXTFLEX Small RNA libraries."
    )

    @add_arg_table! arguments begin
        "--library", "-l"
            help = "The type of sample library, either 'v4' or 'v3'."
            arg_type = String
            default = "v4"
        "--fasta", "-f"
            help = "The full path to a miRNA fasta file, e.g., /home/user/mirna.fa."
            arg_type = String
            default = "data/mirgene.fas"
		"--organism", "-O"
            help = "Abbreviated name of organism, e.g., 'hsa' for human or 'mmu' for mouse. \
            This should match the standard three letter abbreviation found in miRNA databases \
            such as miRBase and MirGeneDB."
            arg_type = String
            default = "hsa"
        "--threads", "-p"
            help = "The number of processors to use for alignment."
            arg_type = Int
            default = 12
    end

    args = parse_args(arguments)

    # Create a Config instance from parsed arguments
    config = Config(
        get(args, "library", "v4"),
        get(args, "fasta", "data/mirgene.fas"),
        get(args, "organism", "hsa"),
        get(args, "threads", 12)
    )

    return config
end

function julia_main()::Cint

    # Parse command line arguments from Config struct
    config = parse_commandline()
    
    library = config.library

    # Make miRNA bowtie2 reference if necessary
    create_target_organism_fasta_file(config.organism, config.fasta)

    #Say hello Trish!
    run(`echo " "`)
    run(`echo -e    "\e[1;31m\t..&&&&&....&&&&&&&&&...&&&&......&&&.........&"`)
    run(`echo -e   "\t.&&&&&&&...&&&&&&&&&...&&&&&.....&&&........&&&"`)
    run(`echo -e   "\t&&&...&&&..&&&...&&&...&&&.&&....&&&.......&&.&&"`)
    run(`echo -e   "\t.&&&&......&&&&&&&&&...&&&..&&...&&&......&&...&&"`)
    run(`echo -e   "\t...&&&&....&&&&&&&&&...&&&...&&..&&&.....&&&&&&&&&"`)
    run(`echo -e   "\t......&&&..&&&...&&....&&&....&&.&&&....&&.......&&"`)
    run(`echo -e   "\t&&..&&&&...&&&...&&&...&&&.....&&&&&...&&.........&&"`)
    run(`echo -e   "\t.&&&&&&....&&&...&&&&..&&&......&&&&..&&...........&&"`)
    run(`echo -e "               (                           )"`)
    run(`echo -e "               (                           )"`)
    run(`echo -e "          ) )( (                           ( ) )( ("`)
    run(`echo -e "       ( ( ( )  ) )                     ( (   (  ) )("`)
    run(`echo -e "      ) )     ,,\\\\\                     ///,,       ) ("`)
    run(`echo -e "   (  ((    (\\\\\\\//                     \\\////)      )"`)
    run(`echo -e "    ) )    (-(..//                       \\\..)-)     ("`)
    run(`echo -e "   (((   ((-(..||                         ||..)-))    ) )"`)
    run(`echo -e "  ) )   ((-(-(.||           '''\..        ||.)-)-))   (("`)
    run(`echo -e "  ((   ((-(-(/(/\\        ''; \e[1;33mo\e[0m\e[1;31m.- '      //\)\)-)-))    )"`)
    run(`echo -e "   )   (-(-(/(/(/\\      '';;;;-\~      //\)\)\)-)-)   (   )"`)
    run(`echo -e "(  (   ((-(-(/(/(/\======,:;:;:;:,======/\)\)\)-)-))   )"`)
    run(`echo -e "    )  '(((-(/(/(/(//////:\e[1;30m%%%%%%%:\e[0m\e[1;31m\\\\\\\\\)\)\)\)-))-)'  ( ("`)
    run(`echo -e "   ((   '((-(/(/(/('uuuu:\e[1;30mWWWWWWWWW\e[0m\e[1;31m:uuuu')\)\)\)-))'    )"`)
    run(`echo -e "     ))  '((-(/(/(/('|||:\e[1;30mwwwwwwwww\e[0m\e[1;31m:|||')\)\)\)-))'    (("`)
    run(`echo -e "  (   ((   '((((/(/('uuu:\e[1;30mWWWWWWWWW\e[0m\e[1;31m:uuu')\)\))))'     ))"`)
    run(`echo -e "        ))   '':::UUUUUU:\e[1;30mwwwwwwwww\e[0m\e[1;31m:UUUUUU:::''     ((   )"`)
    run(`echo -e "          ((      '''''''\uuuuuuuu/'''''''         ))"`)
    run(`echo -e "           ))            'JJJJJJJJJ'           (("`)
    run(`echo -e "             ((           LLLLLLLLLLL         ))"`)
    run(`echo -e "               ))        ///|||||||\\\       (("`)
    run(`echo -e "                 ))     (/(/(/(^)\)\)\)       (("`)
    run(`echo -e "                  ((                           ))"`)
    run(`echo -e "                    ((                       (("`)
    run(`echo -e "                      ( )( ))( ( ( ) )( ) (()"`)
    run(`echo -e "\t\t     \e[1;31m\e[1;4mP  i  p  e  l  i  n  e\e[0m"`)
    run(`echo " "`)
    run(`echo -e "\t\t\e[0;31mBut\e[0m \e[1;31myou\e[0m \e[1;33mcan\e[0m \e[1;32mcall\e[0m \e[0;36mme\e[0m \e[1;35mTrish...\e[0m "`)
    run(`echo " "`)

    fastqs = capture_target_files("_R1_001.fastq.gz")
    sample_names = map(sample -> first(split(sample, "_")), fastqs)
    read_count_dict, dimer_count_dict, q_score_dict = parse_fastqs(fastqs
                                                                    , sample_names
                                                                    , config.library
                                                                    )
    trimmed_fastqs = trim_adapters(fastqs, config.library, sample_names)
    mirna_counts_dfs, sam_files = mirna_discovery_calculation(trimmed_fastqs
                                                            , sample_names
                                                            , config.organism
                                                            , config.threads
                                                            )
    length_files = calculate_read_length_distribution(sam_files, sample_names, config.threads)
    plot_fragment_lengths(length_files, sample_names)
    mirna_counts_files = plot_mirna_counts(mirna_counts_dfs, sample_names)
    metrics_file = calculate_metrics(trimmed_fastqs
                                    , read_count_dict
                                    , sample_names
                                    , dimer_count_dict
                                    , q_score_dict
                                    , config.organism
                                    , config.threads
                                    )
    plot_metrics(metrics_file, library, sample_names)
    violin_plot_metrics(metrics_file)
    full_mirna_names_list, mirna_info, rpm_info = find_common_mirnas(mirna_counts_files
                                                                    , read_count_dict
                                                                    , sample_names
                                                                    )
    write_common_mirna_file(full_mirna_names_list, mirna_info, rpm_info, sample_names)
    if isfile("qpcr_raw_data.csv")
        cq_vs_mirna_count_correlation("qpcr_raw_data.csv"
                                    , "Common_miRNAs_RPM.tsv"
                                    , "Common_miRNAs.tsv"
                                    , sample_names
                                    )
    end
    remove_intermediate_files()
    println(" ")
    println("Analysis Finished")
    println(" ")

    return 0
end

end
