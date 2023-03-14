module small_RNA_analysis

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

"""
    CaptureTargetFiles(Files_To_Capture::String)

List all files in the current directory. 

Check to see if each file contains the target file's string. 

#Example
```julia
julia> CaptureTargetFiles(".txt")
3-element Vector{String}:
 "file1.txt"
 "file2.txt"
 "file3.txt"
```
"""
function CaptureTargetFiles(Files_To_Capture::String)
    return [file for file in readdir() if occursin(Files_To_Capture, file)]
end

"""
Update progress bar on the command line.
"""
function ProgressBarUpdate(Number_of_Records::Int64, Interval::Number
                            , Description::String
                            )
    progress_bar_update = Progress(Number_of_Records
    , dt=Interval
    , barglyphs=BarGlyphs("[=> ]")
    , barlen=100
    , desc=Description
    )

    return progress_bar_update
end

"""
    GetReadQScore!(Line::String, Q_Score_List::Vector{Number})

Calculate the quality score for a given quality string.

Convert quality string to Phred33 score and update q score vector.

#Example
```julia
julia> GetReadQScore!("FGFEEE<FC<FGGGGGGGGGGGFFGFG8<@8CFFGF8EFGGCGFGGGGGGFG", [36.2, 35.9])
3-element Vector{Float64}:
 36.2
 35.9
 35.7
```
"""
function GetReadQScore!(Line::String, Q_Score_List::Vector{Number})

    #Find read quality by converting quality encodings to Unicode numbers
    q_score = sum(codeunits(Line)) / length(Line)

    #Add Q score to list and convert to a Phred33 score
    push!(Q_Score_List, q_score - 33) 

    return Q_Score_List
end

"""
    DetermineLibraryType(Fastqs::Vector{String}
                        , Sample_Names::Vector{SubString{String}}
                        )

Check each library to determine what kind of library it is.

Distinguish between library types by looking at the number of reads matching the expected 
miRNA structure.

Calculate the number of reads in each fastq file.

Create dictionary with each sample name and its read count.

Calculate the percent dimer present in a given fastq.

Calculate the average Q score for a given fastq. Divide sum of Unicode converted quality 
strings by the number of quality strings.

#Example
```julia
julia> DetermineLibraryType(["sample1.fastq.gz","sample2.fastq.gz","sample3.fastq.gz"]
                            , ["sample1", "sample2", "sample3"])
Dict{String, String} with 3 entries:
  "sample1" => "v4"
  "sample2" => "v4"
  "sample3" => "v3"
Dict{String, String} with 3 entries:
  "sample1" => "6390314"
  "sample2" => "5000000"
  "sample3" => "7052928"
```
"""
function DetermineLibraryType(Fastqs::Vector{String}
                                , Sample_Names::Vector{SubString{String}}
                                )

    number_of_records = length(Fastqs)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .5
                                            , "Determining sample library type..."
                                            )

    library_type = Dict{String, String}()
    read_count_dict = Dict{String, Int64}()
    dimer_count_dict = Dict{String, Float64}()
    q_score_dict = Dict{String, Number}()
    sample_tracker = 1

    for fastq_file in Fastqs
        sample_name = Sample_Names[sample_tracker]
        #List to store each reads quality information
        q_score_list = Vector{Number}()
        #List to store each read's information
        unique_start_sequences = Set{SubString{String}}()
        v3_miRNA_library_count = 0
        #Store dimer information
        v3_dimer_count = 0
        v4_dimer_count = 0

        fastq = GZip.open(fastq_file)
        sample_read_count = 1
        line_tracker = 1

        for line in eachline(fastq)
            #Ignore lines without sequence or quality information
            if line_tracker % 4 != 2 && line_tracker % 4 != 0
                line_tracker += 1
            elseif line_tracker % 4 == 2
                #=
                Identify if read contains a "standard" v3 miRNA library, .i.e. 4Ns bookending 
                a 22 bp miRNA. Counts of these are used to determine whether a library is 
                from small RNA v3 or v4. Checks for presence of canonical 0 bp dimer too.
                =#
                v3_library = match(r"^.{4}.{19,24}.{4}(?=TGGAATTCTCGGGTGCCAAGG)", line)
                v3_dimer = startswith(line[9:end], "TGGAATTCTCGGGTGCCAAGG")
                v4_dimer = startswith(line, "TGGAATTCTCGGGTGCCAAGG")

                if !isnothing(v3_library)
                    v3_miRNA_library_count += 1
                end

                if v3_dimer
                    v3_dimer_count += 1
                elseif v4_dimer
                    v4_dimer_count += 1
                end

                sample_read_count += 1
                line_tracker += 1
            elseif line_tracker % 4 == 0
                #Calculate average quality score
                GetReadQScore!(line, q_score_list)
                line_tracker += 1
            end  
        end

        close(fastq)

        #Add read count to dictionary for use later
        read_count_dict[sample_name] = sample_read_count

        #Calculate q score for each sample
        average_q_score = round(mean(q_score_list), digits = 2)
        q_score_dict[sample_name] = average_q_score

        percent_v3_miRNA_library = 100 * v3_miRNA_library_count / sample_read_count
        percent_v3_dimer = 100 * v3_dimer_count / sample_read_count
        percent_v4_dimer = 100 * v4_dimer_count / sample_read_count

        #Percentage cutoff is somewhat arbitrary, though should hold for most cases
        if percent_v3_miRNA_library > 5
            library_type[sample_name] = "v3"
            dimer_count_dict[sample_name] = percent_v3_dimer
        else
            library_type[sample_name] = "v4"
            dimer_count_dict[sample_name] = percent_v4_dimer
        end

        sample_tracker += 1
        next!(progress_bar_update)
    end

    return library_type, read_count_dict, dimer_count_dict, q_score_dict
end

"""
    TrimAdapters(Fastqs::Vector{String}
                    , Library_Type::Dict{String, String}
                    , Sample_Names::Vector{SubString{String}}
                    )

Trim the 3' adapter from each read.

V3 Read 1 Setup:\n
               5' Adapter            -        Insert          -     3' Adapter      
GATCGTCGGACTGTAGAACTCTGAACNNNN - TGTCAGTTTGTCAAATACCCCA - NNNNTGGAATTCTCGGGTGCCAAGG 

V4 Read 1 Setup:\n
               5' Adapter            -        Insert          -     3' Adapter      
AGATCGGAAGAGCGTCGTGTAGGGAAAGA - TGTCAGTTTGTCAAATACCCCA - TGGAATTCTCGGGTGCCAAGG 

The bases on the 3' end are also quality trimmed if their quality score is below 20. Reads 
shorter than 16 bases or that weren't trimmed are discarded.

#Example
```julia
julia> TrimAdapters(["sub_sample1.fastq","sub_sample2.fastq","sub_sample3.fastq"]
                        , ["sample1" => "v4","sample2" => "v4","sample3" => "v3"]
                        , ["sample1", "sample2", "sample3"])
3-element Vector{String}:
 "sample1.cut.fastq"
 "sample2.cut.fastq"
 "sample3.cut.fastq"
```
"""
function TrimAdapters(Fastqs::Vector{String}, Library_Type::Dict{String, String}
                        , Sample_Names::Vector{SubString{String}}
                        )
    number_of_records = length(Fastqs)
    progress_bar_update = ProgressBarUpdate(number_of_records, .5, "Trimming adapters...")

    sample_tracker = 1

    for fastq_file in Fastqs
        sample_name = Sample_Names[sample_tracker]

        if Library_Type[sample_name] == "v4"
            min_length = 16

            wait(run(pipeline(
            `cutadapt 
            --cores=0 
            --discard-untrimmed 
            --quality-cutoff 20 
            --adapter TGGAATTCTCGGGTGCCAAGG 
            --output $sample_name.cut.fastq 
            --minimum-length $min_length
            --too-short-output $sample_name.short.fastq
            $fastq_file`
            , stdout = "$sample_name.cutadapt_information")
            , wait = false
            ))
        elseif Library_Type[sample_name] == "v3"
            min_length = 16

            wait(run(pipeline(
            `cutadapt 
            --cores=0 
            --cut 4
            --discard-untrimmed 
            --quality-cutoff 20 
            --adapter N"{"4"}"TGGAATTCTCGGGTGCCAAGG 
            --output $sample_name.cut.fastq 
            --minimum-length $min_length
            --too-short-output $sample_name.short.fastq
            $fastq_file`
            , stdout = "$sample_name.cutadapt_information")
            , wait = false
            ))
        end

        next!(progress_bar_update)
        sample_tracker += 1
    end

    trimmed_fastq_files::Vector{String} = CaptureTargetFiles(".cut.fastq")

    return trimmed_fastq_files
end

"""
Create dataframe containing all aligned miRNA and the number of times they appear.
"""
function GenerateMiRNACounts(SAM_File::String)
    miRNA_names_list = Vector{String}()
    open(SAM_File, "r") do sam_file
	for line in eachline(sam_file)
	    miRNA_name = match(r"0\s+\Kh[a-z-A-Z0-9]+", line)
	    if !isnothing(miRNA_name)
	        push!(miRNA_names_list, miRNA_name.match)
	    end
	end
    end
    miRNA_counts_dict = countmap(miRNA_names_list)
    miRNA_counts_dataframe = DataFrame([collect(keys(miRNA_counts_dict)), collect(values(miRNA_counts_dict))], [:name, :count])
    sort!(miRNA_counts_dataframe, :count, rev = true)

    return miRNA_counts_dataframe
end

"""
    miRNADiscoveryCalculation(Trimmed_Fastq_Files::Vector{String}
                              , Sample_Names::Vector{SubString{String}}
                              )

Align trimmed fastq files to the single organism miRNA bowtie2 reference. 

The bowtie2 output SAM is input into function to generate miRNA read counts.

#Example
```julia
julia> miRNADiscoveryCalculation(["sample1.cut.fastq","sample2.cut.fastq","sample3.cut.fastq"]
                                , ["sample1", "sample2", "sample3"])
3-element Vector{DataFrame}:
469x2 DataFrame
 Row │ name             count
     │ String           Int64
─────┼────────────────────────
   1 │ hsa-miR-1307-3p      1
   2 │ hsa-miR-425-5p      58
   3 │ hsa-miR-488-5p       2
 [...]

 3-element Vector{String}:
 "sample1.miRNA.sam"
 "sample2.miRNA.sam"
 "sample3.miRNA.sam"
```
"""
function miRNADiscoveryCalculation(Trimmed_Fastq_Files::Vector{String}
                                    , Sample_Names::Vector{SubString{String}}
                                    )
    number_of_records = length(Trimmed_Fastq_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .5
                                            , "Calculating the number of miRNA present..."
                                            )

    sample_tracker = 1
    vector_of_miRNA_counts_dfs = Vector{DataFrame}()

    for fastq_file in Trimmed_Fastq_Files
        sample_name = Sample_Names[sample_tracker]

        #Align to miRBase reference
        wait(run(pipeline(
                `bowtie2
                --norc
                --threads 12
                -x data/miRNA
                -U $fastq_file
                -S $sample_name.miRNA.sam`
                , devnull)
                , wait = false
                ))

        counts_df = GenerateMiRNACounts(string(sample_name, ".miRNA.sam"))
        push!(vector_of_miRNA_counts_dfs, counts_df)
                
        next!(progress_bar_update)
        sample_tracker += 1
    end

    sam_files::Vector{String} = CaptureTargetFiles(".miRNA.sam")

    return vector_of_miRNA_counts_dfs, sam_files
end

"""
    PlotMiRNACounts(miRNA_Counts_Dfs::Vector{DataFrame}
                        , Sample_Names::Vector{SubString{String}}
                        )

Plot counts of the number of unique miRNA each sample aligned to.

#Example
```julia
julia> PlotMiRNACounts([sample1_counts_df,sample2_counts_df,sample3_counts_df]
                        , ["sample1", "sample2", "sample3"])
3-element Vector{String}:
 "sample1_miRNA_counts.csv"
 "sample2_miRNA_counts.csv"
 "sample3_miRNA_counts.csv"
```
"""
function PlotMiRNACounts(miRNA_Counts_Dfs::Vector{DataFrame}
                        , Sample_Names::Vector{SubString{String}}
                        )
    number_of_records = length(miRNA_Counts_Dfs)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .25
                                            , "Counting the miRNA in each sample..."
                                            )

    sample_tracker = 1

    for miRNA_df in miRNA_Counts_Dfs
        #Calculate the total number of reads mapped to miRNA
        total_mapped_reads = sum(miRNA_df[!, :count])

        #Determine counts per million (CPM) for each sample
        CPM::DataFrame = combine(miRNA_df
        , :count => ByRow(miRNA_reads -> round(miRNA_reads / total_mapped_reads * 10^6, digits = 2)
        ))

        #Remove "(unk)" from all miRNA names
        clean_miRNA_names::DataFrame = combine(miRNA_df
                                                , :name => ByRow(name -> replace(name, "(unk)" => ""))
                                                )

        #Reconstitute dataframe with new miRNA names, miRNA counts, and CPM column
        miRNA_df::DataFrame = hcat(clean_miRNA_names, select(miRNA_df, :count), CPM)

        #Rename column with miRNA names and CPM data
        rename!(miRNA_df, :count_function => :CPM, :name_function => :name)

        #Check column containing miRNA counts and adds miRNA above a set threshold
        threshold_of_one = length(filter(>=(1), miRNA_df[!, :count]))
        threshold_of_three = length(filter(>=(3), miRNA_df[!, :count]))
        threshold_of_five = length(filter(>=(5), miRNA_df[!, :count]))
        threshold_of_ten = length(filter(>=(10), miRNA_df[!, :count]))

        sample_name = Sample_Names[sample_tracker]

        #Barplot colors
        colors = [:grey88, :skyblue2, :peachpuff, :lightsalmon]

        #Set x and y-axis values; y-axis is a vector containing each sample's counted miRNAs
        x, y = 1:4, [threshold_of_one, threshold_of_three, threshold_of_five, threshold_of_ten]

        fig = Figure(resolution = (1800, 1200))
        increments = round(threshold_of_one / 10, sigdigits=2)
        ax = Axis(fig[1, 1]
                , xticks = (1:4
                , ["Threshold 1", "Threshold 3", "Threshold 5", "Threshold 10"])
                , yticks = 0:increments:threshold_of_one
                , ylabel = "Number of miRNA"
                , title = string(sample_name, " miRNA Counts")
                )
        CairoMakie.ylims!(ax, 0, threshold_of_one * 1.05)

        #Make sample barplot
        barplot!(ax
                , x
                , y
                , strokewidth = 1
                , color = colors
                , strokecolor = :black
                , bar_labels = :y
                )

        save(string(sample_name, "_miRNA_counts.pdf"), fig)

        #Write output file containing miRNA names, read count, and CPM data
        CSV.write(string(sample_name, "_miRNA_counts.csv"), miRNA_df)
        
        next!(progress_bar_update)
        sample_tracker += 1
    end

    pdf_miRNA_files::Vector{String} = CaptureTargetFiles("_miRNA_counts.pdf")

    #Merge output PDFs into one file
    merge_pdfs([pdf_miRNA_files...], "Small_RNA_miRNA_plots.pdf")
    
    #Remove individual sample plots
    TrashRemoval(pdf_miRNA_files)

    miRNA_counts_files::Vector{String} = CaptureTargetFiles("_miRNA_counts.csv")

    return miRNA_counts_files
end

"""
Use aligned SAM files to create a read length distribution file.
"""
function CalculateReadLengthDistribution(SAM_Files::Vector{String}
                                        , Sample_Names::Vector{SubString{String}}
                                        )
    number_of_records = length(SAM_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .5
                                            , "Calculating Read Length Distribution..."
                                            )

    sample_tracker = 1

    for sam_file in SAM_Files
        sample_name = Sample_Names[sample_tracker]

	    #=
        If the SAM file is empty or almost empty, it will cause errors 
        downstream so the script exits.
        =#
        read_length_distibution = 
        try
            read(pipeline(
            `samtools 
            stats 
            -@ 12
            $sam_file`
            , `grep ^RL`
            , `cut -f 2-`)
            , String
            )
        catch
            println("Error, SAM file contains no read length distribution")
            println("Exiting program...")
            exit()
        end

        next!(progress_bar_update)

        #Write read length data to output file
        output_length_file::IOStream = open("read_lengths_$sample_name.csv", "w")
        write(output_length_file, string("Length", "\t", "Reads", "\n"))
        write(output_length_file, read_length_distibution)

        close(output_length_file)
        sample_tracker += 1
    end

    length_files::Vector{String} = CaptureTargetFiles("read_lengths_")

    return length_files
end

"""
Create barplot of each sample's fragment lengths based on the number of reads found at each length
"""
function PlotFragmentLengths(Length_Files::Vector{String}
                            , Sample_Names::Vector{SubString{String}}
                            )
    number_of_records = length(Length_Files)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , 1
                                            , "Creating Read Length Plot..."
                                            )

    sample_tracker = 1

    for file in Length_Files
        sample_name = Sample_Names[sample_tracker]
        length_file::DataFrame = CSV.read(file, DataFrame)

        #Isolate fragment lengths and the number of reads at each length
        fragment_lengths::Vector{Int64} = length_file[!, :Length]
        reads_per_length::Vector{Int64} = length_file[!, :Reads]
        next!(progress_bar_update)

        plot = Figure(resolution = (1800, 1200))
        next!(progress_bar_update)

        Axis(plot[1, 1]
        , xticks = 15:5:maximum(fragment_lengths)
        , xlabel = "Fragment Lengths"
        , ylabel = "Number of Reads"
        , title = sample_name
        )
        next!(progress_bar_update)

        #Make sample barplot
        barplot!(fragment_lengths
        , reads_per_length
        , fillto = -1
        , color = fragment_lengths
        , strokecolor = :black
        , strokewidth = 1
        )

        save(string(sample_name, "_fragment_lengths.pdf"), plot)

        next!(progress_bar_update)
        sample_tracker += 1
    end

    pdf_length_files::Vector{String} = CaptureTargetFiles("_fragment_lengths.pdf")

    merge_pdfs([pdf_length_files...], "Small_RNA_fragment_length_plots.pdf")

    TrashRemoval(pdf_length_files)
end

"""
    AlignWithBowtie2(Fastq_File::String, Read_Count_Dict::Dict{String, Int64})

Align fastq with specified reference RNA types.

Divide reads aligned with total reads to calculate percent alignment. Store 
alignment information in a dictionary.

#Example
```julia
julia> AlignWithBowtie2("sample1.cut.fastq"
                        , Dict("sample1" => 6390314, "sample2" => 5000000, "sample3" => 7052928))
Dict{String, Float64} with 5 entries:
  "miRNA" => 65.0
  "tRNA" => 4.2
  "piRNA" => 1.2
  "snoRNA" => 1.9
  "rRNA" => 3.3
```
"""
function AlignWithBowtie2(Fastq_File::String, Read_Count_Dict::Dict{String, Int64})
    sample_name = first(split(Fastq_File, "."))

    #RNA types to be aligned to
    reference_RNA_types = ["miRNA", "tRNA", "piRNA", "snoRNA", "rRNA"]

    #Dictonary to hold each RNA types' alignment data
    metrics_dict = Dict{String, Number}()

    for reference_RNA in reference_RNA_types

        #Can skip aligning sample with miRNA reference since that was already done.
        sam_filename = ""
        if reference_RNA != "miRNA"
            sam_filename = string(sample_name, ".metrics.sam")
            #Align with bowtie2
            wait(run(pipeline(
            `bowtie2
            --threads 12
            --norc
            -x data/$reference_RNA
            -U $Fastq_File
            -S $sam_filename`
            , devnull)
            , wait = false
            ))
        else
            sam_filename = string(sample_name, ".miRNA.sam")
        end

        #Count number of read alignments
        aligned_reads_string = read(`
        samtools 
        view 
        -c 
        -F 4 
        $sam_filename`
        , String
        )
        aligned_reads = parse(Float64, aligned_reads_string)

        #Transfer unmapped reads to SAM file
        wait(run(pipeline(`
        samtools 
        view 
        -f 4 
        -o $sample_name.aligned.metrics.sam
        $sam_filename`
        , devnull)
        , wait = false
        ))

        #Convert SAM file back into a fastq
        wait(run(pipeline(`
        samtools 
        fastq 
        -@ 12 
        -n 
        $sample_name.aligned.metrics.sam`
        , stdout = Fastq_File
        ), wait = false
	))

        percent_alignment = round(100 * aligned_reads / Read_Count_Dict[sample_name], digits = 3)
        metrics_dict[reference_RNA] = percent_alignment
    end

    return metrics_dict
end

"""
Calculate the bowtie2 alignment metrics for all samples. Return filename for plotting.
"""
function CalculateMetrics(Trimmed_Fastqs::Vector{String}
                            , Read_Count_Dict::Dict{String, Int64}
                            , Sample_Names::Vector{SubString{String}}
                            , Dimer_Count_Dict::Dict{String, Float64}
                            , Q_Score_Dict::Dict{String, Number}
                            )
    number_of_records = length(Trimmed_Fastqs)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .5
                                            , "Calculating Sample Metrics..."
                                            )

    sample_tracker = 1

    output_metrics_file::IOStream = open("Alignment_Metrics.csv", "w")

    #Write metrics header information
    write(output_metrics_file, string(
        "Sample", ","
        , "Read Count", ","
        , "Aligned Reads", ","
        , "Average Q score", ","
        , "Dimer", ","
        , "<16 bp", ","
        , "miRNA", ","
        , "tRNA", ","
        , "piRNA", ","
        , "snoRNA", ","
        , "rRNA", ","
        , "Unaligned Reads"
        , "\n"
        ))

    for fastq_file in Trimmed_Fastqs
        sample_name = Sample_Names[sample_tracker]

        #Get sample q_score
        average_q_score = Q_Score_Dict[sample_name]

        #Gather the number of canonical dimer reads
        percent_dimer = Dimer_Count_Dict[sample_name]

        #Get number of short fragments from cutadapt output
        short_fragments_string = match(r"Reads that.+\(\K.+(?=%)"
                                        , readchomp(`cat $sample_name.cutadapt_information`)
                                        ).match
        short_fragments = abs(parse(Float64, short_fragments_string) - percent_dimer)

        #Function to align each fastq to each RNA type and calculate the number of reads mapped
        metrics_dict::Dict{String, Number} = AlignWithBowtie2(fastq_file, Read_Count_Dict)

        aligned_reads = round(values(metrics_dict) |> sum, digits = 2)
        unaligned_reads = round(100 - aligned_reads - percent_dimer - short_fragments, digits = 2)
        
        #Write metrics data to output file
        write(output_metrics_file, string(
        sample_name, ","
        , Read_Count_Dict[sample_name], ","
        , aligned_reads, ","
        , average_q_score, ","
        , round(percent_dimer, digits = 2), ","
        , round(short_fragments, digits = 2), ","
        , metrics_dict["miRNA"], ","
        , metrics_dict["tRNA"], ","
        , metrics_dict["piRNA"], ","
        , metrics_dict["snoRNA"], ","
        , metrics_dict["rRNA"], ","
        , unaligned_reads
        , "\n"
        )) 

        next!(progress_bar_update)
        sample_tracker += 1
    end

    close(output_metrics_file)
    
    #Return the metrics file name
    return "Alignment_Metrics.csv"
end


"""
Create violin plot of each sample's RNA aligment metrics
"""
function ViolinPlotMetrics(Metrics_File::String)

    #Import sample metrics file
    metrics_file::DataFrame = CSV.read(Metrics_File, DataFrame)

    #Create new dataframe without sample names and add sample and column names to arrays
    select!(metrics_file, Not([:Sample, Symbol("Read Count")]))
    column_names::Vector{String} = names(metrics_file)
    num_of_cols = last(size(metrics_file))

    number_of_records = length(column_names)
    progress_bar_update = ProgressBarUpdate(number_of_records, .25
    , "Creating Violin Plots..."
    )
    
    #Gather each metric's column data for plotting
    metrics_columns(a) = [number for number in metrics_file[!, a]]

    #Make empty figure and set x and y axes
    #x-axis is vector of integers from 1 to the number of metrics
    #y-axis is vector of vectors containing the data from each metric
    fig = Figure(resolution = (1800, 1200))
    x, y = 1:num_of_cols, [metrics_columns(column) for column in 1:num_of_cols]
    
    #Create violin plot with data from all samples analyzed
    for column in 1:num_of_cols
        ax = Axis(fig[1, column]
        , yticks = 0:5:100
        , ylabel = "Percent"
        , title = column_names[column]
        )
        CairoMakie.ylims!(ax, 0, 100)

        #Make violin plot with combined sample data
        CairoMakie.violin!(fig[1, column]
        , repeat([x[column]]
        , first(size(metrics_file))
	)
        , y[column]
        , show_median=true
        )
        next!(progress_bar_update)

        #Save file once all columns have been added
        if column == num_of_cols
            save("Violin_plot_metrics.png", fig)
        end

        next!(progress_bar_update)
    end
end

"""
Create barplot of each sample's RNA aligment metrics.
"""
function PlotMetrics(Metrics_File::String, Library_Type::Dict{String, String}
                    , Sample_Names::Vector{SubString{String}}
                    )

    metrics_file::DataFrame = CSV.read(Metrics_File, DataFrame)

    select!(metrics_file, Not([:Sample, Symbol("Read Count")]))
    column_names::Vector{String} = names(metrics_file)
    num_of_cols = last(size(metrics_file))

    colors = [:snow3, :honeydew2, :lightpink, :bisque2, :deepskyblue, :aquamarine3
    , :peachpuff, :paleturquoise1, :orange1, :lightsalmon2
    ]
    
    sample_rows(a) = [values for values in metrics_file[a, :]]

    number_of_records = length(Sample_Names)
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .25
                                            , "Creating Metrics Plots..."
                                            )
    
    for sample in eachindex(Sample_Names)
        sample_name = Sample_Names[sample]

        #Set axis values; y-axis is a vector containing each sample's calculated metrics
        x, y = 1:num_of_cols, sample_rows(sample)
        fig = Figure(resolution = (1800, 1200))

        ax = Axis(fig[1, 1]
        , xticks = (1:num_of_cols
        , column_names)
        , yticks = 0:10:100
        , ylabel = "Percent"
        , title = string(sample_name, " Metrics - ", Library_Type[sample_name])
        )

        #Set y-axis limits
        CairoMakie.ylims!(ax, 0, 105)

        barplot!(ax
        , x
        , y
        , strokewidth = 1
        , strokecolor = :black
        , color = colors
        , bar_labels = :y
        )

        next!(progress_bar_update)

        save(string(sample_name, "_metrics.pdf"), fig)
    end

    pdf_metrics_files::Vector{String} = CaptureTargetFiles("_metrics.pdf")

    merge_pdfs([pdf_metrics_files...], "Small_RNA_metrics_plots.pdf")
    
    TrashRemoval(pdf_metrics_files)
end

"""
Find all miRNAs the samples have in common. Calculate the reads per million reads (RPM).
"""
function FindCommonMiRNAs(miRNA_Counts_Files::Vector{String}
                            , Read_Count_Dict::Dict{String, Int64}
                            , Sample_Names::Vector{SubString{String}}
                            )
    #Vector of dictionaries containing the miRNA counts from each sample
    miRNA_info = Vector{Dict{String, Int64}}()
    RPM_info = Vector{Dict{String, Number}}()
    full_miRNA_names_list = Vector{String}()

    for file in miRNA_Counts_Files
        miRNA_file::IOStream = open(file, "r")
        sample_name = first(split(file, "_"))

        #Dictionaries to hold each miRNA and its associated read count
        sample_miRNA_counts_dictionary = Dict{String, Int64}()
        sample_miRNA_RPM_dictionary = Dict{String, Number}()

        #Skip header line
        readline(miRNA_file)

        for line in eachline(miRNA_file)
            split_line = split(line, ",")
            miRNA_count = parse(Int64, split_line[2])
            miRNA_name = first(split_line)
            RPM = round(miRNA_count / (Read_Count_Dict[sample_name] / 10^6), digits = 4)

            sample_miRNA_counts_dictionary[miRNA_name] = miRNA_count
            sample_miRNA_RPM_dictionary[miRNA_name] = RPM
        end

        push!(miRNA_info, sample_miRNA_counts_dictionary)
        push!(RPM_info, sample_miRNA_RPM_dictionary)
        
        append!(full_miRNA_names_list, keys(sample_miRNA_counts_dictionary))

        close(miRNA_file)
    end

    #Count number each miRNA's occurances
    miRNA_names_dict = countmap(full_miRNA_names_list)

    #Find miRNA present in all samples, .i.e. have counts equal to the number of samples analyzed
    miRNAs_in_common = filter(miRNA -> last(miRNA) === length(Sample_Names), miRNA_names_dict) |> keys

    return miRNAs_in_common, miRNA_info, RPM_info
end

"""
Write miRNAs all samples have in common to output file.
"""
function WriteCommonMiRNAFile(miRNA_Names::Base.KeySet{String, Dict{String, Int64}}
                                , miRNA_Info::Vector{Dict{String, Int64}}
                                , RPM_Info::Vector{Dict{String, Number}}
                                , Sample_Names::Vector{SubString{String}}
                                )
    counter = 1
    
    for sample in Sample_Names
        common_miRNA_file::IOStream = open(string(sample, "_common_miRNAs.tsv"), "w")
        common_miRNA_file_RPM::IOStream = open(string(sample, "_common_miRNAs_RPM.tsv"), "w")

        if counter == 1
            write(common_miRNA_file, string("miRNA", "\t", sample, "\n"))
            write(common_miRNA_file_RPM, string("miRNA", "\t", sample, "\n"))
            for miRNA in miRNA_Names
                write(common_miRNA_file, string(miRNA, "\t", miRNA_Info[counter][miRNA], "\n"))
                write(common_miRNA_file_RPM, string(miRNA, "\t", RPM_Info[counter][miRNA], "\n"))
            end
        else
            write(common_miRNA_file, string(sample, "\n"))
            write(common_miRNA_file_RPM, string(sample, "\n"))
            for miRNA in miRNA_Names
                write(common_miRNA_file, string(miRNA_Info[counter][miRNA], "\n"))
                write(common_miRNA_file_RPM, string(RPM_Info[counter][miRNA], "\n"))
            end
        end

        counter += 1

        close(common_miRNA_file)
        close(common_miRNA_file_RPM)
    end

    common_files::Vector{String} = CaptureTargetFiles("_common_miRNAs.tsv")
    common_RPM_files::Vector{String} = CaptureTargetFiles("_common_miRNAs_RPM.tsv")

    #Combine all the separate common miRNA files into one file
    run(pipeline(`paste $common_files`, stdout = "Common_miRNAs.tsv"))
    run(pipeline(`paste $common_RPM_files`, stdout = "Common_RPM_miRNAs.tsv"))

    #Run principal component analysis and UMAP on common miRNA
    PlotClustering("Common_miRNAs.tsv")

    TrashRemoval(common_files)
    TrashRemoval(common_RPM_files)
end

"""
Plot principal component analysis (PCA) and miRNA counts UMAP of all libraries if possible.
"""
function PlotClustering(Common_miRNA_File::String)
    common_miRNA_counts = DataFrame(CSV.File(Common_miRNA_File))
    sample_names = DataFrames.names(common_miRNA_counts, Not([:miRNA]))
    miRNA_names = common_miRNA_counts[!, :miRNA]

    #There must be at least two samples in order to perform the principal component analysis
    #Must also have at least one miRNA in common
    if last(size(common_miRNA_counts)) > 2 && first(size(common_miRNA_counts)) > 1
        common_miRNA_counts_matrix = Matrix{Float64}(select(common_miRNA_counts, Not([:miRNA])))

        #Make PCA covariance matrix for PCA plot
        PCA_matrix = fit(PCA, common_miRNA_counts_matrix; maxoutdim = 20)
        transformed_counts = predict(PCA_matrix, common_miRNA_counts_matrix)

        #Make PCA covariance matrix for UMAP plot
        UMAP_PCA_matrix = fit(PCA, common_miRNA_counts_matrix'; maxoutdim = 20)
        transposed_transformed_counts = predict(UMAP_PCA_matrix, common_miRNA_counts_matrix')

        #Plot PCA only if there are at least two principal components
        if first(size(transformed_counts)) >= 2
            sample_name_df = DataFrame("samples" => sample_names)
            miRNA_names_df = DataFrame("miRNA" => miRNA_names)

            #Write PCA information to file for cluster tracking
            PCA_values = hcat(sample_name_df, DataFrame(transformed_counts[1:2, :]', ["PC1", "PC2"]))
            CSV.write("PCA_information.csv", PCA_values)
            
            #Plot PCA
            Plots.scatter(size = (1200, 800), dpi = 300, titlefont = (16, "Computer Modern")
                            , xlabel = "PC1", ylabel = "PC2", title = "Common miRNA: PCA"
                            , transformed_counts[1, :], transformed_counts[2, :]
                            , left_margin = 23mm, right_margin = 8mm, bottom_margin = 8mm
                            , leg = false)
            savefig("common_miRNA_PCA.png")
        end
        
        #Plot UMAP only if there are at least three principal components
        if first(size(transposed_transformed_counts)) >= 3
            sample_name_df = DataFrame("samples" => sample_names)
            miRNA_names_df = DataFrame("miRNA" => miRNA_names)

            #Create low dimensional embedding for UMAP
            near_neighbors = first(size(common_miRNA_counts_matrix)) - 1
            if near_neighbors < 15
                embedding = umap(transposed_transformed_counts
                                , 2
                                ; n_neighbors = near_neighbors
                                , min_dist = 0.05)
            else
                embedding = umap(transposed_transformed_counts
                                , 2
                                ; n_neighbors = 15
                                , min_dist = 0.05)
            end

            #Write UMAP information to file for cluster tracking
            UMAP_values = hcat(miRNA_names_df, DataFrame(embedding', ["UMAP1", "UMAP2"]))
            CSV.write("UMAP_information.csv", UMAP_values)

            #=
            Plot UMAP of common miRNAs. Color UMAP clusters based on K-medoids clustering; 
            choose number of clusters based on the number of principal components.
            =#
            distance_matrix = pairwise(Euclidean(), embedding, embedding)
            kmedoids_cluster_colors = kmedoids(distance_matrix, first(size(transposed_transformed_counts)))
            Plots.scatter(size = (1200, 800), embedding[1, :], embedding[2, :]
                                        , title="Common miRNA: UMAP", left_margin = 13mm
                                        , bottom_margin = 10mm, dpi = 300
                                        , marker_z = kmedoids_cluster_colors.assignments
                                        , color = :lighttest, xlabel = "UMAP1", ylabel = "UMAP2"
                                        , titlefont = (16, "Computer Modern"), leg = false
                                        , markersize = 9, markerstrokewidth = 0.1)
            savefig("common_miRNA_UMAP.png")
        end
    end
end

"""
    CqVsMiRNACountCorrelation(qPCR_Data_File::String
                                , Common_miRNA_RPM_File::String
                                , Common_miRNA_File::String
                                , Sample_Names::Vector{SubString{String}}
                                )

The common miRNA counts are used to calculate the Pearson's correlation coefficient.

Filter the qPCR data* to remove measurements without a quantification cycle (Cq) value.
Compare the filtered list with the miRNA present in the common miRNA files and filter each 
common miRNA file. Calculate the average Cq value for a given miRNA from each replicate
and determine the Pearson's correlation coefficient.

Samples with a Pearson's coefficient < -.5 are plotted using linear regression. 

#Example
```julia
julia> CqVsMiRNACountCorrelation("qpcr_raw_data.csv", "Common_RPM_miRNAs.csv"
                                    , "Common_miRNAs.csv", ["sample1", "sample2", "sample3"])
3-element Vector{String}:
 "sample1"
 "sample2"
 "sample3"

 500x5 DataFrame
Row  │ miRNA            Cq       sample1  sample2  sample3  
     │ String31         Float64  Int64         Int64         Int64         
─────┼─────────────────────────────────────────────────────────────────────
   1 │ hsa-miR-425-5p   25.9225           220           132           198  
   2 │ hsa-miR-92b-3p   25.38            1318          1164          1362   
   3 │ hsa-miR-22-3p    24.215           2987          2878          2482
 [...]
```
"""
function CqVsMiRNACountCorrelation(qPCR_Data_File::String
                                    , Common_miRNA_RPM_File::String
                                    , Common_miRNA_File::String
                                    , Sample_Names::Vector{SubString{String}}
                                    )
    number_of_records = 4
    progress_bar_update = ProgressBarUpdate(number_of_records
                                            , .25
                                            , "Calculating qPCR data vs miRNA count correlation data..."
                                            )

    #Import qPCR data file and both common miRNA files
    qPCR_data_table_full = DataFrame(CSV.File(qPCR_Data_File))
    common_miRNA_counts = DataFrame(CSV.File(Common_miRNA_File))
    RPM_dataframe = DataFrame(CSV.File(Common_miRNA_RPM_File))
    next!(progress_bar_update)

    #Filter out measurements with no Cq value
    qPCR_data_table::DataFrame = filter(:Cq => !=("NA"), qPCR_data_table_full)

    #Gather intersection between all common miRNA and miRNA list after removing filtered miRNA with no Cq value
    common_miRNA = common_miRNA_counts[!, :miRNA] ∩ qPCR_data_table[!, :miRname]

    #Need to remove the rows with the miRNA filtered out above from the common miRNA RPM file
    filtered_RPM_dataframe = filter(:miRNA => in(common_miRNA), RPM_dataframe)
    next!(progress_bar_update)

    #Find Cqs for all common miRNAs, convert to Floats, and calculate the replicates' average Cq
    average_cq_vector::Vector{Vector{Float64}} = map(
        miRNA -> filter(:miRname => ==(miRNA), qPCR_data_table)[!, :Cq] 
        .|> Cq -> parse(Float64, Cq)
        , common_miRNA
        )
    average_cqs_values::Vector{Float64} = map(
        cq_vector -> sum(cq_vector) |> num -> num / length(cq_vector), average_cq_vector
        )

    next!(progress_bar_update)

    #Create new dataframe with the names of all common miRNAs and their Cq values
    cq_miRNA_dataframe = zip(common_miRNA, average_cqs_values) |> DataFrame

    #Add Cq values and miRNA names to RPM data
    combined_cq_RPM_dataframe::DataFrame = hcat(cq_miRNA_dataframe, filtered_RPM_dataframe[!, 2:end])
    rename!(combined_cq_RPM_dataframe, :1 => :miRNA, :2 => :Cq)

    #Pearson's correlation coefficient analysis between the qPCR values and the log10 of miRNA counts
    correlation_coeffs::Vector{Float64} = map( 
        sample -> cor(combined_cq_RPM_dataframe[!, :Cq], map(log10, combined_cq_RPM_dataframe[!, sample]))
        , Sample_Names
        )

    correlation_dataframe = DataFrame(hcat(Sample_Names, correlation_coeffs), [:sample, :corr_coeff])

    #Find samples with correlation values < -.6
    negative_correlation_samples::Vector{String} = Sample_Names[correlation_coeffs.<(-.6)]

    #Write correlation values to output file
    CSV.write("Sample_correlation_coefficients.csv", correlation_dataframe)
    next!(progress_bar_update)

    return negative_correlation_samples, combined_cq_RPM_dataframe
end

"""
Make linear regression plot for samples correlated with the qPCR data.
"""
function PlotLinearRegressionCqVsMiRNAReadCount(Negative_Correlation_Samples::Vector{String}
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
    
        #Make new dataframe with Cq values and log10 converted miRNA counts
        sample_cq_dataframe::DataFrame = hcat(DataFrame(name = x_values), DataFrame(Cq = y_values))
    
        #Create linear regression model
        linear_model = lm(@formula(Cq ~ name), sample_cq_dataframe)
    
        #Use linear model to predict the Cq values expected based on the miRNA counts
        predictions = predict(linear_model, sample_cq_dataframe, interval = :confidence, level = 0.95)
    
        #Create scatter plot with Cq vs the log of miRNA read counts
        linear_plot = @df sample_cq_dataframe Plots.scatter(:name, :Cq, leg = false, markersize = 7)
    
        #Sort data on log10 miRNA counts; reverse order since the values are negatively correlated
        sorted_predictions::DataFrame = predictions[sortperm(predictions[!, :prediction], rev = true), : ]
        sorted_x_values::Vector{Float64} = sort(x_values)
    
        #Extends plot y-axis by 1/7 each side
        dY::Float64 = first(diff([extrema(y_values)...])) / 7  
        y1, y2 = (-dY, dY) .+ extrema(y_values)
    
        #Goodness of fit
        M::Matrix{Float64} = [ones(length(x_values)) x_values]
        a, b = M \ y_values
        Ŝ = round(sqrt(sum((y_values - b*x_values .- a).^2) / (length(x_values) - 2)), digits = 3)

        #Make plot Cq vs log10 miRNA read counts with R-squared, Goodness of fit, and Pearson's coefficient in title
        Plots.plot!(size = (1200, 800), title = string(sample, " Linear Regression", "\n", " R-Squared = "
        , round(r2(linear_model), digits = 3), "\n", "Goodness of Fit = ", Ŝ, "\n", "Pearson \U03c1 = "
        , round(cor(x_values, y_values), digits = 3)), xlabel = "log10 miRNA Read Count"
        , ylabel = "Average Cq", dpi = 300, titlefont = (14, "Computer Modern")
        , ylims = (y1, y2), guidefont = (12, :black), tickfont = (11, :black), left_margin = 13mm, bottom_margin = 10mm
        , linear_plot, sorted_x_values, sorted_predictions.prediction, linewidth = 4
        , ribbon = (sorted_predictions.prediction .- sorted_predictions.lower
        , sorted_predictions.upper .- sorted_predictions.prediction
        ))
    
        savefig(string(sample, "_miRNA_correlation.png"))
    end
end

"""
Spot remove unecessary intermediate files.
"""
function TrashRemoval(Files_to_Delete::Vector{String})
    for file in Files_to_Delete
        rm(file)
    end
end

"""
Remove all intermediate files.
"""
function GarbageCollection()
    Files_to_Delete = Set(vcat(
    CaptureTargetFiles(".bam")
    ,CaptureTargetFiles(".cut.fastq")
    ,CaptureTargetFiles("sub_")
    ,CaptureTargetFiles("read_lengths_")
    ,CaptureTargetFiles(".miRNA.")
    ,CaptureTargetFiles("v22_cluster")
    ,CaptureTargetFiles(".cutadapt_information")
    ,CaptureTargetFiles(".sam")
    ))

    for file in Files_to_Delete
        rm(file)
    end
end

function julia_main()::Cint

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

    Fastqs = CaptureTargetFiles("_R1_001.fastq.gz")
    Sample_Names = map(sample -> first(split(sample, "_")), Fastqs)
    Library_Type, Read_Count_Dict, Dimer_Count_Dict, Q_Score_Dict = DetermineLibraryType(Fastqs
                                                                                        , Sample_Names)
    Trimmed_Fastqs = TrimAdapters(Fastqs, Library_Type, Sample_Names)
    miRNA_Counts_Dfs, SAM_Files = miRNADiscoveryCalculation(Trimmed_Fastqs
                                                                                        , Sample_Names)

    Length_Files = CalculateReadLengthDistribution(SAM_Files, Sample_Names)
    PlotFragmentLengths(Length_Files, Sample_Names)
    miRNA_Counts_Files = PlotMiRNACounts(miRNA_Counts_Dfs, Sample_Names)
    Metrics_File = CalculateMetrics(Trimmed_Fastqs
                                            , Read_Count_Dict
                                            , Sample_Names
                                            , Dimer_Count_Dict
                                            , Q_Score_Dict)
    PlotMetrics(Metrics_File, Library_Type, Sample_Names)
    ViolinPlotMetrics(Metrics_File)
    Full_miRNA_Names_List, miRNA_Info, RPM_Info = FindCommonMiRNAs(miRNA_Counts_Files
                                                                    , Read_Count_Dict
                                                                    , Sample_Names)
    WriteCommonMiRNAFile(Full_miRNA_Names_List, miRNA_Info, RPM_Info, Sample_Names)
    Correlated_Samples, Cq_RPM_Dataframe = CqVsMiRNACountCorrelation("qpcr_raw_data.csv"
                                                                    , "Common_RPM_miRNAs.tsv"
                                                                    , "Common_miRNAs.tsv"
                                                                    , Sample_Names
                                                                    )
    PlotLinearRegressionCqVsMiRNAReadCount(Correlated_Samples, Cq_RPM_Dataframe)
    GarbageCollection()

    return 0
end

end
