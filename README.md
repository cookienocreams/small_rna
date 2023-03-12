# NEXTFLEX_SRNA_Analysis
Generates alignment metrics and plots for NEXTFLEX small RNA libraries

Script dependencies
- bowtie2
- samtools
- cutadapt
- ghostscript - install with sudo apt install ghostscript

Create executable to run on local machine:
Step 1: using Pkg; Pkg.add("PackageCompiler")

Step 2: using PackageCompiler

Step 3: Pkg.generate("small_RNA_analysis")

Step 4: cd("small_RNA_analysis")

Step 5: Pkg.activate("./")

Step 6: Pkg.add.(["CairoMakie", "CSV", "DataFrames", "ProgressMeter", "MultivariateStats", "StatsBase", "StatsPlots", "Statistics", "MLBase", "GLM", "Measures", "GZip", "UMAP", "Clustering", "Distances"])

Step 7: Replace script in small_RNA_analysis/src with small_RNA_analysis.jl in this respository.

Step 8: PackageCompiler.create_app("/path/to/small_RNA_analysis", "/home/user/sRNA_app", incremental=true, precompile_execution_file="/path/to/small_RNA_analysis/src/small_RNA_analysis.jl")

The app can be run using the small_RNA_analysis executabile in the /home/user/sRNA_app/bin folder in a folder containing fastqs to be analyzed. Note that currently the data and qpcr_raw_data.csv files must be downloaded and placed into the analysis folder too.
