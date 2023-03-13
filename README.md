# NEXTFLEX_SRNA_Analysis
Generates alignment metrics and plots for NEXTFLEX small RNA libraries

## Script dependencies
- bowtie2
- samtools
- cutadapt
- ghostscript - install with `sudo apt install ghostscript`

## Installation instructions

Create executable to run on local machine using the `julia` library `PackageCompiler`:

```julia
using Pkg
Pkg.add("PackageCompiler")
```
Note: To use PackageCompiler a C-compiler needs to be available.

```julia
using PackageCompiler
Pkg.generate("small_RNA_analysis")
```
This will create a small RNA module which contains a `Project.toml` file and a `src` folder.
Change directory into the newly created package and activate the new environment.

```julia
cd("small_RNA_analysis")
Pkg.activate("./")
```
The next step is to install all libraries and their dependencies.

```julia
Pkg.add.(["CairoMakie", "CSV", "DataFrames", "ProgressMeter", "MultivariateStats", "StatsBase", "StatsPlots", "Statistics", "MLBase", "GLM", "Measures", "GZip", "UMAP", "Clustering", "Distances", "PDFmerger"])
```
Once the libraries are installed, copy the `small_RNA_analysis.jl` in this respository into the `small_RNA_analysis/src` folder, replacing the auto-generated file.
The last step is to create the precompiled executable. Make sure to set the correct paths for your machine.

```julia
PackageCompiler.create_app("/path/to/small_RNA_analysis", "/home/user/sRNA_app", incremental=true, precompile_execution_file="/path/to/small_RNA_analysis/src/small_RNA_analysis.jl", include_lazy_artifacts=true)
```

The app can be run using the small_RNA_analysis executabile in the `/home/user/sRNA_app/bin` folder in a folder containing fastqs to be analyzed. Note that currently the `data` and `qpcr_raw_data.csv`* files must be downloaded and placed into the analysis folder too.

```bash
cd fastqs
/home/user/sRNA_app/bin/small_RNA_analysis
```

*qPCR data from here: Maguire, S.et al. (2020). A low-bias and sensitive small RNA library preparation method using randomized splint ligation. Nucleic Acids Research, 48(14). https://doi.org/10.1093/nar/gkaa480. The data is for human brain samples, but can be swapped for qPCR data for other sample types.
