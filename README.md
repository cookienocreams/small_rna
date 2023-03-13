# NEXTFLEX_SRNA_Analysis
Generates alignment metrics and plots for NEXTFLEX small RNA libraries

## Script dependencies
- bowtie2
- samtools
- cutadapt

## Installation instructions

Create executable to run on local machine using the `julia` library `PackageCompiler`:

You will need to have Julia installed on your computer before starting. Julia can be installed from here: https://julialang.org/downloads/

Download the git repository using git or manually and change into the repository folder.
```bash
git clone https://github.com/cookienocreams/NEXTFLEX_SRNA_Analysis.git small_RNA_analysis
cd small_RNA_analysis
```
Activate the downloaded `julia` environment.
```julia
using Pkg
Pkg.activate("./")
```
The next step is to install all libraries and their dependencies.
```julia
Pkg.instantiate()
```

The last step is to create the precompiled executable. Make sure to set the correct paths for your machine.

```julia
using PackageCompiler
PackageCompiler.create_app("./", "/home/user/sRNA_app", incremental=true, precompile_execution_file="./src/small_RNA_analysis.jl", include_lazy_artifacts=true)
```

The app can be run using the small_RNA_analysis executable in the `/home/user/sRNA_app/bin` folder in a folder containing fastqs to be analyzed. Note that currently the `data` and `qpcr_raw_data.csv`* files must be downloaded and placed into the folder with the fastqs to be analyzed too.

```bash
cd fastqs
/home/user/sRNA_app/bin/small_RNA_analysis
```

*qPCR data from here: Maguire, S.et al. (2020). A low-bias and sensitive small RNA library preparation method using randomized splint ligation. Nucleic Acids Research, 48(14). https://doi.org/10.1093/nar/gkaa480. The data is for human brain samples, but can be swapped for qPCR data for other sample types.
