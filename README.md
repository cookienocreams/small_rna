# NEXTFLEX Small RNA Analysis

[![CI](https://github.com/cookienocreams/small_rna/actions/workflows/CI.yml/badge.svg)](https://github.com/cookienocreams/small_rna/actions/workflows/CI.yml)

## Script dependencies
- bowtie2
- samtools
- cutadapt
  
## Installation instructions

Create executable to run on local machine using the `julia` library `PackageCompiler`:

You will need to have Julia installed on your computer before starting. Julia can be installed from here: https://julialang.org/downloads/

 ## Note: This currently only supports the Upcoming release: v1.10.0-beta2
 There is a bug in the Plots library that prevents packages from being complied using PackageCompiler.
 It should be fixed in v1.10, once that releases. Versions < v1.8.5 may also work.

Download the git repository using git or manually and change into the repository folder.
```bash
git clone https://github.com/cookienocreams/small_rna.git small_rna
cd small_rna
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
PackageCompiler.create_app("./", "/home/user/small_rna_app", incremental=true, precompile_execution_file="./src/small_rna.jl", include_lazy_artifacts=true)
```

There are numerous options that can be changed if desired. Use `-h` or `--help` flags to see options.

```bash
/home/user/small_rna_app/bin/small_rna --help
```

## Basic usage

The app can be run using the `small_rna` executable in the `/home/user/small_rna_app/bin` folder in a folder containing fasta to be analyzed. Note that currently the data reference files can be downloaded and placed into the folder with the fastqs to be analyzed for faster setup.

```bash
cd fastqs
/home/user/small_rna_app/bin/small_rna -l v4
```
