# NEXTFLEX Small RNA Analysis

[![CI](https://github.com/cookienocreams/small_rna/actions/workflows/CI.yml/badge.svg)](https://github.com/cookienocreams/small_rna/actions/workflows/CI.yml)
[![Documentation](https://github.com/cookienocreams/small_rna/actions/workflows/documentation.yaml/badge.svg)](https://github.com/cookienocreams/small_rna/actions/workflows/documentation.yaml)
[![pages-build-deployment](https://github.com/cookienocreams/small_rna/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/cookienocreams/small_rna/actions/workflows/pages/pages-build-deployment)

## Script dependencies
- bowtie2
- samtools
- cutadapt
  
## Installation instructions

Create executable to run on local machine using the `julia` library `PackageCompiler`:

You will need to have Julia installed on your computer before starting. Julia can be installed from here: https://julialang.org/downloads/

 ## Note: This currently only compiles under the Upcoming release: v1.10.0-beta2
 There is a bug in the Plots library that prevents packages from being complied using PackageCompiler, see [Plots.lj issue #825](https://github.com/JuliaLang/PackageCompiler.jl/issues/825).
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
Which outputs:
```
usage: Small RNA Analysis [-l LIBRARY] [-f FASTA] [-O ORGANISM]
                        [-p THREADS] [-h]

Basic analysis of NEXTFLEX Small RNA libraries.

optional arguments:
  -l, --library LIBRARY
                        The type of sample library, either 'v4' or
                        'v3'. (default: "v4")
  -f, --fasta FASTA     The full path to a miRNA fasta file, e.g.,
                        /home/user/mirna.fa. (default:
                        "data/mirgene.fas")
  -O, --organism ORGANISM
                        Abbreviated name of organism, e.g., 'hsa' for
                        human or 'mmu' for mouse. This should match
                        the standard three letter abbreviation found
                        in miRNA databases such as miRBase and
                        MirGeneDB. (default: "hsa")
  -p, --threads THREADS
                        The number of processors to use for alignment.
                        (type: Int64, default: 12)
  -h, --help            show this help message and exit

```

## Basic usage

The app can be run using the `small_rna` executable in the `/home/user/small_rna_app/bin` folder in a folder containing fasta to be analyzed. Note that currently the data reference files can be downloaded and placed into the folder with the fastqs to be analyzed for faster setup.

```bash
cd fastqs
/home/user/small_rna_app/bin/small_rna -l v4
```
