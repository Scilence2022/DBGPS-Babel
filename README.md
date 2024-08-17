
- [Introduction](#introduction)
- [System Requirements](#system-requirements)
- [Installation and run](#install-and-run)
- [Example usage](#example-usage)
- [License](#license)


# Introduction
DBGPS-greedy-path is a greedy strand assembler designed for Babel-DNA encryption [https://github.com/Scilence2022/Babel-DNA-Encry/]. It assembles all possible strand sequences when sequencing file(s) and the index range are provided. 

# System Requirements
## Hardware requirements
This package requires only a standard computer with enough RAM to support the k-mer counting.

## Software requirements
### OS Requirements
This package is supported for *Linux*. The package has been tested on the following systems:
+ Ubuntu 20.04.6 LTS
+ GNU Make 4.2.1


# Install and run
```sh
git clone https://github.com/Scilence2022/DBGPS-Babel.git
cd DBGPS-Babel
make

#Check the usage instructions
./DBGPS-greedy-path -h 

Usage: DBGPS-greedy-path [options] <input file> 
                       [Supporting formats: *fq, *fa, *fq.gz, *fa.gz]
Options:
  -k INT     k-mer size [31]
  -i INT     length of index [16] bp
  -l INT     data encoding length [128] bp
  -t INT     number of threads [3]
  -c INT     k-mer coverage cut-off for exclusion of noise k-mers [5]
  -d INT     Switch on k-mer coverage testing mode.  [0]
  -a INT     Initial index [3111111]
  -b INT     End index [3145156]

```

# Example usage

Sample sequencing data is available at: https://doi.org/10.6084/m9.figshare.20424126.v1

```sh
./DBGPS-greedy-path -k 21 ZDNA_1.clean.fq.gz ZDNA_2.clean.fq.gz > assembled.strands
#remove strands with long 'AAAAAAA'
grep -v "AAAAAAA" assembled.strands > assembled.strands.fixAAA
```

The default parameters of this implementation has been optimized for Babel-DNA application. Please specify the actual index length(bp) and the data encoding length(bp) when applying it to your own data. Please use the -i option to specify the length(bp) of the index and use the -a & -b option to specifiy the index range, respectively. 

# License

This project is released under the **GPL License V3**.
