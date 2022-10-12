# Brief overview
A pipeline for mapping sequencing data to viral genome reference and visualization of depth of coverage.
GenMapViz is a bioinformatics pipeline that I use to align sequencing reads (fastq files) to a reference genome outpus visual representation of coverage depth. I made this pipeline with intention of working with viral genomes and is therefore optimized for small genomes. 

# Installation 

The pipeline is meant to be used an Unix/Linux operative system.  The pipiline does not need tobe installed  does not need to be installed. All that is needed is to have the main python script and the scripts directory in the same location.  However, there are several dependancies that need to be pre-installed for it to work. 

# Dependancies

## Python and Biopython

The best way to install most of the tools is to install [anaconda](https://www.anaconda.com/) and use `conda` and/or `pip` to install specific tools.

* [Python 3.x](https://www.python.org/downloads/)
* [Biopython](https://biopython.org/wiki/Download)

## Tools for alignment and variant calling

* [Burrows-Wheeler Aligner (bwa)](http://bio-bwa.sourceforge.net/) : Alignment of sequencing reads to a reference genome.
* [fastp](https://github.com/OpenGene/fastp) : Quality Control - Trimming out low quality reads and remove adatpers.
* [SAMtools](http://www.htslib.org/) : Manipulating alignments data.
* [bcftools](http://www.htslib.org/) : Calling variants. Should be installed together with samtools.
* [bedtools](https://github.com/arq5x/bedtools2) : Controlling interval data  like coordinates of regions on the genome.


## R packages

These are better installed and managed by [Bioconductor](https://www.bioconductor.org/), a software project for management and analysis of genomic and other high throughput data.

* [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html) : Vizualization of genomic data.
* [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) : Handling genomics locations like genes and other features to be annotated.
* [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) : Storage and manipulation of genomic intervals and variables defined along a genome.
* [data.table](https://www.rdocumentation.org/packages/data.table/versions/1.14.2) : Can be installed by `install.packages("data.table")`.
* [tools]() : Can be installed by `install.packages("tools")`.

__Note__: As I was installing Gviz, the process failed several times because of uninstalled depandancies. I had to read the Bioconductor output, find all depandancies, installed them and restarted the process to have it finally installed. Pacience was the key.


# Usage

```
usage: GenMapViz.py [-h] -s SAMPLE [-v VIRUS] [-1 FASTQ1] [-2 FASTQ2] [-r REFFILE] [-g GFF] -o OUTDIR
                    [-a ALIGNER] [-t THREADS] [-qc] [-zs ZOOM_START] [-ze ZOOM_END]

options:
  -h, --help            show this help message and exit
  -s SAMPLE, --sample SAMPLE
                        Sample Name. Will be the base name for result files.
  -v VIRUS, --virus VIRUS
                        Name of the virus
  -1 FASTQ1, --fastq1 FASTQ1
                        First file of raw reads.
  -2 FASTQ2, --fastq2 FASTQ2
                        Second file of raw reads.
  -r REFFILE, --refFile REFFILE
                        Path to the reference genome file (in fasta format).
  -g GFF, --gff GFF     Path to the gene file (in GFF3 format). Can also be a tab delimited file with 3
                        columns: start, end, symbol.
  -o OUTDIR, --outdir OUTDIR
                        Path to output directory.
  -a ALIGNER, --aligner ALIGNER
                        Short read aligner to be used. options = bwa, bowtie2, bbmap. Default: bwa (bwa
                        mem will be used).
  -t THREADS, --threads THREADS
                        Number of threads. Default: 1 or half of available cpus.
  -qc, --qualCntr       Perfom quality control prior to alignment (fastp will be used).
  -zs ZOOM_START, --zoom_start ZOOM_START
                        The start of the area to be zoomed. Default: 1/3 of the reference length.
  -ze ZOOM_END, --zoom_end ZOOM_END
                        The end of the area to be zoomed. Default: 2/3 of the reference length.
```

## To do:

* Add variant call option
* Add variant visualization to the coverage graph

