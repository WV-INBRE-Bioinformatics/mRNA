# mRNA
Scripts for analysis of RNA-Seq, specifically of messenger RNA

---

align_mRNA.sh
=============

Aligns paired-end reads of mRNA samples stored in fastq.gz files, producing bam files as output, using [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) version 2.0.4.

Assumptions:
------------

1. `source_dir` contains a collection of pairs of fastq.gz files, with samples represented by a "read 1" and "read 2" files. "Read 1" files end with `_R1.fastq.gz` and the matching "read 2" file is identical except that it ends `_R2.fastq.gz`. 
2. Sample names are represented by the "read 1" file name with `_R1.fastq.gz` stripped from the end. Output files will be the sample name with the extension `.bam` added.

Requirements:
-------------

1. [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) version 2.0.4.
2. [samtools](http://www.htslib.org/) version 1.3.1.

Syntax:
-------

    align.sh [options] source_dir genome_file output_dir
    
    
