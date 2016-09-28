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
    
---    

run_de.Rscript
==============

Performs gene-level differential expression analysis starting from a set of bam files generated from alignment of RNA-Seq reads of mRNA.

This uses R version 3.3.1 and the bioconductor packages DESeq2, GenomicAlignments, GenomicFeatures, AnnotationDbi, GenomicRanges,  Rsamtools, pheatmap, RColorBrewer, SummarizedExperiment, and BiocParallel along with their dependencies Biobase, Biostrings, XVector, GenomeInfoDb, IRanges, S4Vectors, and BiocGenerics. These packages may be installed with 

    source("http://bioconductor.org/biocLite.R")
    biocLite("Rsamtools")
    biocLite("GenomicFeatures")
    biocLite("GenomicAlignments")
    biocLite("BiocParallel")
    biocLite("DESeq2")
    biocLite("pheatmap")
    biocLite("RColorBrewer")
    biocLite("AnnotationDbi")
    
Additionally, the annotation package for the organism to which the reads are aligned is required. E.g. for human data, the package org.Hs.eg.db is required, for mouse, org.Mm.eg.db, etc. These can similarly be installed with

    biocLite("org.Hs.eg.db")
    
etc. 

Assumptions:
------------

This script loads a collection of specifically-named bam files from a directory, counts gene hits according to a specified 
annotation file in gtf format, and performs differential expression analysis using DESeq2 and a specified, essentially arbitrary
statistical model.

Each sample must be represented by a single bam file, and these bam files must be named according to the values ascribed to model variables for the corresponding sample. By default values are separated by `-` in the file name, though this separator can be changed 
via the `-s` or `--filenameSplit` option. After removing the `.bam` extension from filename and splitting the remainder of the filenames by this delimiter, all files should resolve to the same number of tokens, say `p`. There should be `p` occurences of the `-a` or `--attribute` option to the script invocation, with each one naming the variable whose value is specified by the filename, in order. The `-m` or `--model` option specifies a statistical model in R syntax, and must only include variables defined by an `-a` or `--attribute` option (though it is not necessary to include all variables).

Examples:

For a simple one-variable experiment, the files might be named

 - control-1.bam
 - control-2.bam
 - control-3.bam
 - control-4.bam
 - treated-1.bam
 - treated-2.bam
 - treated-3.bam
 - treated-4.bam

and then the script would include the options

    -a treatment -a replicate -m "~ treatment"
    
The following is an example of a two-variable experiment, with treatment and tissue type the variables, in which we were interested in the interaction term:

 - control-normal-1.bam
 - control-normal-2.bam
 - control-normal-3.bam
 - control-normal-4.bam
 - control-tumor-1.bam
 - control-tumor-2.bam
 - control-tumor-3.bam
 - control-tumor-4.bam
 - treated-normal-1.bam
 - treated-normal-2.bam
 - treated-normal-3.bam
 - treated-normal-4.bam
 - treated-tumor-1.bam
 - treated-tumor-2.bam
 - treated-tumor-3.bam
 - treated-tumor-4.bam

and then the options would include

    -a treatment -a tissue -a replicate -m "~ treatment + tissue + treatment:tissue"
    
The script requires a single directory containing the bam files (specified by `-b` or `--bamDirectory`), a GTF-format file containing
annotations corresponding to the genome to which the reads are aligned (specified by `-g` or `--gtfFile`), an output file to save the R data image, attributes and model as described above, and a file to which to save significantly differentially expressed genes in CSV format (`-v` or `--csvFile`). Other options are described in the "syntax" section below.
    
Requirements:
-------------

Requires R version 3.3.1 and bioconductor packages described above. The output from `sessionInfo()` from a typical run is shown below:

    R version 3.3.1 (2016-06-21)
    Platform: x86_64-redhat-linux-gnu (64-bit)
    Running under: Scientific Linux 7.2 (Nitrogen)

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] parallel  stats4    methods   stats     graphics  grDevices utils    
    [8] datasets  base     

    other attached packages:
     [1] org.Hs.eg.db_3.3.0         RColorBrewer_1.1-2        
     [3] pheatmap_1.0.8             DESeq2_1.12.4             
     [5] BiocParallel_1.6.6         GenomicAlignments_1.8.4   
     [7] SummarizedExperiment_1.2.3 GenomicFeatures_1.24.5    
     [9] AnnotationDbi_1.34.4       Biobase_2.32.0            
    [11] Rsamtools_1.24.0           Biostrings_2.40.2         
    [13] XVector_0.12.1             GenomicRanges_1.24.3      
    [15] GenomeInfoDb_1.8.7         IRanges_2.6.1             
    [17] S4Vectors_0.10.3           BiocGenerics_0.18.0       

    loaded via a namespace (and not attached):
     [1] Rcpp_0.12.7         plyr_1.8.4          bitops_1.0-6       
     [4] tools_3.3.1         zlibbioc_1.18.0     digest_0.6.10      
     [7] rpart_4.1-10        biomaRt_2.28.0      annotate_1.50.0    
    [10] RSQLite_1.0.0       gtable_0.2.0        lattice_0.20-34    
    [13] Matrix_1.2-7.1      DBI_0.5-1           gridExtra_2.2.1    
    [16] genefilter_1.54.2   rtracklayer_1.32.2  cluster_2.0.4      
    [19] locfit_1.5-9.1      nnet_7.3-12         grid_3.3.1         
    [22] data.table_1.9.6    XML_3.98-1.4        survival_2.39-5    
    [25] foreign_0.8-67      latticeExtra_0.6-28 Formula_1.2-1      
    [28] geneplotter_1.50.0  ggplot2_2.1.0       Hmisc_3.17-4       
    [31] scales_0.4.0        splines_3.3.1       xtable_1.8-2       
    [34] colorspace_1.2-6    labeling_0.3        acepack_1.3-3.3    
    [37] RCurl_1.95-4.8      munsell_0.4.3       chron_2.3-47   


Syntax:
-------

    run_de.Rscript [options]
    
Options:

 - -b --bamDirectory Directory of input bam files (required)
 - -g --gtffile gtf file of gene annotations for organism to which reads are aligned (required)
 - -O --outputFile file to which R data image is saves (required)
 - -m --model R formula for statistical model for testing (required)
 - -a --attribute name of variable that may appear in model, values are determined from filenames (required for each variable)
 - -s --filenameSplit delimiter for splitting filenames to determine variables (default `'-'`)
 - -o --organism organism specifier (Hs for human, Mm for mouse, etc) (required)
 - -c --sampleClusterFilename File name for sample cluster map image (default sampleCluster.png)
 - -p --samplePCAFilename File name for sample PCA plot (default samplePCA.png)
 - -n --numHeatmapGenes Number of genes used for heatmap (default 250)
 - -h --geneHeatmapFilename File name for heatmap image (default geneHeatmap.png)
 - -H --geneHeatmapHeight Height (in pixels) of heatmap (default 1920)
 - -v --csvFile File name for CSV file of significant genes and associated data (required)
 - ?  --help Display help
 
