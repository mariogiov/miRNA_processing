smRNA / miRNA Processing
========================

Usage
-----

Usage:

    pipeline.sh
                [-r <genome_reference_file> (FASTA)>]
                [-g <genome_feature_file (GTF/GFF)>]
                [-m <mirbase_file (FASTA)>]
                [-o <output_directory>]
                [-n <cores>]
                [-f (overwrite existing files)]
                [-k (keep temp files)]
                <sequence_file> [<additional_sequence_files> <will_be_merged> <before_processing>]


What It Does
------------
Aligns single-read (not paired-end) RNA samples to reference genomes, annotates them, and counts gene frequencies. Also aligns against miRBase. Visualizes read lengths after trimming.
More specifically:

 - Merges multiple FASTQ files into a single FASTQ file (simple concatenation)
 - Trims adapter sequences using cutadapt
 - Aligns to a reference file specified by the user using bowtie2
 - Tallies the number of times reads align to each feature in an annotation file specified by the user
 - Aligns to the miRBase microRNA database (http://www.mirbase.org)
 - Visualizes the read lengths after trimming on a histogram created with matplotlib


What It Does Not Do
-------------------
 - Automatically generate a report or statistics. This must be done manually, although an example report in RST format is supplied


Required Software
-----------------
 - bowtie2
 - HTSeq (python module)
 - matplotlib (python module)


Input
-----
 - Single-read sample data in FASTQ format
 - Reference genome for alignment in fasta format (with bowtie2 indexes)
 - Annotation file in gff format
 - miRBase file for alignment in fasta format (with bowtie2 indexes)


Output
------
 - Alignment files in BAM format
 - Annotated alignment files in BAM format
 - Feature frequency counts in CSV format
 - Visualization of read lengths in JPG format
 - A log of the commands executed and their output


Example Usage
-------------
::

    bash pipeline.sh -r path/to/reference_file.fa -g path/to/annotation_file.gff -m path/to/miRBase.fa -o output_dir/ -n 8 -f -k data_file1.fastq data_file2.fq data_file3_fastq

