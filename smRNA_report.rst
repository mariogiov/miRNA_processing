.. image:: /home/mario/miRNA_processing/sll_logo.jpg
   :height: 429
   :width: 1500

smRNA Analysis Report
=====================

SciLifeLab Stockholm
--------------------

Date
~~~~
2013-12-17

Project Name
~~~~~~~~~~~~
C.Dixelius_13_01

UPPNEX Project ID
~~~~~~~~~~~~~~~~~
b2013216

Samples
~~~~~~~
15 samples total

    - dicerm
    - pig 0
    - 88069m
    - dicerd24
    - dicerd48
    - dicerd72
    - dicerp24
    - dicerp48
    - dicerp72
    - 88069p24
    - 88069p48
    - 88069p72
    - 88069D24
    - 88069D48
    - 88069D72

|

The miRNA Pipeline Step-by-Step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Your samples were analyzed as follows:

1. **Merging:** If samples in your project were split across multiple lanes, flowcells, or runs, the fastq files were merged before processing.
2. **Trimming:** Every read was scanned for adapter sequences corresponding to the adapters used during library preparation, and if found, the adapter sequence and everything downsteam of it were clipped off. Additionally, reads were clipped at any base with a quality score below Q10.
3. **Reference Alignment:** Trimmed reads were aligned to the reference file 'Phytophthora_infestans.ASM14294v1.20.dna.toplevel.fa' and the reference file 'PGSC_DM_v4-03_pseudomolecules' using bowtie2.
4. **Gene Counting:** The number of times each read aligned to a feature in the annotation file 'asm1429v1.gff' was counted using htseq-count.
5. **miRBase Alignment:** Trimmed reads were aligned to the microRNA database (see http://http://www.mirbase.org/).
6. **Visualization:** A distribution of trimmed read lengths was plotted.

|

Directory Structure
~~~~~~~~~~~~~~~~~~~
Processing a project produces three directories for each sample:

- aligned
    - *<sample_name>_aln_miRBase.bam*: an alignment file in BAM format of sample reads aligned to miRBase.
    - *<sample_name>_aln_PGSC_DM_v4-03_pseudomolecules.bam*: an alignment file in BAM format of sample reads aligned to the draft potato genome assembly.

- annotated
    - *<sample_name>_aln_Phytophthora_infestans<...>_annotated_asm14294v1.bam*: an alignment file in BAM format of sample reads aligned to the Phytophthora_infestans genome annotated with feature information about the site to which they mapped.

    - *<sample_name>_aln_Phytophthora_infestans<...>_annotated.csv*: a csv file listing each gene feature from the annotation file and how many reads mapped to it.

- logs
    - *<sample_name>_<date-time>.log*: a text file with output from the sample processing showing the commands executed and their results.

- visualization
    - *<sample_name>_readlength-histogram.png*: a JPG-formatted histogram of the trimmed read lengths.

|

Also included are the reference files (FASTA, GFF) used for alignment and annotation.

|

Sample Alignment Statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Alignment to reference 'Phytophthora_infestans.ASM14294v1.20.dna.toplevel'*
_______________________________________________________________________________

:dicerm:
	 93.76%
:pig 0:
	 12.84%
:88069m:
	 95.11%
:dicerd24:
	 17.58%
:dicerd48:
	 16.22%
:dicerd72:
	 7.95%
:dicerp24:
	 15.67%
:dicerp48:
	 8.26%
:dicerp72:
	 13.09%
:88069p24:
	 9.07%
:88069p48:
	 15.59%
:88069p72:
	 12.63%
:88069D24:
	 15.18%
:88069D48:
	 17.29%
:88069D72:
	 28.30%


*Alignment to reference 'PGSC_DM_v4-03_pseudomolecules'* (Potato draft assembly)
_________________________________________________________________________________

:dicerm:
    14.52%
:pig 0:
    97.92%
:88069m:
    8.20%
:dicerd24:
    91.15%
:dicerd48:
    94.71%
:dicerd72:
    96.84%
:dicerp24:
    95.58%
:dicerp48:
    97.24%
:dicerp72:
    95.72%
:88069p24:
    97.50%
:88069p48:
    94.62%
:88069p72:
    94.37%
:88069D24:
    92.79%
:88069D48:
    92.65%
:88069D72:
    81.22%


*Alignment to reference 'miRBase'*
__________________________________

:dicerm:
    0.22%
:pig 0:
    2.66%
:88069m:
    0.10%
:dicerd24:
    2.12%
:dicerd48:
    2.49%
:dicerd72:
    2.60%
:dicerp24:
    2.00%
:dicerp48:
    2.32%
:dicerp72:
    2.30%
:88069p24:
    2.11%
:88069p48:
    2.13%
:88069p72:
    2.44%
:88069D24:
    2.32%
:88069D48:
    2.12%
:88069D72:
    2.26%

|

Sample Read Length Histograms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*01_dicerm*
-----------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_101_rpi1/visualization/1_131014_BH1211ADXX_P570_101_rpi1_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*02_pig_0*
----------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_102F_rpi2/visualization/1_131014_BH1211ADXX_P570_102F_rpi2_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*03_88069m*
-----------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_103F_rpi3/visualization/1_131014_BH1211ADXX_P570_103F_rpi3_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*04_dicerd24*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_104F_rpi4/visualization/1_131014_BH1211ADXX_P570_104F_rpi4_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*05_dicerd48*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_105F_rpi5/visualization/1_131014_BH1211ADXX_P570_105F_rpi5_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*06_dicerd72*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_106F_rpi6/visualization/1_131014_BH1211ADXX_P570_106F_rpi6_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*07_dicerp24*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_107F_rpi7/visualization/1_131014_BH1211ADXX_P570_107F_rpi7_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*08_dicerp48*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_108F_rpi8/visualization/1_131014_BH1211ADXX_P570_108F_rpi8_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*09_dicerp72*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_109F_rpi9/visualization/1_131014_BH1211ADXX_P570_109F_rpi9_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*10_88069p24*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_110F_rpi10/visualization/1_131014_BH1211ADXX_P570_110F_rpi10_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*11_88069p48*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_111F_rpi11/visualization/1_131014_BH1211ADXX_P570_111F_rpi11_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*12_88069p72*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_112F_rpi12/visualization/1_131014_BH1211ADXX_P570_112F_rpi12_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*13_88069D24*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_113F_rpi13/visualization/1_131014_BH1211ADXX_P570_113F_rpi13_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*14_88069D48*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_114F_rpi14/visualization/1_131014_BH1211ADXX_P570_114F_rpi14_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


*15_88069D72*
-------------
.. image:: /pica/v1/b2013064/INBOX/C.Dixelius_13_01/processed/P570_115F_rpi15/visualization/1_131014_BH1211ADXX_P570_115F_rpi15_1_merged_trimmed_readlength-histogram.jpg
   :width: 5in
   :height: 4in


