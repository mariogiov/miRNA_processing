#!/bin/bash  
#
#Create a project directory containing a script directory with all scripts and a reference directory with bowtie refs and gff file.
# as well as a seqdata directory with the raw sequencing data.
#For the moment bowtie has to be under bowtie-1.0.0/ in the project directory.
#
# Take Input and contruct directory structure        
BASE=$PWD
DATE=$(date "+%Y%m%d")
INFILE=$1
SUFFIX=$2
REFBASE=$3
BARCODED=$4
NEWSUFFIX=$2
LOGFILE='../logs/'$(basename "$INFILE" $SUFFIX)"_"$DATE".log"
#
#  project/seqdata, aligned, annotated, scripts, logs
#
mkdir ../logs ../seqdata ../aligned ../annotated ../references #, ../plots

cd ../seqdata
# Move barcode to header
if [ $BARCODED == 'barcoded' ]
then
	NEWSUFFIX='_barcoded.fastq'
	OUTFILE_BARCODING=$(basename "$INFILE" $SUFFIX)$NEWSUFFIX
	if [ ! -f $OUTFILE_BARCODING ]
	then
    	python ../scripts/oop.py 'fq_barcoding' $INFILE $OUTFILE_BARCODING
	else
    	echo 'Barcoding already performed'
    fi
    INFILE=$OUTFILE_BARCODING
fi

# Remove adaptor and throw away too short reads
SUFFIX=$NEWSUFFIX
NEWSUFFIX='_trimmed.fastq'
INFILE_CUTADAPT=$INFILE
OUTFILE_CUTADAPT=$(basename "$INFILE_CUTADAPT" $SUFFIX)$NEWSUFFIX
ADAPTOR="TGGAATTCTCGGGTGCCAAGG"
#Think it's the TruSeq adapter   "TGGAATTCTCGGGTGCCAAGG"
LENGTH='18'

if [ ! -f $OUTFILE_CUTADAPT ]
then
    cutadapt -a $ADAPTOR --match-read-wildcards -O 5 -m $LENGTH --too-short-output=/dev/null -o $OUTFILE_CUTADAPT $INFILE_CUTADAPT | tee $LOGFILE
else
	echo 'Trimming already performed.'
fi

cd ../aligned
#Align to genome
#seed = 18, 1 missmatch, keep best hit
#convert from SAM to BAM
SUFFIX=$NEWSUFFIX
NEWSUFFIX='_aln_'$REFBASE'.sam'
REFERENCE_GENOME='../references/'$GREFBASE
REFERENCE='../references/'$REFBASE
INFILE_ALN='../seqdata/'$OUTFILE_CUTADAPT
OUTFILE_ALN=$(basename "$INFILE_ALN" $SUFFIX)$NEWSUFFIX

if [ ! -f $OUTFILE_ALN ]
then
    #../bowtie-1.0.0/  needed for testing localy
    bowtie -n 1 -l 18 -k 1 --best -p 8 $REFERENCE $INFILE_ALN -S > $OUTFILE_ALN #| tee $LOGFILE #| samtools view -S -b - > $OUTFILE_MIRBASE 
else
	echo 'Alignment to mirbase already performed.'
fi

# When it's working!!!!!
#rm $OUTFILE_CUTADAPT

# Move barcode to header
if [ $BARCODED == 'barcoded' ]
then
	SUFFIX=$NEWSUFFIX
	NEWSUFFIX='_aln_'$REFBASE'_barcoded.sam'
	OUTFILE_BSAM=$(basename "$OUTFILE_ALN" $SUFFIX)$NEWSUFFIX
	if [ ! -f $OUTFILE_BSAM ]
	then
    	python ../scripts/oop.py 'sam_barcoding' $OUTFILE_ALN $OUTFILE_BSAM
    	rm $OUTFILE_ALN
    	cp $OUTFILE_BSAM $OUTFILE_ALN
	else
    	echo 'Barcoding of sam already performed'
    fi
fi

cd ../annotated
# Convert from BAM to SAM | annotate hits to genes, treating them as nonstranded, mode = union
SUFFIX=$NEWSUFFIX
NEWSUFFIX='_count_'$REFBASE'.csv'
ALNFILE='../aligned/'$OUTFILE_ALN
GTFFILE='../references/'$REFBASE'.gff'
COUNTFILE=$(basename "$ALNFILE" '_aln_'$REFBASE'.sam')$NEWSUFFIX

if [ ! -f $COUNTFILE ]
then
	echo "doesn't exist"
	#samtools view OUTFILE_GENOME | 
	htseq-count -t 'exon' -s no -q -i 'ID' $ALNFILE $GTFFILE > $COUNTFILE 
	#| sort -n -k 2 -r | tee $COUNTFILE
fi

# Perform barcoded gene count
if [ $BARCODED == 'barcoded' ]
then
	SUFFIX=$NEWSUFFIX
	NEWSUFFIX='_count_'$REFBASE'_barcoded.csv'
	OUTFILE_BSAM=$(basename "$COUNTFILE" $SUFFIX)$NEWSUFFIX
	if [ ! -f $OUTFILE_BSAM ]
	then
    	python ../scripts/oop.py 'sam_annotation' $COUNTFILE $OUTFILE_BSAM $GTFFILE
	else
    	echo 'Barcoding of sam already performed'
    fi
fi



# Plot Size distribution of each library
#ALNFILE
#GTFFILE
#COUNTFILE

#python initialplots.py > 


  
