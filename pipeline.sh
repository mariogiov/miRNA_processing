#!/bin/bash

# USAGE
# <scriptname> input_fastq_file REFERENCE_GENOME ANNOTATION_FILE
function print_usage { echo "Usage: -s <SEQUENCE_FILE> -r <GENOME_REFERENCE FILE (fasta)> -g <GENOME_FEATURE_FILE (gtf/gff)> -d <WORKING_DIRECTORY>." >&2 ; }


#
#Create a project directory containing a script directory with all scripts and a reference directory with bowtie refs and gff file.
# as well as a seqdata directory with the raw sequencing data.
#For the moment bowtie has to be under bowtie-1.0.0/ in the project directory.
#


# GET INPUT
while getopts ":s:r:g:d:" opt; do
    case $opt in
        s)
            SEQ_FILE=$OPTARG
            ;;
        r)
            GENOME_REF=$OPTARG
            ;;
        g)
            GFF_FILE=$OPTARG
            ;;
        d)
            WORK_DIR=$OPTARG
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            print_usage
            exit 1;
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            print_usage
            exit 1;
            ;;
    esac
done


# CHECK INPUTS

if [[ ! $SEQ_FILE ]]; then
    echo "Must specify input sequence file (-s)."
    print_usage
    exit 1;
fi

if [[ $GENOME_REF ]] || [[ $GFF_FILE ]]; then
    [[ $GENOME_REF ]] || echo "Warning: no genome reference file (-r) specified; will not perform alignment."
    [[ $GFF_FILE ]] || echo "Warning: no genome feature file (-g) specified; will not perform annotation or feature counting."
else
    echo "No genome reference file (-r) or genome feature file (-g) specified; nothing to do. Exiting."
    print_usage
    exit
fi

for file in $SEQ_FILE $GENOME_REF $GFF_FILE; do
    if [[ ! -e $file ]]; then
        echo "File does not exist: $file"
        exit
    elif [[ ! -r $file ]]; then
        echo "File cannot be read: $file"
        exit
    fi
done

if [[ ! $WORK_DIR ]]; then
    echo "Info:    no working directory (-d) specified; using '$PWD'"
    WORK_DIR=$PWD
fi
echo "WORK_DIR: $WORK_DIR"

# INPUT / DIRECTORY TREE CONSTRUCTION

DATE=$(date "+%Y%m%d")

INPUTFILE=$(basename $SEQ_FILE)
INPUTFILE_ABSPATH=$(readlink -f $SEQ_FILE)
INPUTFILE_DIR=$(dirname $INPUTFILE_ABSPATH)
INPUTFILE_BASE="${INPUTFILE%.*}"
INPUTFILE_EXTENSION="${INPUTFILE##*.}"

#echo "INPUTFILE:            $INPUTFILE"
#echo "INPUTFILE_ABSPATH:    $INPUTFILE_ABSPATH"
#echo "INPUTFILE_DIR:        $INPUTFILE_DIR"
#echo "INPUTFILE_BASE:       $INPUTFILE_BASE"
#echo "INPUTFILE_EXTENSION:  $INPUTFILE_EXTENSION"

if [[ $GENOME_REF ]]; then
    REFERENCE=$(basename $GENOME_REF)
    REFERENCE_ABSPATH=$(readlink -f $GENOME_REF)
    REFERENCE_DIR=$(dirname $REFERENCE_ABSPATH)
    REFERENCE_BASE="${REFERENCE%.*}"
    REFERENCE_EXTENSION="${REFERENCE##*.}"
fi

#echo "REFERENCE:            $REFERENCE"
#echo "REFERENCE_ABSPATH:    $REFERENCE_ABSPATH"
#echo "REFERENCE_DIR:        $REFERENCE_DIR"
#echo "REFERENCE_BASE:       $REFERENCE_BASE"
#echo "REFERENCE_EXTENSION:  $REFERENCE_EXTENSION"


if [[ $GFF_FILE ]]; then
    GFF_FILE=$(basename $GFF_FILE)
    GFF_FILE_ABSPATH=$(readlink -f $GFF_REF)
    GFF_FILE_DIR=$(dirname $GFF_FILE_ABSPATH)
    GFF_FILE_BASE="${GFF_FILE%.*}"
    GFF_FILE_EXTENSION="${GFF_FILE##*.}"
fi

#echo "GFF_FILE:            $GFF_FILE"
#echo "GFF_FILE_ABSPATH:    $GFF_FILE_ABSPATH"
#echo "GFF_FILE_DIR:        $GFF_FILE_DIR"
#echo "GFF_FILE_BASE:       $GFF_FILE_BASE"
#echo "GFF_FILE_EXTENSION:  $GFF_FILE_EXTENSION"

# Directories
LOG_DIR=$WORK_DIR"/logs/"
SEQDATA_DIR=$WORK_DIR"/seqdata/"
ALIGNED_DIR=$WORK_DIR"/aligned/"
ANNOTATED_DIR=$WORK_DIR"/annotated/"
LOGFILE=$LOG_DIR/$INPUFILET_BASE$DATE".log"
for dir in $LOG_DIR $SEQDATA_DIR $ALIGNED_DIR $ANNOTATED_DIR; do
    if [[ -e $dir ]]; then
        # TODO add -f flag to force overwrite
        echo "Directory '$dir' exists. Exiting."
        exit
    elif [[ $(mkdir -p $dir) -ne 0 ]]; then
        echo "Cannot create directory $dir; exiting."
        exit
    fi
done

#echo "LOG_DIR:          $LOG_DIR"
#echo "SEQDATA_DIR:      $SEQDATA_DIR"
#echo "ALIGNED_DIR:      $ALIGNED_DIR"
#echo "ANNOTATED_DIR:    $ANNOTATED_DIR"
#echo "LOGFILE:          $LOGFILE"

exit

# MODULE LOADING
# Modules, activate the module command
case "$0" in
          -sh|sh|*/sh)  modules_shell=sh ;;
       -ksh|ksh|*/ksh)  modules_shell=ksh ;;
       -zsh|zsh|*/zsh)  modules_shell=zsh ;;
    -bash|bash|*/bash)  modules_shell=bash ;;
esac
module() { eval `/usr/local/Modules/$MODULE_VERSION/bin/modulecmd $modules_shell $*`; } 
module load cutadapt
module load bowtie2

# I'm not sure what this code is for -- I'm guessing if your samples aren't demultiplexed?
# Move barcode to header
#if [ $BARCODED == 'barcoded' ]
#then
#	NEWSUFFIX='_barcoded.fastq'
#	OUTFILE_BARCODING=$(basename "$INFILE" $SUFFIX)$NEWSUFFIX
#	if [ ! -f $OUTFILE_BARCODING ]
#	then
#    	python ../scripts/oop.py 'fq_barcoding' $INFILE $OUTFILE_BARCODING
#	else
#    	echo 'Barcoding already performed'
#    fi
#    INFILE=$OUTFILE_BARCODING
#fi
#

# ADAPTER TRIMMING
# Remove 3' adapter sequences, discarding reads shorter than $LENGTH
PRESUFFIX='_trimmed'
INFILE_CUTADAPT=$INPUTFILE_ABSPATH
OUTFILE_CUTADAPT=$INPUTFILE_DIR$INPUTFILE_BASE$PRESUFFIX$INPUTFILE_EXTENSION
# TruSeq adapter sequence "TGGAATTCTCGGGTGCCAAGG"
ADAPTER="TGGAATTCTCGGGTGCCAAGG"
LENGTH='18'

if [ ! -f $OUTFILE_CUTADAPT ]
then
    cutadapt -a $ADAPTER --match-read-wildcards -O 5 -m $LENGTH --too-short-output=/dev/null -o $OUTFILE_CUTADAPT $INFILE_CUTADAPT | tee $LOGFILE
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


  
