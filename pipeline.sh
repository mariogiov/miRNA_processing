#!/bin/bash

# TODO this probably breaks with gzipped files what with the extension determination
#      maybe easiest just to decompress them first?
# TODO add -f flag to overwrite existing files/directory


function print_usage { echo "Usage: -s <SEQUENCE_FILE> [-f (overwrite existing files)] [-r <GENOME_REFERENCE FILE (fasta)>] [-g <GENOME_FEATURE_FILE (gtf/gff)>] [-w <WORKING_DIRECTORY>] [-m <MIRBASE_FILE>] [-n <CORES>]" >&2 ; }

# GET INPUT
while getopts ":s:r:g:w:m:n:f" opt; do
    case $opt in
        s)
            SEQ_FILE=$OPTARG
            ;;
        r)
            GENOME_REF=$OPTARG
            ;;
        g)
            FEATURES_FILE=$OPTARG
            ;;
        w)
            WORK_DIR=$OPTARG
            ;;
        m)
            MIRBASE_FILE=$OPTARG
            ;;
        n)
            NUM_CORES=$OPTARG
            ;;
        f)
            FORCE_OVERWRITE=1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." 1>&2
            print_usage
            exit 1;
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            print_usage
            exit 1;
            ;;
    esac
done


# CHECK INPUTS

if [[ ! $SEQ_FILE ]]; then
    echo "Must specify input sequence file (-s)." 1>&2
    print_usage
    exit 1;
fi

[[ $GENOME_REF ]] || echo "Warning: no genome reference file (-r) specified; will not perform alignment." 1>&2
[[ $FEATURES_FILE ]] || echo "Warning: no genome feature file (-g) specified; will not perform annotation or feature counting." 1>&2

for file in $SEQ_FILE $GENOME_REF $FEATURES_FILE $MIRBASE_FILE; do
    if [[ ! -e $file ]]; then
        echo "File does not exist: $file" 1>&2
        exit 1
    elif [[ ! -r $file ]]; then
        echo "File cannot be read: $file" 1>&2
        exit 1
    fi
done

SYS_CORES=$(nproc)
if [[ $NUM_CORES =~ ^[0-9]+$ ]]; then
    if [[ $NUM_CORES -gt $SYS_CORES ]]; then
        echo "Number of cores specified ($NUM_CORES) greater than number of cores available ($SYS_CORES). Setting to maximum $SYS_CORES." 1>&2
        NUM_CORES=$SYS_CORES
    fi
else
    echo "Number of cores must be a positive integer between 1 and $SYS_CORES. Setting number of cores to 1." 1>&2
    NUM_CORES=1
fi


if [[ $WORK_DIR ]]; then
    WORK_DIR=$(readlink -f $WORK_DIR)
else
    echo "Info:    no working directory (-d) specified; using '$PWD'" 1>&2
    WORK_DIR=$PWD
fi

SCRIPT_SELF_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# INPUT / DIRECTORY TREE CONSTRUCTION

INPUTFILE=$(basename $SEQ_FILE)
INPUTFILE_ABSPATH=$(readlink -f $SEQ_FILE)
INPUTFILE_DIR=$(dirname $INPUTFILE_ABSPATH)
INPUTFILE_BASE="${INPUTFILE%.*}"
INPUTFILE_EXTENSION="${INPUTFILE##*.}"

# Decompress compressed data automatically
#if [[ $(file $INPUTFILE_ABSPATH | grep gzip) ]]; then
#    echo "Input file is compressed; decompressing to $SEQ_DATA
#fi

if [[ $GENOME_REF ]]; then
    REFERENCE=$(basename $GENOME_REF)
    REFERENCE_ABSPATH=$(readlink -f $GENOME_REF)
    REFERENCE_DIR=$(dirname $REFERENCE_ABSPATH)
    REFERENCE_BASE="${REFERENCE%.*}"
    REFERENCE_EXTENSION="${REFERENCE##*.}"
fi

if [[ $FEATURES_FILE ]]; then
    FEATURES_FILE_ABSPATH=$(readlink -f $FEATURES_FILE)
    FEATURES_FILE=$(basename $FEATURES_FILE)
    FEATURES_FILE_DIR=$(dirname $FEATURES_FILE_ABSPATH)
    FEATURES_FILE_BASE="${FEATURES_FILE%.*}"
    FEATURES_FILE_EXTENSION="${FEATURES_FILE##*.}"
fi

if [[ $MIRBASE_FILE ]]; then
    MIRBASE_ABSPATH=$(readlink -f $MIRBASE_FILE)
    MIRBASE_FILE=$(basename $MIRBASE_FILE)
    MIRBASE_DIR=$(dirname $MIRBASE_ABSPATH)
    MIRBASE_BASE="${MIRBASE_FILE%.*}"
    MIRBASE_EXTENSION="${MIRBASE_FILE##*.}"
fi

# Directories
DATE=$(date "+%Y%m%d_%X")
LOG_DIR=$WORK_DIR"/logs/"
SEQDATA_DIR=$WORK_DIR"/seqdata/"
ALIGNED_DIR=$WORK_DIR"/aligned/"
ANNOTATED_DIR=$WORK_DIR"/annotated/"
VIS_DIR=$WORK_DIR"/visualization/"
LOG_FILE=$LOG_DIR/$INPUTFILE_BASE"_"$DATE".log"
for dir in $LOG_DIR $SEQDATA_DIR $ALIGNED_DIR $ANNOTATED_DIR $VIS_DIR; do
    if [[ $(mkdir -p $dir) -ne 0 ]]; then
        echo "Cannot create directory $dir; exiting." | tee 1>&2
        exit 1
    fi
done


# MODULE LOADING
source $HOME/.virtualenvs/python-276/bin/activate

# Modules, activate the module command
#if [[ $(hostname -s | grep milou) ]]; then
case "$(basename $SHELL)" in
          -sh|sh|*/sh)  modules_shell=sh ;;
       -ksh|ksh|*/ksh)  modules_shell=ksh ;;
       -zsh|zsh|*/zsh)  modules_shell=zsh ;;
    -bash|bash|*/bash)  modules_shell=bash ;;
esac
module() { eval `/usr/local/Modules/$MODULE_VERSION/bin/modulecmd $modules_shell $*`; } 
#fi
export PATH=$HOME/bin:$HOME/.local/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lib/
export PYTHONPATH=$HOME/lib/python2.7/

module unload python
module load python/2.7.4
module load bowtie2/2.1.0
module load cutadapt
module load htseq
module load samtools
#export PATH=/proj/a2010002/nobackup/sw/mf/bioinfo-tools/samtools/0.1.19:$PATH

# Jimmie barcoding experiments
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


echo -e "\n\nADAPTER TRIMMING\n================" | tee -a $LOG_FILE 1>&2
# Remove 3' adapter sequences, discarding reads shorter than $MIN_SEQ_LENGTH
INFILE_CUTADAPT=$INPUTFILE_ABSPATH
OUTFILE_CUTADAPT=$SEQDATA_DIR"/"$INPUTFILE_BASE"_trimmed."$INPUTFILE_EXTENSION
# TruSeq adapter sequence "TGGAATTCTCGGGTGCCAAGG"
ADAPTER="TGGAATTCTCGGGTGCCAAGG"
MIN_SEQ_LENGTH='18'
if [[ ! -f $OUTFILE_CUTADAPT ]] || [[ $FORCE_OVERWRITE ]]; then
    echo -e "Starting cutadapt trimming at $(date)." | tee -a $LOG_FILE 1>&2
    CL="cutadapt -f fastq -a $ADAPTER --match-read-wildcards -O 5 -m $MIN_SEQ_LENGTH --too-short-output=/dev/null -o $OUTFILE_CUTADAPT $INFILE_CUTADAPT"
    echo "Cutadapt command: $CL" | tee -a $LOG_FILE 1>&2
    eval $CL 2>&1 | tee -a $LOG_FILE
else
    echo -e "cutadapt adapter trimming already performed on $(basename $INFILE_CUTADAPT):\noutput file $(basename $OUTFILE_CUTADAPT)\nexists." | tee -a $LOG_FILE 1>&2
fi


echo -e "\n\nALIGNMENT TO GENOME REFERENCE\n=============================" | tee -a $LOG_FILE 1>&2
INFILE_ALN=$OUTFILE_CUTADAPT
OUTFILE_ALN=$ALIGNED_DIR"/"$INPUTFILE_BASE"_aln_"$REFERENCE_BASE".sam"
if [[ ! -f $OUTFILE_ALN ]]  || [[ $FORCE_OVERWRITE ]]; then
    echo -e "Started bowtie2 alignment to genome reference at $(date)." | tee -a $LOG_FILE 1>&2
    # Need some kinda tricky redirection I guess for this to send stderr to a file but stdout to samtools?
    CL="bowtie2 -N 1 -L 18 -p $NUM_CORES -x $REFERENCE_DIR/$REFERENCE_BASE $INFILE_ALN | samtools view -S -b - > $OUTFILE_ALN 2>>$LOG_FILE"
    echo "Executing alignment command: $CL" | tee -a $LOG_FILE
    eval $CL
else
    echo -e "Alignment of $(basename $INFILE_ALN) to $(basename $REFERENCE_ABSPATH) already performed:\n$(basename $OUTFILE_ALN)\nexists." | tee -a $LOG_FILE 1>&2
fi

# Jimmie barcoding experiments
# Move barcode to header
#if [ $BARCODED == 'barcoded' ]
#then
#	SUFFIX=$NEWSUFFIX
#	NEWSUFFIX='_aln_'$REFBASE'_barcoded.sam'
#	OUTFILE_BSAM=$(basename "$OUTFILE_ALN" $SUFFIX)$NEWSUFFIX
#	if [ ! -f $OUTFILE_BSAM ]
#	then
#    	python ../scripts/oop.py 'sam_barcoding' $OUTFILE_ALN $OUTFILE_BSAM
#    	rm $OUTFILE_ALN
#    	cp $OUTFILE_BSAM $OUTFILE_ALN
#	else
#    	echo 'Barcoding of sam already performed'
#    fi
#fi
#

echo -e "\n\nANNOTATION\n==========" | tee -a $LOG_FILE 1>&2
# Convert from BAM to SAM | annotate hits to genes, treating them as nonstranded, mode = union
OUTFILE_ALN_BASE=$(basename "${OUTFILE_ALN%.*}")
OUTFILE_COUNTS=$ANNOTATED_DIR"/"$OUTFILE_ALN_BASE"_counts_"$FEATURES_FILE_BASE".csv"
ANNOTATED_FILE=$ANNOTATED_DIR"/"$OUTFILE_ALN_BASE"_annotated_"$FEATURES_FILE_BASE".sam"
if ( [[ ! -f $OUTFILE_COUNTS ]] || [[ $FORCE_OVERWRITE ]] ) && [[ ! -f $ANNOTATED_FILE ]]; then
    echo -e "Started annotation of alignment at $(date)." | tee -a $LOG_FILE 1>&2
    CL="samtools view $OUTFILE_ALN | htseq-count -o $ANNOTATED_FILE -t exon -s no -q -i 'ID' - $FEATURES_FILE_ABSPATH | sort -n -k 2 -r > $OUTFILE_COUNTS"
    echo -e "Executing command line:\n$CL" | tee -a $LOG_FILE 1>&2
    eval $CL 2>>$LOG_FILE
else
    echo -e "htseq-count already performed: counts file $OUTFILE_COUNTS\nand annotated sam file $ANNOTATED_FILE\nboth exist." | tee -a $LOG_FILE 1>&2
fi


echo -e "\n\nVISUALIZATION\n=============" | tee -a $LOG_FILE 1>&2
# Plot read length distribution with matplotlib
VIS_FILE_NAME=$VIS_DIR/$(basename $OUTFILE_CUTADAPT .fastq)"_readlength-histogram.png"
if [[ ! -f $VIS_FILE_NAME ]] || [[ $FORCE_OVERWRITE ]]; then
    echo -e "Creating read length distribution plot for $OUTFILE_CUTADAPT." | tee -a $LOG_FILE 1>&2
    CL="python $SCRIPT_SELF_DIR/plots.py -i $OUTFILE_CUTADAPT -d $VIS_DIR"
    echo "Executing: $CL" | tee -a $LOG_FILE 1>&2
    eval $CL | tee -a $LOG_FILE 1>&2
else
    echo -e "Visualization already performed: plot file\n$VIS_FILE_NAME\nexists."
fi

# MIRBASE ALIGNMENT
if [[ $MIRBASE_FILE ]]; then
    echo -e "\n\nALIGNMENT TO MIRBASE\n====================" | tee -a $LOG_FILE 1>&2
    INFILE_ALN=$OUTFILE_CUTADAPT
    OUTFILE_ALN=$ALIGNED_DIR"/"$INPUTFILE_BASE"_aln_"$MIRBASE_BASE".sam"
    if [[ ! -f $OUTFILE_ALN ]] || [[ $FORCE_OVERWRITE ]]; then
        echo -e "Started bowtie2 alignment to miRBASE at $(date)." | tee -a $LOG_FILE 1>&2
        # Need some kinda tricky redirection I guess for this to send stderr to a file but stdout to samtools?
        CL="bowtie2 -N 1 -L 18 -p $NUM_CORES -x $MIRBASE_DIR/$MIRBASE_BASE $INFILE_ALN | samtools view -S -b - > $OUTFILE_ALN 2>>$LOG_FILE"
        echo "Executing alignment command: $CL" | tee -a $LOG_FILE
        eval $CL
    else
        echo -e "Alignment of $(basename $INFILE_ALN) to $MIRBASE_FILE already performed:\n$(basename $OUTFILE_ALN)\nexists." | tee -a $LOG_FILE 1>&2
    fi
fi

#if [ ! -f $COUNTFILE ]
#then
#	echo "doesn't exist"
	#samtools view OUTFILE_GENOME | 
#	htseq-count -t 'exon' -s no -q -i 'ID' $ALNFILE $GTFFILE > $COUNTFILE 
	#| sort -n -k 2 -r | tee $COUNTFILE
#fi

# Perform barcoded gene count
#if [ $BARCODED == 'barcoded' ]
#then
#	SUFFIX=$NEWSUFFIX
#	NEWSUFFIX='_count_'$REFBASE'_barcoded.csv'
#	OUTFILE_BSAM=$(basename "$COUNTFILE" $SUFFIX)$NEWSUFFIX
#	if [ ! -f $OUTFILE_BSAM ]
#	then
#    	python ../scripts/oop.py 'sam_annotation' $COUNTFILE $OUTFILE_BSAM $GTFFILE
#	else
#    	echo 'Barcoding of sam already performed'
#    fi
#fi
#


# Plot Size distribution of each library
#ALNFILE
#GTFFILE
#COUNTFILE

#python initialplots.py >
