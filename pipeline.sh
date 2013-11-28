#!/bin/bash

# TODO the concatenation of input fastq files to merge will cause the python script to break if there is a frameshift e.g. from extra headers
# TODO it would be good to log the file rigamarole to know when decompression fails or if files are not readable, etc.
# TODO consider checking if target files for decompression already exist? currently they are overwritten

# DEFINE FUNCTIONS

# print_usage()
function print_usage { echo -e "\nUsage:\t$0\n\t\t[-r <genome_reference_file> (FASTA)>]\n\t\t[-g <genome_feature_file (GTF/GFF)>]\n\t\t[-m <mirbase_file (FASTA)>]\n\t\t[-w <working_directory>]\n\t\t[-f (overwrite existing files)]\n\t\t[-n <cores>]\n\t\t<sequence_file> [<additional_sequence_files> <will_be_merged> <before_processing>]" >&2 ; }

# extension_is_fastq()
NOT_FASTQ_ERROR_TEXT="input file is not in fastq format (doesn't end with .fastq or .fq)"
function extension_is_fastq () {
    TIF_EXT="${1##*.}"
    if ( [[ $TIF_EXT == "fq" ]] || [[ $TIF_EXT == "fastq" ]] ); then
        echo $TIF_EXT
        return 0
    else
        return 1
    fi
}

# exists_is_readable()
FILE_NOT_EXISTS_OR_NOT_READABLE_ERROR_TEXT="file does not exist or is not readable"
function exists_is_readable () {
#for file in $SEQ_FILE $GENOME_REF $FEATURES_FILE $MIRBASE_FILE; do
    if [[ ! -e $1 ]]; then
        #echo "Error: file does not exist: $1" 1>&2
        return 1
    elif [[ ! -r $1 ]]; then
        #echo "Error: file cannot be read: $1" 1>&2
        return 1
    else
        echo "$1"
        return 0
    fi
#done
}

# decompress_file()
function decompress_file () {
    if [[ $(file $1 | grep gzip ) ]]; then
        # If file is compressed
        echo -e "Info:\t\tinput file \"$1\" is compressed; decompressing..." 1>&2
        DF_OUTPUT_NAME=$(basename "${1%.*}")
        if [[ $2 ]]; then
            # If output directory is passed, decompress to alternate output directory
            if [[ ! $(mkdir -p $2) -eq 0 ]]; then
                echo -e "Warning:\t\toutput directory \"$2\" could not be created. Decompressing to PWD \"$PWD\"." 1>&2
                OUTPUT_DIR=$PWD
            else
                OUTPUT_DIR=$(readlink -f $2)
            fi
            echo -en "Info:\t\tgunzipping file \"$1\" to \"$OUTPUT_DIR/$DF_OUTPUT_NAME\"..." >&2
            gunzip -c $1 >> $OUTPUT_DIR/$DF_OUTPUT_NAME
            if [[ ! $? -eq 0 ]]; then
                echo -e "Error:\t\tcouldn't decompress file $1 to output directory \"$2\"." 1>&2
                return 1
            else
                echo " done." 1>&2
            fi
            echo $OUTPUT_DIR/$DF_OUTPUT_NAME
        else
            # Else unzip in place
            gzip -d $1
            if [[ ! $? -eq 0 ]]; then
                echo -e "Error:\t\tcouldn't decompress file $1." 1>&2
                return 1
            fi
            echo $DF_OUTPUT_NAME
        fi
        return 0
    else
        # No compression needed -- file not compressed
        echo $1
        return 0
    fi
}




# GET INPUT
while getopts ":r:g:w:m:n:fh" opt; do
    case $opt in
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
        h)
            print_usage
            exit
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


# CONSTANTS
DATETIME=$(date "+%Y%m%d_%X")
SYS_CORES=$(nproc)

# CHECK FOR PRESENCE OF POSITIONAL ARGUMENTS (further checks later)
if [[ $OPTIND > ${#@} ]]; then
    echo -e "Fatal:\t\tno sequence files passed as positional arguments (must be passed at the end of the line)." 1>&2
    print_usage
    exit 1
fi

# CHECK INPUT FROM FLAGS
[[ $GENOME_REF ]] || echo -e "Warning:\tno genome reference file (-r) specified; will not perform alignment." 1>&2
[[ $FEATURES_FILE ]] || echo -e "Warning:\tno genome feature file (-g) specified; will not perform annotation or feature counting." 1>&2

# verify the work directory or use default value
if [[ $WORK_DIR ]]; then
    WORK_DIR=$(readlink -f $WORK_DIR)
else
    echo -e "Info:\t\tno working directory (-d) specified; using '$PWD/miRNA-run_$DATETIME'" 1>&2
    WORK_DIR=$PWD"/miRNA-run_"$DATETIME
fi

# determine the number of threads to use
if [[ ! $NUM_CORES ]]; then
    echo -e "Info:\t\tnumber of cores not specified; setting to 1." 1>&2
    NUM_CORES=1
else
    if [[ $NUM_CORES =~ ^[0-9]+$ ]]; then
        if [[ $NUM_CORES -gt $SYS_CORES ]]; then
           echo -e "Warning:\tnumber of cores specified ($NUM_CORES) greater than number of cores available ($SYS_CORES). Setting to maximum $SYS_CORES." 1>&2
           NUM_CORES=$SYS_CORES
        fi
    else
        echo -e "Warning:\tnumber of cores must be a positive integer between 1 and $SYS_CORES. Setting number of cores to 1." 1>&2
       NUM_CORES=1
    fi
fi






# CREATE DIRECTORY TREE
LOG_DIR=$WORK_DIR"/logs/"
SEQDATA_DIR=$WORK_DIR"/seqdata/"
ALIGNED_DIR=$WORK_DIR"/aligned/"
[[ $FEATURES_FILE ]] && ANNOTATED_DIR=$WORK_DIR"/annotated/"
VIS_DIR=$WORK_DIR"/visualization/"
#LOG_FILE=$LOG_DIR/$INPUTFILE_BASE"_"$DATE".log"
for dir in $LOG_DIR $SEQDATA_DIR $ALIGNED_DIR $VIS_DIR $ANNOTATED_DIR; do
    if [[ $(mkdir -p $dir) -ne 0 ]]; then
        echo -e "Fatal:\t\tcannot create directory $dir; exiting." 1>&2
        exit 1
    fi
done




# GET/CHECK POSITIONAL ARGS (INPUT FILES)
if [[ $OPTIND == ${#@} ]]; then
    TMP_FILE="${@:$OPTIND:1}"
    if [[ ! $( exists_is_readable $TMP_FILE ) ]]; then
        echo -e "Fatal:\t\t\"$TMP_FILE\": "$FILE_NOT_EXISTS_OR_NOT_READABLE_ERROR_TEXT 1>&2
        exit 1
    fi

    TMP_FILE=$(decompress_file $TMP_FILE $SEQDATA_DIR)
    if [[ ! $TMP_FILE ]]; then
        TMP_FILE="${@:$OPTIND:1}"
        echo -e "Fatal:\t\tunhandled error when decompressing file \"$TMP_FILE\"; exiting." 1>&2
        exit 1
    fi

    if [[ ! $( extension_is_fastq $( basename $TMP_FILE) ) ]]; then
        echo -e "Fatal:\t\t\"$( basename $TMP_FILE)\": "$NOT_FASTQ_ERROR_TEXT 1>&2
        exit 1
    fi

    SEQ_FILE=$TMP_FILE
    echo "One file, \$SEQ_FILE is $SEQ_FILE" 1>&2

else
    # If more than one positional argument is passed, they must be merged
    for (( i=$OPTIND; i <= ${#@}; i++ )) {
        TMP_FILE="${@:$i:1}"
        TMP_FILE=$(decompress_file $TMP_FILE $SEQDATA_DIR)
        if [[ ! $( exists_is_readable $TMP_FILE ) ]]; then
            echo -e "Error:\t\tskipping file \"$TMP_FILE\": "$FILE_NOT_EXISTS_OR_NOT_READABLE_ERROR_TEXT 1>&2
            continue
        fi
        if [[ ! $( extension_is_fastq $TMP_FILE ) ]]; then
            echo -e "Error:\t\tskipping file \"$TMP_FILE\": "$NOT_FASTQ_ERROR_TEXT 1>&2
            continue
        fi

        if [[ ! $MERGE_FILE_NAME ]]; then
            MERGE_FILE_NAME=$(basename "${TMP_FILE%.*}_merged.${TMP_FILE##*.}")
            MERGE_FILE_ABSPATH=$SEQDATA_DIR"/"$MERGE_FILE_NAME
        fi

        echo -en "Info:\t\tConcatenating \"$TMP_FILE\" to merge file \"$MERGE_FILE_ABSPATH\"..." 1>&2
        cat $TMP_FILE >> $MERGE_FILE_ABSPATH
        echo " done." 1>&2
    }
    if [[ $MERGE_FILE_ABSPATH ]]; then
        SEQ_FILE=$MERGE_FILE_ABSPATH
        #echo "Multiple files, final merged file is at $SEQ_FILE" >&2
    else
        echo -e "Fatal:\t\tno valid input files found. Exiting." 1>&2
    fi
fi


# INPUT FILE NAME HANDLES
SCRIPT_SELF_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

INPUTFILE_ABSPATH=$(readlink -f $SEQ_FILE)
INPUTFILE=$(basename $SEQ_FILE)
INPUTFILE_DIR=$(dirname $INPUTFILE_ABSPATH)
INPUTFILE_BASE="${INPUTFILE%.*}"
INPUTFILE_EXTENSION="${INPUTFILE##*.}"

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


# LOG FILE
LOG_FILE=$LOG_DIR/$INPUTFILE_BASE"_"$DATETIME".log"
echo $0 $@ > $LOG_FILE


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
MIN_QUALITY_SCORE='10'
MIN_SEQ_LENGTH='18'
if [[ ! -f $OUTFILE_CUTADAPT ]] || [[ $FORCE_OVERWRITE ]]; then
    echo -e "Starting cutadapt trimming at $(date)." | tee -a $LOG_FILE 1>&2
    CL="cutadapt -f fastq -a $ADAPTER -q $MIN_QUALITY_SCORE --match-read-wildcards -O 5 -m $MIN_SEQ_LENGTH --too-short-output=/dev/null -o $OUTFILE_CUTADAPT $INFILE_CUTADAPT"
    echo "Cutadapt command: $CL" | tee -a $LOG_FILE 1>&2
    eval $CL 2>&1 | tee -a $LOG_FILE 1>&2
else
    echo -e "cutadapt adapter trimming already performed on $(basename $INFILE_CUTADAPT):\noutput file $(basename $OUTFILE_CUTADAPT)\nexists." | tee -a $LOG_FILE 1>&2
fi


echo -e "\n\nALIGNMENT TO GENOME REFERENCE\n=============================" | tee -a $LOG_FILE 1>&2
INFILE_ALN=$OUTFILE_CUTADAPT
OUTFILE_ALN=$ALIGNED_DIR"/"$INPUTFILE_BASE"_aln_"$REFERENCE_BASE".bam"
if [[ ! -f $OUTFILE_ALN ]]  || [[ $FORCE_OVERWRITE ]]; then
    echo -e "Started bowtie2 alignment to genome reference at $(date)." | tee -a $LOG_FILE 1>&2
    # Need some kinda tricky redirection I guess for this to send stderr to a file but stdout to samtools?
    CL="bowtie2 -N 1 -L 18 -p $NUM_CORES -x $REFERENCE_DIR/$REFERENCE_BASE $INFILE_ALN | samtools view -S -b - > $OUTFILE_ALN 2>>$LOG_FILE"
    echo "Executing alignment command: $CL" | tee -a $LOG_FILE
    eval $CL 2>>$LOG_FILE
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
    OUTFILE_ALN=$ALIGNED_DIR"/"$INPUTFILE_BASE"_aln_"$MIRBASE_BASE".bam"
    if [[ ! -f $OUTFILE_ALN ]] || [[ $FORCE_OVERWRITE ]]; then
        echo -e "Started bowtie2 alignment to miRBASE at $(date)." | tee -a $LOG_FILE 1>&2
        # Need some kinda tricky redirection I guess for this to send stderr to a file but stdout to samtools?
        CL="bowtie2 -N 1 -L 18 -p $NUM_CORES -x $MIRBASE_DIR/$MIRBASE_BASE $INFILE_ALN | samtools view -S -b - > $OUTFILE_ALN 2>>$LOG_FILE"
        echo "Executing alignment command: $CL" | tee -a $LOG_FILE
        eval $CL 2>>$LOG_FILE
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
