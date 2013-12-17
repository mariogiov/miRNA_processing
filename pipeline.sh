#!/bin/bash

# TODO the concatenation of input fastq files to merge will cause the python script to break if there is a frameshift e.g. from extra headers
# TODO consider checking if target files for decompression already exist? currently they are overwritten

# DEFINE FUNCTIONS

# print_usage()
function print_usage { echo -e  "\nUsage:\t$0\n" \
                                "\t\t[-r <genome_reference_file> (FASTA)>]\n" \
                                "\t\t[-g <genome_feature_file (GTF/GFF)>]\n" \
                                "\t\t[-m <mirbase_file (FASTA)>]\n" \
                                "\t\t[-o <output_directory>]\n" \
                                "\t\t[-n <cores>]\n" \
                                "\t\t[-f (overwrite existing files)]\n" \
                                "\t\t[-k (keep temp files)]\n" \
                                "\t\t<sequence_file> [<additional_sequence_files> <will_be_merged> <before_processing>]\n" >&2 ;
                     }

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
    if [[ ! -e $1 ]]; then
        return 1
    elif [[ ! -r $1 ]]; then
        return 1
    else
        echo "$1"
        return 0
    fi
}

# decompress_file()
function decompress_file () {
    if [[ $(file $1 | grep gzip ) ]]; then
        # If file is compressed
        echo -e "INFO:\t\tInput file \"$1\" is compressed; decompressing..." 1>&2
        DF_OUTPUT_NAME=$(basename "${1%.*}")
        if [[ $2 ]]; then
            # If output directory is given as an argument, decompress to this alternate output directory
            if [[ ! $(mkdir -p $2) -eq 0 ]]; then
                echo -e "WARNING:\tOutput directory \"$2\" could not be created. Decompressing to PWD \"$PWD\"." 1>&2
                OUTPUT_DIR=$PWD
            else
                OUTPUT_DIR=$(readlink -m $2)
            fi
            echo -en "INFO:\t\tGunzipping file \"$1\" to \"$OUTPUT_DIR/$DF_OUTPUT_NAME\"..." >&2
            gunzip -c $1 >> $OUTPUT_DIR/$DF_OUTPUT_NAME
            if [[ ! $? -eq 0 ]]; then
                echo -e "ERROR:\t\tCouldn't decompress file $1 to output directory \"$2\"." 1>&2
                return 1
            else
                echo " done." 1>&2
            fi
            echo $OUTPUT_DIR/$DF_OUTPUT_NAME
        else
            # Else unzip in place
            gzip -d $1
            if [[ ! $? -eq 0 ]]; then
                echo -e "ERROR:\t\tCouldn't decompress file $1." 1>&2
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

function fastqc_analysis() {
    FN_FQA_INPUT_FILE=$1
    FN_FQA_OUTPUT_DIR=$2
    if [[ ! -r $FN_FQA_INPUT_FILE ]]; then
        echo -e "\nWARNING:\tFastqc input file \"$FN_FQA_INPUT_FILE\" does not exists or cannot be read. Skipping fastqc analysis." 1>&2
        return 1
    elif [[ ! -e $FN_FQA_OUTPUT_DIR ]]; then
        if [[ ! $(mkdir -p $FN_FQA_OUTPUT_DIR) -eq 0 ]]; then
            echo -e "\nWARNING:\tFastqc output directory \"$FN_FQA_OUTPUT_DIR\" could not be created. Skipping fastqc analysis." 1>&2
            return 1
        else
            CL="fastqc -o $FN_FQA_OUTPUT_DIR -f fastq $FN_FQA_INPUT_FILE"
            echo -e "\nINFO:\t\tExecuting fastqc command line:\n\t\t$CL" 1>&2
            eval $CL 1>&2
            if [[ ! $? -eq 0 ]]; then
                echo -e "\nWARNING:\tFastqc analysis of input file \"$FN_FQA_INPUT_FILE\" failed." 1>&2
                return 1
            else
                echo $(readlink -e $FN_FQA_OUTPUT_DIR)
                return 0
            fi
        fi
    else
        echo -e "WARNING:\tFastqc output directory \"$FN_FQA_OUTPUT_DIR\" already exists. Skipping fastqc analysis." 1>&2
        return 1
    fi
}


# GET INPUT
while getopts ":r:g:o:m:n:fkh" opt; do
    case $opt in
        r)
            GENOME_REF=$OPTARG
            ;;
        g)
            FEATURES_FILE=$OPTARG
            ;;
        o)
            OUTPUT_DIR=$OPTARG
            ;;
        m)
            MIRBASE_FILE=$OPTARG
            ;;
        n)
            NUM_CORES=$OPTARG
            ;;
        t)
            TMP_DIR=$OPTARG
            ;;
        f)
            FORCE_OVERWRITE=1
            ;;
        k)
            KEEP_TMPFILES=1
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
DATETIME=$(date "+%Y%m%d_%H%M%S")
[[ $SLURM_CPUS_ON_NODE ]] && SYS_CORES=$SLURM_CPUS_ON_NODE || SYS_CORES=$(nproc --all)

# CHECK FOR PRESENCE OF POSITIONAL ARGUMENTS (further checks later)
if [[ $OPTIND > ${#@} ]]; then
    echo -e "FATAL:\t\tNo sequence files passed as positional arguments (must be passed at the end of the line)." 1>&2
    print_usage
    exit 1
fi

# CHECK INPUT FROM FLAGS
[[ $GENOME_REF ]]    || echo -e "WARNING:\tNo genome reference file (-r) specified; will not perform alignment to reference genome." 1>&2
[[ $FEATURES_FILE ]] || echo -e "WARNING:\tNo genome feature file (-g) specified; will not perform annotation or feature counting." 1>&2
[[ $MIRBASE_FILE ]]  || echo -e "WARNING:\tNo miRBase reference file (-m) specified; will not perform alignment to miRBase genome." 1>&2

# VERIFY OUTPUT DIRECTORY
if [[ $OUTPUT_DIR ]]; then
    OUTPUT_DIR=$(readlink -m $OUTPUT_DIR)
    if [[ ! $(mkdir -p $OUTPUT_DIR) -eq 0 ]]; then
        echo -e "FATAL:\t\tCannot create output directory $OUTPUT_DIR; exiting." 1>&2
        exit 1
    fi
else
    echo -e "INFO:\t\tNo working directory (-d) specified; using '$PWD/miRNA-run_$DATETIME'" 1>&2
    OUTPUT_DIR=$PWD"/miRNA-run_"$DATETIME
fi




# VERIFY TMP DIRECTORY IF PASSED OR USE DEFAULT SYSTEM TMP
if [[ $TMP_DIR ]]; then
    NEW_TMP_DIR=$(mkdir $(readlink -m $TMP_DIR)"/tmp/")
    if [[ !$NEW_TMP_DIR ]]; then
        echo -e "WARNING:\tUnable to create temporary directory in user-supplied directory $TMP_DIR; falling back to output directory \"$OUTPUT_DIR\"" 1>&2
    else
        TMP_DIR=$NEW_TMP_DIR
    fi
else
    # Try using environment variables to locate system tmp
    if [[ $TMPDIR ]]; then
        TMP_DIR=$(mktemp -d $TMPDIR"/tmp.XXX")
    elif [[ $SNIC_TMP ]]; then
        TMP_DIR=$(mktemp -d $SNIC_TMP"/tmp.XXX")
    else
        echo -e "INFO:\t\tCould not determine temporary directory location from environment variables: falling back to output directory \"$OUTPUT_DIR\"" 1>&2
    fi
fi
# this is kind of like a 'finally' clause in case all the other attempts fail
if [[ ! $TMP_DIR ]]; then
    TMP_DIR=$OUTPUT_DIR"/tmp/"
    mkdir -p $TMP_DIR
    if [[ ! $TMP_DIR ]]; then
        echo -e "FATAL:\t\tUnable to create temporary working directory \"$TMP_DIR\" -- please specify on the command line with the -t flag." 1>&2
        exit
    fi
fi

# DETERMINE THE NUMBER OF CORES TO USE
if [[ ! $NUM_CORES ]]; then
    echo -e "INFO:\t\tNumber of cores not specified; setting to 1." 1>&2
    NUM_CORES=1
else
    if [[ $NUM_CORES =~ ^[0-9]+$ ]]; then
        if [[ $NUM_CORES -gt $SYS_CORES ]]; then
           echo -e "WARNING:\tNumber of cores specified ($NUM_CORES) greater than number of cores available ($SYS_CORES). Setting to maximum $SYS_CORES." 1>&2
           NUM_CORES=$SYS_CORES
        fi
    else
        echo -e "WARNING:\tNumber of cores must be a positive integer between 1 and $SYS_CORES. Setting number of cores to 1." 1>&2
       NUM_CORES=1
    fi
fi



# CREATE DIRECTORY TREE
LOG_DIR=$OUTPUT_DIR"/logs/"
SEQDATA_DIR=$TMP_DIR"/seqdata/"
FASTQC_DIR=$OUTPUT_DIR"/fastqc/"
( [[ $GENOME_REF ]] || [[ $MIRBASE_FILE ]] ) && ALIGNED_DIR=$OUTPUT_DIR"/aligned/"
[[ $FEATURES_FILE ]] && ANNOTATED_DIR=$OUTPUT_DIR"/annotated/"
VIS_DIR=$OUTPUT_DIR"/visualization/"
for dir in $LOG_DIR $SEQDATA_DIR $ALIGNED_DIR $VIS_DIR $ANNOTATED_DIR; do
    if [[ $(mkdir -p $dir) -ne 0 ]]; then
        echo -e "FATAL:\t\tCannot create directory $dir; exiting." 1>&2
        exit 1
    fi
done


# CREATE TEMPORARY LOG FILE
# must create temporary log file until we get the definitive name of the input file
TMP_LOG_FILE=$(mktemp $LOG_DIR"run_"$DATETIME"-XXX.prelog")
START_TIME=$(date +%s)
echo -e "Beginning script execution at $(date)." > $TMP_LOG_FILE
echo -e "Script invoked as: $0 $@\n" >> $TMP_LOG_FILE


# GET/CHECK POSITIONAL ARGS (INPUT FILES)
if [[ $OPTIND == ${#@} ]]; then
    TMP_FILE="${@:$OPTIND:1}"
    if [[ ! $( exists_is_readable $TMP_FILE 2>>$TMP_LOG_FILE ) ]]; then
        echo -e "FATAL:\t\t\"$TMP_FILE\": "$FILE_NOT_EXISTS_OR_NOT_READABLE_ERROR_TEXT | tee -a $TMP_LOG_FILE 1>&2
        exit 1
    fi

    TMP_FILE=$(decompress_file $TMP_FILE $SEQDATA_DIR 2>>$TMP_LOG_FILE)
    if [[ ! $TMP_FILE ]]; then
        TMP_FILE="${@:$OPTIND:1}"
        echo -e "FATAL:\t\tUnhandled error when decompressing file \"$TMP_FILE\"; exiting." | tee -a $TMP_LOG_FILE 1>&2
        exit 1
    fi

    if [[ ! $( extension_is_fastq $( basename $TMP_FILE) ) ]]; then
        echo -e "FATAL:\t\t\"$( basename $TMP_FILE)\": "$NOT_FASTQ_ERROR_TEXT | tee -a $TMP_LOG_FILE 1>&2
        exit 1
    fi

    SEQ_FILE=$TMP_FILE
    #echo "Single input file \"$SEQ_FILE\"" 1>&2

else
    # If more than one positional argument is passed, they must be merged
    for (( i=$OPTIND; i <= ${#@}; i++ )) {
        TMP_FILE="${@:$i:1}"
        TMP_FILE=$(decompress_file $TMP_FILE $SEQDATA_DIR 2>>$TMP_LOG_FILE)
        if [[ ! $( exists_is_readable $TMP_FILE 2>>$TMP_LOG_FILE) ]]; then
            echo -e "ERROR:\t\tSkipping file \"$TMP_FILE\": "$FILE_NOT_EXISTS_OR_NOT_READABLE_ERROR_TEXT | tee -a $TMP_LOG_FILE 1>&2
            continue
        fi
        if [[ ! $( extension_is_fastq $TMP_FILE 2>>$TMP_LOG_FILE ) ]]; then
            echo -e "ERROR:\t\tSkipping file \"$TMP_FILE\": "$NOT_FASTQ_ERROR_TEXT | tee -a $TMP_LOG_FILE 1>&2
            continue
        fi

        if [[ ! $MERGE_FILE_NAME ]]; then
            MERGE_FILE_NAME=$(basename "${TMP_FILE%.*}_merged.${TMP_FILE##*.}")
            MERGE_FILE_ABSPATH=$SEQDATA_DIR"/"$MERGE_FILE_NAME
        fi

        echo -en "INFO:\t\tConcatenating \"$TMP_FILE\" to merge file \"$MERGE_FILE_ABSPATH\"..." | tee -a $TMP_LOG_FILE 1>&2
        cat $TMP_FILE >> $MERGE_FILE_ABSPATH 2>>$TMP_LOG_FILE
        if [[ $? -eq 0 ]]; then
        #if [[ $(cat $TMP_FILE >> $MERGE_FILE_ABSPATH 2>>$TMP_LOG_FILE) ]]; then
            echo " done." | tee -a $TMP_LOG_FILE 1>&2
        else
            echo -en "FATAL:\t\tCould not merge input files. Exiting." | tee -a $TMP_LOG_FILE 1>&2
            exit 1
        fi
    }
    if [[ $MERGE_FILE_ABSPATH ]]; then
        SEQ_FILE=$MERGE_FILE_ABSPATH
        #echo "Multiple files, final merged file is at $SEQ_FILE" >&2
    else
        echo -e "FATAL:\t\tNo valid input files found. Exiting." | tee -a $TMP_LOG_FILE 1>&2
    fi
fi


# INPUT FILE NAME HANDLES
SCRIPT_SELF_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


INPUTFILE_ABSPATH=$(readlink -m $SEQ_FILE)
INPUTFILE=$(basename $SEQ_FILE)
INPUTFILE_DIR=$(dirname $INPUTFILE_ABSPATH)
INPUTFILE_BASE="${INPUTFILE%.*}"
INPUTFILE_EXTENSION="${INPUTFILE##*.}"

if [[ $GENOME_REF ]]; then
    REFERENCE_ABSPATH=$(readlink -m $GENOME_REF)
    REFERENCE=$(basename $REFERENCE_ABSPATH)
    REFERENCE_DIR=$(dirname $REFERENCE_ABSPATH)
    REFERENCE_BASE="${REFERENCE%.*}"
    REFERENCE_EXTENSION="${REFERENCE##*.}"
fi

if [[ $FEATURES_FILE ]]; then
    FEATURES_FILE_ABSPATH=$(readlink -m $FEATURES_FILE)
    FEATURES_FILE=$(basename $FEATURES_FILE)
    FEATURES_FILE_DIR=$(dirname $FEATURES_FILE_ABSPATH)
    FEATURES_FILE_BASE="${FEATURES_FILE%.*}"
    FEATURES_FILE_EXTENSION="${FEATURES_FILE##*.}"
fi

if [[ $MIRBASE_FILE ]]; then
    MIRBASE_ABSPATH=$(readlink -m $MIRBASE_FILE)
    MIRBASE_FILE=$(basename $MIRBASE_FILE)
    MIRBASE_DIR=$(dirname $MIRBASE_ABSPATH)
    MIRBASE_BASE="${MIRBASE_FILE%.*}"
    MIRBASE_EXTENSION="${MIRBASE_FILE##*.}"
fi


# CREATE PERMANENT LOG FILE
LOG_FILE=$LOG_DIR/$INPUTFILE_BASE"_"$DATETIME".log"
mv $TMP_LOG_FILE $LOG_FILE



# MODULE LOADING
source $HOME/.virtualenvs/python-276/bin/activate

# Modules, activate the module command
case "$(basename $SHELL)" in
          -sh|sh|*/sh)  modules_shell=sh ;;
       -ksh|ksh|*/ksh)  modules_shell=ksh ;;
       -zsh|zsh|*/zsh)  modules_shell=zsh ;;
    -bash|bash|*/bash)  modules_shell=bash ;;
esac
module() { eval `/usr/local/Modules/$MODULE_VERSION/bin/modulecmd $modules_shell $*`; } 
export PATH=$HOME/bin:$HOME/.local/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lib/
export PYTHONPATH=$HOME/lib/python2.7/

module unload python
module load python/2.7.4
module load bowtie2/2.1.0
module load cutadapt
module load fastqc
module load htseq
module load samtools




echo -e "\n\nADAPTER TRIMMING\n================" | tee -a $LOG_FILE 1>&2
# Remove 3' adapter sequences, discarding reads shorter than $MIN_SEQ_LENGTH
INFILE_CUTADAPT=$INPUTFILE_ABSPATH
OUTFILE_CUTADAPT=$SEQDATA_DIR"/"$INPUTFILE_BASE"_trimmed."$INPUTFILE_EXTENSION
# TruSeq adapter sequence "TGGAATTCTCGGGTGCCAAGG"
ADAPTER="TGGAATTCTCGGGTGCCAAGG"
MIN_QUALITY_SCORE='10'
MIN_SEQ_LENGTH='18'
if [[ ! -f $OUTFILE_CUTADAPT ]] || [[ $FORCE_OVERWRITE ]]; then
    echo -e "INFO:\t\tStarting cutadapt trimming at $(date)." | tee -a $LOG_FILE 1>&2
    CL="cutadapt -f fastq -a $ADAPTER -q $MIN_QUALITY_SCORE --match-read-wildcards -O 5 -m $MIN_SEQ_LENGTH --too-short-output=/dev/null -o $OUTFILE_CUTADAPT $INFILE_CUTADAPT"
    echo -e "INFO:\t\tExecuting cutadapt command:\n\t\t$CL" | tee -a $LOG_FILE 1>&2
    eval $CL 2>&1 | tee -a $LOG_FILE 1>&2
    [[ ${PIPESTATUS[0]} -eq 0 ]] ||  echo -e "ERROR:\t\tcutadapt trimming failed." | tee -a $LOG_FILE 1>&2
else
    echo -e "INFO:\t\tcutadapt adapter trimming already performed on $(basename $INFILE_CUTADAPT):\noutput file $(basename $OUTFILE_CUTADAPT)\nexists." | tee -a $LOG_FILE 1>&2
fi





FASTQC_TRIMMED_DIR=$FASTQC_DIR
echo -e "\n\nTRIMMED READS FASTQC ANALYSIS\n=============================" | tee -a $LOG_FILE 1>&2
if [[ ! -f $FASTQC_DIR ]]; then
    echo -e "INFO:\t\tStarted fastqc analysis of trimmed reads at $(date)." | tee -a $LOG_FILE 1>&2
    CL="fastqc_analysis $OUTFILE_CUTADAPT $FASTQC_TRIMMED_DIR"
    echo -e "INFO:\t\tExecuting fastqc command:\n\t\t$CL" | tee -a $LOG_FILE 1>&2
    eval $CL 2>&1 | tee -a $LOG_FILE 1>&2
    [[ ${PIPESTATUS[0]} -eq 0 ]] || echo -e "ERROR:\t\tfastqc analysis of $OUTFILE_CUTADAPT failed." | tee -a $LOG_FILE 1>&2
else
    echo -e "INFO:\t\tfastqc analysis already performed on $(basename $OUTFILE_CUTADAPT):\noutput directory \"$FASTQ_TRIMMED_DIR\" exists." | tee -a $LOG_FILE 1>&2
fi





if [[ $GENOME_REF ]]; then
    echo -e "\n\nALIGNMENT TO GENOME REFERENCE\n=============================" | tee -a $LOG_FILE 1>&2
    INFILE_ALN=$OUTFILE_CUTADAPT
    OUTFILE_ALN=$ALIGNED_DIR"/"$INPUTFILE_BASE"_aln_"$REFERENCE_BASE".bam"
    if [[ ! -f $OUTFILE_ALN ]]  || [[ $FORCE_OVERWRITE ]]; then
        echo -e "INFO:\t\tStarted bowtie2 alignment to genome reference at $(date)." | tee -a $LOG_FILE 1>&2
        # Need some kinda tricky redirection I guess for this to send stderr to a file but stdout to samtools?
        # popen3 would work
        CL="bowtie2 -N 1 -L 18 -p $NUM_CORES -x $REFERENCE_DIR/$REFERENCE_BASE $INFILE_ALN | samtools view -S -b - > $OUTFILE_ALN"
        echo -e "INFO:\t\tExecuting alignment command:\n\t\t$CL" | tee -a $LOG_FILE 1>&2
        eval $CL 2>&1 | tee -a $LOG_FILE 1>&2
        [[ ${PIPESTATUS[0]} -eq 0 ]] || echo -e "ERROR:\t\tAlignment of $(basename $INFILE_ALN) to $(basename $REFERENCE_ABSPATH) failed." | tee -a $LOG_FILE 1>&2
    else
        echo -e "INFO:\t\tAlignment of $(basename $INFILE_ALN) to $(basename $REFERENCE_ABSPATH) already performed:\n$(basename $OUTFILE_ALN)\nexists." | tee -a $LOG_FILE 1>&2
    fi
fi



# Convert from BAM to SAM | annotate hits to genes, treating them as nonstranded, mode = union
if [[ $FEATURES_FILE ]]; then
    echo -e "\n\nANNOTATION\n==========" | tee -a $LOG_FILE 1>&2
    if [[ ! $OUTFILE_ALN ]]; then
        echo -e "ERROR:\t\tNo input alignment file available to annotate: please specify genome reference file via flag \"-r\"." | tee -a $LOG_FILE 1>&2
    else
        OUTFILE_ALN_BASE=$(basename "${OUTFILE_ALN%.*}")
        OUTFILE_COUNTS=$ANNOTATED_DIR"/"$OUTFILE_ALN_BASE"_counts_"$FEATURES_FILE_BASE".csv"
        ANNOTATED_FILE_SAM=$ANNOTATED_DIR"/"$OUTFILE_ALN_BASE"_annotated_"$FEATURES_FILE_BASE".sam"
        ANNOTATED_FILE_BAM="${ANNOTATED_FILE_SAM%.*}"".bam"
        # This binary logic gets a little wild but the idea is that if the annotation hasn't been done, or if it has but the bam file wasn't produced,
        # or if they all have but the force_overwrite (-f) flag was passed, process the files.
        # At least I -think- that's what it does.
        if ( ( [[ ! -f $OUTFILE_COUNTS ]] || [[ ! -f $ANNOTATED_FILE_SAM ]] ) || [[ ! -f $ANNOTATED_FILE_BAM ]] ) || [[ $FORCE_OVERWRITE ]]; then
            echo -e "INFO:\t\tStarted annotation of alignment at $(date)." | tee -a $LOG_FILE 1>&2
            CL="samtools view -h $OUTFILE_ALN | htseq-count -o $ANNOTATED_FILE_SAM -t exon -s no -q -i 'ID' - $FEATURES_FILE_ABSPATH | sort -n -k 2 -r > $OUTFILE_COUNTS"
            echo -e "INFO:\t\tExecuting annotation command:\n\t\t$CL" | tee -a $LOG_FILE 1>&2
            eval $CL 2>>$LOG_FILE
            [[ $? -eq 0 ]] || echo -e "ERROR:\t\tAnnotation failed." | tee -a $LOG_FILE 1>&2
        fi
        if [[ ! -f $ANNOTATED_FILE_BAM ]] || [[ $FORCE_OVERWRITE ]]; then
            if [[ ! $REFERENCE_ABSPATH ]]; then
                echo -e "WARNING:\tCannot convert SAM file to BAM format without reference file. Please specify the reference using the \"-r\" flag." | tee -a $LOG_FILE 1>&2
            else
                echo -e "INFO:\t\tConverting SAM file to BAM format." | tee -a $LOG_FILE 1>&2
                CL="samtools view -S -bT $REFERENCE_ABSPATH $ANNOTATED_FILE_SAM >> $ANNOTATED_FILE_BAM"
                echo -e "INFO:\t\tExecuting conversion command:\n\t\t$CL" | tee -a $LOG_FILE 1>&2
                eval $CL 2>>$LOG_FILE
                if [[ ! $? -eq 0 ]]; then
                    echo -e "ERROR:\t\tConversion of SAM file to BAM format failed." | tee -a $LOG_FILE 1>&2
                    # remove empty file on failed conversion but leave SAM file
                    rm $ANNOTATED_FILE_BAM
                else
                    # if we produce a new (annotated) bam file, delete the sam file and the original alignment file
                    rm $ANNOTATED_FILE_SAM $OUTFILE_ALN
                fi
            fi
        else
            echo -e "INFO:\t\tHtseq-count already performed: counts file $OUTFILE_COUNTS\nand annotated bam file $ANNOTATED_FILE_BAM\nboth exist." | tee -a $LOG_FILE 1>&2
        fi
    fi
fi




if [[ -f $OUTFILE_CUTADAPT ]]; then
    echo -e "\n\nVISUALIZATION\n=============" | tee -a $LOG_FILE 1>&2
    # Plot read length distribution with matplotlib
    VIS_FILE_NAME=$VIS_DIR/$(basename $OUTFILE_CUTADAPT .fastq)"_readlength-histogram.png"
    if [[ ! -f $VIS_FILE_NAME ]] || [[ $FORCE_OVERWRITE ]]; then
        echo -e "INFO:\t\tCreating read length distribution plot for $OUTFILE_CUTADAPT." | tee -a $LOG_FILE 1>&2
        CL="python $SCRIPT_SELF_DIR/plots.py -i $OUTFILE_CUTADAPT -d $VIS_DIR"
        echo "Executing: $CL" | tee -a $LOG_FILE 1>&2
        eval $CL | tee -a $LOG_FILE 1>&2
        [[ ${PIPESTATUS[0]} -eq 0 ]] || echo "ERROR:\t\tVisualization of read lengths failed." | tee -a $LOG_FILE 1>&2
    else
        echo -e "INFO:\t\tVisualization already performed: plot file\n$VIS_FILE_NAME\nexists."
    fi
else
    echo -e "WARNING:\tCould not find trimmed reads for plotting; skipping visualization." | tee -a $LOG_FILE 1>&2
fi





if [[ $MIRBASE_FILE ]]; then
    echo -e "\n\nALIGNMENT TO MIRBASE\n====================" | tee -a $LOG_FILE 1>&2
    INFILE_ALN=$OUTFILE_CUTADAPT
    OUTFILE_ALN=$ALIGNED_DIR"/"$INPUTFILE_BASE"_aln_miRBase.bam"
    if [[ ! -f $OUTFILE_ALN ]] || [[ $FORCE_OVERWRITE ]]; then
        echo -e "INFO:\t\tStarted bowtie2 alignment to miRBase at $(date)." | tee -a $LOG_FILE 1>&2
        # Need some kinda tricky redirection I guess for this to send stderr to a file but stdout to samtools?
        CL="bowtie2 -N 1 -L 18 -p $NUM_CORES -x $MIRBASE_DIR/$MIRBASE_BASE $INFILE_ALN | samtools view -S -b - > $OUTFILE_ALN 2>>$LOG_FILE"
        echo -e "INFO:\t\tExecuting alignment command:\n\t\t$CL" | tee -a $LOG_FILE
        eval $CL 2>>$LOG_FILE
        [[ ${PIPESTATUS[0]} -eq 0 ]] || echo -e "ERROR:\t\tAlignment of $(basename $INFILE_ALN) to $MIRBASE_FILE failed." | tee -a $LOG_FILE 1>&2
    else
        echo -e "INFO:\t\tAlignment of $(basename $INFILE_ALN) to $MIRBASE_FILE already performed:\n$(basename $OUTFILE_ALN)\nexists." | tee -a $LOG_FILE 1>&2
    fi
fi

# Remove tmp directory
if [[ $KEEP_TMPFILES ]]; then
    echo -e "\nINFO:\t\tKeeping temp files: moving to \"$OUTPUT_DIR/tmp\"..." | tee -a $LOG_FILE 1>&2
    if [[ ! "$OUTPUT_DIR/tmp" -eq $TMP_DIR ]]; then
        mkdir -p $OUTPUT_DIR/tmp/
        mv -v $TMP_DIR $OUTPUT_DIR/tmp/ 2>&1 | tee -a $LOG_FILE 1>&2
    fi
else
    echo -e "\nINFO:\t\tRemoving tmp files." | tee -a $LOG_FILE 1>&2
    rm -rf $TMP_DIR
fi

END_TIME=$(date +%s)
let DIFF_TIME=( $END_TIME-$START_TIME )/60

echo -e "INFO:\t\tRun finished at $(date); total processing time $DIFF_TIME minutes." | tee -a $LOG_FILE 1>&2
