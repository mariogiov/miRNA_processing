
# GET INPUT
while getopts ":s:f" opt; do
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

# Check for positional args
echo "OPTIND $OPTIND"
echo "length ${#@}"

# List all the positional args
# See http://www.ibm.com/developerworks/opensource/library/l-bash-parameters/index.html
#     http://www.devhands.com/2010/01/handling-positional-and-non-positional-command-line-arguments-from-a-shell-script/
for ((i=$OPTIND; i<=${#@}; i++)) {
    #echo "$OPTIND: ${@:$OPTIND:1}"
    echo "$i: ${@:$i:1}"
}
