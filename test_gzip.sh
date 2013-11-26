# decompress_file()
function decompress_file () {
    echo "Fn: \$1 is $1" >&2
    echo "fn: file \$1 is $(file $1)" >&2
    if [[ $(file $1 | grep gzip ) ]]; then
        echo -e "Info:\tinput file "$1" is compressed; decompressing..." | tee -a $LOG_FILE 1>&2
        gzip -d $1
        R_V="${1%.*}"
        echo "Fn: \$R_V is $R_V" >&2
        echo $R_V
    fi
}

RETURN_VAL=$( decompress_file $1 )
if [[ $RETURN_VAL ]]; then
    [[ -e $RETURN_VAL ]] && echo "RETURN_VAL exists: $(readlink -e $RETURN_VAL)" || echo "RETURN_VAL not a file."
fi
