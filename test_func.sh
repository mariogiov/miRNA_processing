function test_is_fastq () {
    echo "Fn: \$1 is $1" >&2
    TIF_EXT="${1##*.}"
    echo "Fn: \$TIF_EXT $TIF_EXT" >&2
    if [[ $TIF_EXT == "fq" ]] || [[ $TIF_EXT == "fastq" ]]; then
        echo "Fn: True" >&2
        echo $TIF_EXT
        return 0
    else
        echo "Fn: False" >&2
        return 1
    fi
}

echo $1
echo "$(test_is_fastq $1)"
if [[ $( test_is_fastq $1 ) ]]; then
    echo "file is fastq"
else
    echo "file is not fastq"
fi
