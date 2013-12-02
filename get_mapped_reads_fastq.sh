#!/bin/bash

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
module load samtools
module load picard


input_file_fullpath=$(readlink -f $1)
input_file_basename=$(basename $input_file_fullpath)
output_dir_fullpath=$(dirname $input_file_fullpath)
output_vis_dir=$(dirname $(readlink -f $output_dir_fullpath))/visualization/
output_base="${input_file_basename%.*}""_mappedreads"
output_samfile=$output_dir_fullpath/$output_base".sam"
output_fastqfile=$output_dir_fullpath/$output_base".fastq"


echo "Creating output sam file \"$output_samfile\""
# -h outputs sam header, required by Picard
[[ ! -e $output_samfile ]] && samtools view -h -F 4 $input_file_fullpath >> $output_samfile || echo "Sam file \"$output_samfile\" already exists; skipping sam file generation."
echo "Done with sam file generation."

echo -n "Creating output fastqfile \"$output_fastqfile\""
[[ ! -e $output_fastqfile ]] && java -Xmx250M -jar $PICARD_HOME/SamToFastq.jar I=$output_samfile F=$output_fastqfile || echo "Fastq file \"$output_fastqfile\" already exists; skipping fastq file generation."
echo "Done with fastq file generation."

echo "Creating output plot under $output_vis_dir"
python plots.py -i $output_fastqfile -d $output_vis_dir
echo "Done with output plot generation."
