#!/usr/bin/python
"""
small RNA pipeline
"""

from __future__ import print_function

# I like to import my modules in alphabetical order but that's not an official python thing
import argparse
import datetime
import HTSeq
import os
import shlex
import sys
import subprocess

def main(genomeref_file, annotation_file, mirbase_file, output_dir, num_cores, force_overwrite, keep_tmp, tmp_dir, input_fastq_list):
    ## Docstrings are useful as they are what comes up when you do help(name_of_function) in a python interactive session.
    ## There is even a specific type of Docstring that we encourage folks to use around here that makes the autogeneration of PDF-
    ## or HTML-based docs extremely simple, and then gives us a free readthedocs.org website too!
    """
    Add docstring
    """

    ## here is that object we were talking about
    run_directory = RunDirectory(output_dir, tmp_dir, keep_tmp, force_overwrite)

    num_cores = get_validated_cores(num_cores)

    working_fastq = merge_input_fastq_files(input_fastq_list, run_directory)

    # Somehow activate the modules command, or use python functions

    # cutadapt

    # fastqc

    # bowtie2 alignment

    # annotate (htseq-count)

    # visualizations

    # mirbase alignment

    # remove or save tmp directory
    pass

def is_compressed(input_file):
    """
    Checks to see if the file is gzip compressed, returns T/F
    """
    # TODO add bzip2 support


    ## There are many ways to substitute variables into a string in Python,
    ## but this is the recommended way (and in python >3 the only way, I think).
    ## I think this is also covered well in Learning Python
    ## See http://infohost.nmt.edu/tcc/help/pubs/python/web/new-str-format.html
    cmd = shlex.split("gzip -d -t {}".format(input_file))
    try:
        ## subprocess is worth taking a look at in more depth,
        ## but the most commonly-used functions are
        ##   check_call
        ##   check_output
        ##   Popen
        ## A good guide: http://pymotw.com/2/subprocess/
        ## Official docs: https://docs.python.org/2/library/subprocess.html
        file_output = subprocess.check_call(cmd)

        # True and False are special names in python but only when written precisely that way
        # I think they are shortcuts to 1 and 0, in reality -- try True + 1 or False - 1
        return True

    ## This try-except paradigm below is very "Pythonic" --
    ## the idea is to ask for forgiveness (try-except) rather than permission (if-then)
    ## This is covered really well in that Learning Python book
    except subprocess.CalledProcessError:
    # Technically this is too broad a check (file not exists, file not readable)
    # but that is dealt with elsewhere
        return False

## We can talk about python parameters a bit more, but basically you can set default values in this way
def decompress_file(input_file, output_dir=os.getcwd(), return_pipe=False):
    """
    Decompresses an input file in gzip format.
    Returns the abspath to the decompressed output file
    or a PIPE if requested.
    """
    # TODO add bzip2 support

    ## You can evaluate the return value of a function as if it were that value
    if is_compressed(input_file):
        cmd = shlex.split(dcmp_str.format(fastq_file))
        # Remove .gz, .gzip, .bz2 extensions if present
        if return_pipe:
            ## See the guide referenced above for a treatment of buffer control via Popen
            return subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout
        else:
            basename = strip_gzip_ext(os.path.basename(fastq_file))

            ## os.path.join is the Pythonic way of doing this -- better than just
            ## appending a "/" as the function is platform-independent
            ## Actually the os module has lots of good functionality
            ## (os.path.isfile, os.mkdir, os.path.exists, etc.)
            output_file = os.path.join(run_directory.tmp_dir, basename)

            ## This with _fxn_ as _name_ usage is called a "context manager" and it is fantastic.
            ## The basic idea is that the context manager will execute some code on entry,
            ## and some code on exit no matter what occurs (which errors you encounter).
            ## It is the preferred way to deal with open files like this so you can
            ## be certain that the file will be closed when you exit the scope.
            with open(output_file, 'w') as f:
                subprocess.Popen(cmd, stdout=f)
            return output_file
    else:
        if return_pipe:

            ## shlex.split intelligently splits a string up into its composite pieces
            ## based on the assumption that they are a shell command, and returns
            ## a list (as required by subprocess.Popen).
            ## http://pymotw.com/2/shlex/
            cmd = shlex.split("cat {}".format(input_file))
            return subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout
        else:
            return os.path.realpath(input_file)

def strip_gzip_ext(file_name):
    """
    Returns the file without the .gz or .gzip suffix, if present
    """
    # TODO add bzip2 support

    base, ext = os.path.splitext(file_name)
    ## This is probably going to make you want to barf but you need to import a special module
    ## to do string regex matches in python, and it's called re.
    ## You can also use some shell-based matching with e.g. glob.glob and fnmatch.filter
    ## but those are shell-based wildcards and not really perl-style regexes
    ## https://developers.google.com/edu/python/regular-expressions
    ## http://www.pythonregex.com/
    if re.match('\.gz(ip)?$|.bz2', ".bz2", ext):
        file_name = base
    return file_name


def merge_input_fastq_files(input_fastq_list, run_directory):
    """
    Merge multiple fastq files into one fastq file.
    Returns the path to the final merged fastq file,
    or the original file in the case of just one input file.
    """
    if len(input_fastq_list) is 1:
        fastq_file = input_fastq_list[0]
        # decompress if compressed
        return decompress_file(fastq_file, output_dir=run_directory.tmp_dir)
    else:
        ## more string formatting
        merged_filename = "MERGED_{}".format(strip_gzip_ext(os.path.basename(input_fastq_list[0])))
        merged_filepath = os.path.join(run_directory.tmp_dir, merged_filename)
        with open(merged_filepath) as output_file:
            for fastq_file in input_fastq_list:
                stream = decompress_file(fastq_file, return_pipe=True)
                output_file.write(stream.stdout.read())
        return merged_filepath

def get_validated_cores(num_cores):
    """Get the number of cores available on the current node."""
    ## as Python is whitespace-sensitive, you need these \ at the end to continue onto the next line (as with bash)
    sys_cores =     os.getenv('SLURM_CPUS_ON_NODE') \
                or  subprocess.check_output(['nproc', '--all']) \
                or  1
    sys_cores = int(sys_cores)

    if not num_cores:
        num_cores = sys_cores
    if num_cores > sys_cores:
        ## This **locals business just expands the dictionary (hash) returned by the locals() function
        ## and the local variables num_Cores and sys_cores wind up matching the string formatting names
        print(  "Requested number of cores ({num_cores}) greater than number of system cores ({sys_cores}); " \
                "using {sys_cores} instead.".format(**locals()), file=sys.stderr)
        num_cores = sys_cores
    if num_cores < 1:
        print(  "Requested number of cores ({num_cores}) must be a postive integer; " \
                "using 1 instead.".format(**locals()), file=sys.stderr)
        num_cores = 1
    return num_cores

class RunDirectory(object):
    """
    Keeps track of the various directories used for work and output.
    """
    
    ## Here I was thinking about making this a context manager
    ## but I don't really think that's such a great idea
    ## I think I had just learned how to write context managers

    # TODO I wonder if this is worth making into a context manager
    #       e.g. create tmp on entry, remove on exit

    def __init__(self, output_dir, tmp_dir, keep_tmp, force_overwrite):
        self.output_dir         = self.create_output_dir(output_dir)
        self.tmp_dir            = self.create_tmp_dir(tmp_dir)
        self.force_overwrite    = force_overwrite
        self.keep_tmp           = keep_tmp
        self.res_dirs           = []
        self.tmp_dirs           = []

    def create_output_dir(self, output_dir):
        """
        Create the output directory passed in by the user
        or one following the format ./"smRNA_run_(datetime)/"
        Returns the absolute path to the output directory.
        """
        if not output_dir:
            output_dir = os.path.join(  os.getcwd(),
                                        "smRNA_run_{}".format(datetime.date.strftime(
                                        datetime.datetime.now(), format="%Y%m%d_%X")))
        else:
            output_dir = os.path.realpath(output_dir)

        try:
            print('Creating output directory "{}"'.format(output_dir), file=sys.stderr)
            os.makedirs(output_dir)
        except OSError as e:
            if e.errno is 17:
            # Output directory already exists
                pass
        self.output_dir = output_dir
        return self.output_dir


    def create_tmp_dir(self, tmp_dir):
        """
        Create a temporary working directory; use the value passed by the user
        else system tmp if determinable else output directory.
        Returns the absolute path to the tmp directory.
        """
        # if keep_tmp, still write to tmp_dir so as not to hammer network drives, etc. (?)
        if not tmp_dir:
            # try using environment vars to locate system tmp
            tmp_dir = os.getenv('TMPDIR') or os.getenv('SNIC_TMP') or self.output_dir
        tmp_dir = os.path.join(os.path.realpath(tmp_dir), "tmp")

        try:
            print('Creating tmp directory "{}"'.format(tmp_dir), file=sys.stderr)
            os.makedirs(tmp_dir)
        ## This "except Error as e" thing labels the error object as e so you can use it later
        except OSError as e:
            if e.errno is 17:
            # Tmp dir already exists, which should be fine
                pass

        self.tmp_dir = tmp_dir
        return self.tmp_dir

    def create_dir(self, dir_name, in_tmp=False):
        """
        Create a new directory within the working or tmp directory.
        """
        if in_tmp:
            full_dir_path = os.path.join(self.tmp_dir, dir_name)
        else:
            full_dir_path = os.path.join(self.output_dir, dir_name)

        print('Creating directory "{}"'.format(dir_name), file=sys.stderr)
        ## makedirs is like mkdir -p (create intervening directories)
        os.makedirs(full_dir_path)

        if in_tmp:
            self.tmp_dirs.append(full_dir_path)
        else:
            self.dirs.append(full_dir_paht)

        return full_dir_path


    def remove_tmp_dir(self):
        """
        Remove the tmp directory.
        """
        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Execute the small RNA pipeline.")
    # TODO should allow multiple reference genomes but then need to determine how annotation files and reference files are linked
    parser.add_argument("-r", "--genome-reference-file", dest="genomeref_file",
                        help="The genome reference file against which to align.")
    parser.add_argument("-g", "--genome-feature-file", dest="annotation_file",
                        help="GTF/GFF genome feature file to use for annotation (must match reference file).")
    parser.add_argument("-m", "--mirbase-file", dest="mirbase_file",
                        help="The miRBase reference file.")
    parser.add_argument("-o", "--output-dir", dest="output_dir",
                        help="The output directory.")
    parser.add_argument("-n", "--num-cores", dest="num_cores", type=int,
                        help="The number of cores to use for alignment.")
    parser.add_argument("-t", "--tmp-dir", dest="tmp_dir",
                        help="Optionally specify a temporary directory.")
    parser.add_argument("-f", "--force-overwrite", dest="force_overwrite", action="store_true", default=False,
                        help="Force overwrite of existing files.")
    parser.add_argument("-k", "--keep-tmp", dest="keep_tmp", action="store_true", default=False,
                        help="Keep temporary files after processing.")
    parser.add_argument("input_fastq_list", nargs="+", required=True)

    kwargs = vars(parser.parse_args())
    main(**kwargs)
