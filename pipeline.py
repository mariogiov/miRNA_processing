#!/usr/bin/python
"""
small RNA pipeline
"""

from __future__ import print_function

import argparse
import datetime
import HTSeq
import os
import shlex
import sys
import subprocess

def main(genomeref_file, annotation_file, mirbase_file, output_dir, num_cores, force_overwrite, keep_tmp, tmp_dir, input_fastq_list):
    """
    Add docstring
    """

    rundirectory = RunDirectory(output_dir, tmp_dir, keep_tmp, force_overwrite)
    # TODO Initialize logger -- see how Guillermo tends to do this -- need to write to local file, stderr, and logstash
    # logger = 


    # Verify / get cores to use
    num_cores = get_validated_cores(num_cores)

    # Merge input fastq files if > 1
    working_fastq = merge_input_fastq_files(input_fastq_list)

    # Somehow activate the modules command, or use python functions

    # cutadapt

    # fastqc

    # bowtie2 alignment

    # annotate (htseq-count)

    # visualizations

    # mirbase alignment

    # remove or save tmp directory
    pass

def get_validated_cores(num_cores):
    sys_cores =     os.getenv('SLURM_CPUS_ON_NODE') \
                or  subprocess.check_output(['nproc', '--all']) \
                or  1
    sys_cores = int(sys_cores)

    if not num_cores:
        num_cores = sys_cores
    if num_cores > sys_cores:
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
    # TODO I wonder if this is worth making into a context manager
    #       e.g. create tmp on entry, remove on exit

    def __init__(self, output_dir, tmp_dir, keep_tmp, force_overwrite):
        self.output_dir         = self.create_output_dir(output_dir)
        self.tmp_dir            = self.create_tmp_dir(tmp_dir)
        self.force_overwrite    = force_overwrite
        self.keep_tmp           = keep_tmp

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
        # TODO how to deal with keep_tmp -- write to output_dir from the start or move after?
        if not tmp_dir:
            # try using environment vars to locate system tmp
            tmp_dir = os.getenv('TMPDIR') or os.getenv('SNIC_TMP') or self.output_dir
        tmp_dir = os.path.join(os.path.realpath(tmp_dir), "tmp")

        try:
            print('Creating tmp directory "{}"'.format(tmp_dir), file=sys.stderr)
            os.makedirs(tmp_dir)
        except OSError as e:
            if e.errno is 17:
            # Tmp dir already exists, which should be fine
                pass

        self.tmp_dir = tmp_dir
        return self.tmp_dir

    def remove_tmp_dir(self):
        """
        Remove the tmp directory.
        """

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
