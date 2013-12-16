#!/usr/bin/python
"""
small RNA pipeline
"""

from __future__ import print_function

import argparse
import os

def main(genomeref_file, annotation_file, mirbase_file, output_dir, num_cores, force_overwrite, keep_tmp, tmp_dir, input_fastq_list):
    # Create output directory tree or fail
    # Initialize logger
    # Set up tmp directory
    # Verify / get cores to use

    # Merge input fastq files if > 1
    # Somehow activate the modules command, or use python functions

    # cutadapt

    # fastqc

    # bowtie2 alignment

    # annotate (htseq-count)

    # visualizations

    # mirbase alignment

    # remove or save tmp directory
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
    parser.add_argument("-n", "--num-cores", dest="num_cores",
                        help="The number of cores to use for alignment.")
    parser.add_argument("-f", "--force-overwrite", dest="force_overwrite", action="store_true", default=False,
                        help="Force overwrite of existing files.")
    parser.add_argument("-k", "--keep-tmp", dest="keep_tmp", action="store_true", default=False,
                        help="Keep temporary files after processing.")
    parser.add_argument("-t", "--tmp-dir", dest="tmp_dir",
                        help="Optionally specify a temporary directory.")
    parser.add_argument("input_fastq_list", nargs="+")

    kwargs = vars(parser.parse_args())
    main(**kwargs)
