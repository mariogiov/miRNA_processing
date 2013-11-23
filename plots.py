#!/usr/bin/python

import argparse
import os

from itertools import islice
from matplotlib import pyplot as plt

def hist_readlen(seq_file, output_name=None, output_dir=None):
    """Plots a histogram of read lengths from a fastq file.
    """

    seq_file_name, seq_file_ext = os.path.splitext( os.path.basename(seq_file) )

    # Find the first header line
    first_headerline = 0
    with open(seq_file, 'r') as f:
        for num, line in enumerate(f):
            if line.startswith('@'):
                first_headerline = num
                break
    # Offset 1 from the header line, get every 4th
    with open(seq_file, 'r') as f:
        lengths = [ len(seq) for seq in islice(f, first_headerline+1, None, 4) ]

    # Figure business
    fig = plt.figure()
    axes = fig.add_subplot(111)
    plt.hist(lengths, bins=len(lengths))
    axes.set_xlabel("Read Lengths (After Adapter Trimming)")
    axes.set_ylabel("# Sequences with Length X")
    axes.set_title("Read Length Histogram\n{}".format(seq_file_name))

    if not output_dir:
        output_dir = os.getcwd()
    if not output_name:
        output_name = seq_file_name + "_readlength-histogram"
    output_path = os.path.join(output_dir, output_name)
    fig.savefig(output_path)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Create a histogram of read lengths from a given fastq file.")
    parser.add_argument("-i", "--input-file", dest="seq_file", required=True, help="The fastq file from which to read.")
    parser.add_argument("-d", "--output-dir", dest="output_dir", default=os.getcwd(), help="The directory in which to save the output file.")
    parser.add_argument("-n", "--output-name", dest="output_name", default=None,  help="The name to use for the output file. Default same name as input file.")

    args = vars(parser.parse_args())
    #hist_readlen(input_file], args[output_name], args[output_dir])
    hist_readlen(**args)