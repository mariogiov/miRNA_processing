#!/usr/bin/python

import argparse
import os

from itertools import islice
from matplotlib import pyplot as plt

def hist_readlen(seq_file, output_name=None, output_directory=None):
    """Plots a histogram of read lengths from a fastq file.
    """
    hist_lengths = []
    first_headerline = 0
    with open(seq_file, 'r') as f:
        for num, line in enumerate(f):
            if line.startswith('@'):
                first_headerline = num
    with open(seq_file, 'r') as f:
        lengths = [ len(seq) for seq in islice(f, first_headerline, None, 4) ]
    fig = plt.figure()
    axes = fig.add_subplot(111)
    plt.hist(lengths)
    if not output_directory:
        output_directory = os.getcwd()
    if not output_name:
        output_name, _ = os.path.splitext(seq_file)
    output_path = os.path.join(output_directory, output_name)
    fig.save(output_path)

if __name__=="__main__":
    parser = argparse.ArgumentParser("Create a histogram of read lengths from a given fastq file.")
    parser.add_argument("-i", "--input-file", dest="input_file", required=True, help="The fastq file from which to read.")
    parser.add_argument("-d", "--output-dir", dest="output_dir", default=os.getcwd(), help="The directory in which to save the output file.")
    parser.add_argument("-n", "--output-name", dest="output_name", default=None,  help="The name to use for the output file. Default same name as input file.")

    args = vars(parser.parse_args())
    hist_readlen(input_file, output_name, output_dir)
