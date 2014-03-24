#!/usr/bin/python
"""
Plot a histogram of readlengths from a fastq file.
"""

import argparse
import gzip
import os
import re

from collections import defaultdict
from itertools import islice
from matplotlib import pyplot as plt

def hist_readlen(seq_file, metrics=False, output_name=None, output_dir=None):
    """Plots a histogram of read lengths from a fastq file.
    """
    seq_file, output_dir = [ os.path.realpath(path) for path in seq_file, output_dir ]
    seq_file_name, seq_file_ext = os.path.splitext( os.path.basename(seq_file) )
    seq_file_basename = seq_file_name.replace("_trimmed", "")
    fastq_data_handle = FastQParser(seq_file)
    lengths_dict = defaultdict(int)
    for seq in fastq_data_handle:
        seq_length = seq[1].strip()
        lengths_dict[len(seq_length)] += 1
    # Figure business
    fig = plt.figure()
    axes = fig.add_subplot(111)
    plt.bar(lengths_dict.keys(), lengths_dict.values(), align='center')
    #plt.hist(lengths, bins=len(lengths))
    axes.set_xlabel("Read Lengths (After Adapter Trimming)")
    axes.set_ylabel("# Sequences with Length X")
    axes.set_yticklabels( [ "{}k".format(int(x) / 1000) for x in axes.get_yticks() ] )
    axes.set_title("Read Length Histogram\n{}".format(seq_file_basename))
    output_format="png"
    if not output_dir:
        output_dir = os.getcwd()
    if not output_name:
        output_name = seq_file_name + "_readlength-histogram" + "." + output_format
    output_path = os.path.join(output_dir, output_name)
    fig.savefig(output_path)

         
class FastQParser:
    """Parser for fastq files, possibly compressed with gzip.
    Iterates over one record at a time. A record consists
    of a list with 4 elements corresponding to 1) Header,
    2) Nucleotide sequence, 3) Optional header, 4) Qualities
    """
    def __init__(self,file,filter=None):
        self.fname = file
        self.filter = filter
        fh = open(file,"rb")
        if file.endswith(".gz"):
            self._fh = gzip.GzipFile(fileobj=fh)
        else:
            self._fh = fh
        self._records_read = 0
        self._next = self.setup_next()
    def __iter__(self):
        return self
    def next(self):
        return self._next(self)
    def setup_next(self):
        """Return the function to return the next record
        """
        if self.filter is None or len(self.filter.keys()) == 0:
            def _next(self):
                self._records_read += 1
                return [self._fh.next().strip() for n in range(4)]
        else:
            def _next(self):
                while True:
                    record = [self._fh.next().strip() for n in range(4)]
                    header = parse_header(record[0])
                    skip = False
                    for k, v in self.filter.items():
                        if k in header and header[k] not in v:
                            skip = True
                            break
                    if not skip:
                        self._records_read += 1
                        return record
        return _next
    def name(self):
        return self.fname
    def rread(self):
        return self._records_read
    def seek(self,offset,whence=None):
        self._fh.seek(offset,whence)
    def close(self):
        self._fh.close()

if __name__=="__main__":
    parser = argparse.ArgumentParser("Create a histogram of read lengths from a given fastq file.")
    parser.add_argument("-i", "--input-file", dest="seq_file", required=True, help="The fastq file from which to read.")
    parser.add_argument("-m", "--metrics", action="store_true", help="Produce a metrics file.")
    parser.add_argument("-d", "--output-dir", default=os.getcwd(), help="The directory in which to save the output file.")
    parser.add_argument("-n", "--output-name", default=None,  help="The name to use for the output file (not including extension). Default same name as input file.")

    args = vars(parser.parse_args())
    hist_readlen(**args)
