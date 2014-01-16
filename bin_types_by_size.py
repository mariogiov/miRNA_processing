#!/usr/bin/python
"""
A module to produce that plot I've been talking about.
"""
from __future__ import print_function

import argparse
import collections
import counts
import HTSeq
import numpy as np
import os
import shlex
import sys
import yaml

#from cStringIO import StringIO
from subprocess import Popen, PIPE

def main(aln_file, ann_file, ref_file, output_dir=None, id_attr='ID', feature_type=None):
    if not output_dir:
        output_dir = os.path.join(os.getcwd(), "annotated")
        print("No output directory specified; using {output_dir}".format(**locals()), file=sys.stderr)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Get feature types from GFF file
    if not feature_type:
        gff = HTSeq.GFF_Reader(ann_file)
        feature_types = sorted(collections.Counter((f.type for f in gff)))
    else:
        feature_types = [feature_type]
    feature_counts = {}
    # Annotate
    for feature_type in feature_types:
        # TODO here it probably makes more sense to pull the gff file into memory to keep from reloading it each time
        length_counts = annotate_count_aln(aln_file, ann_file, output_dir, feature_type)
        feature_counts[feature_type] = length_counts
    in_base, _ = os.path.splitext(os.path.basename(aln_file))
    # quick hack to dump output for later use
    tmp_output_file = os.path.join(output_dir, "feature_counts-{}.tsv".format(in_base))
    with open(tmp_output_file, 'w') as f:
        f.write( yaml.dump(feature_counts) )


# TODO this whole function could be a generator comprehension
def get_features_given_type(gff_file, feature_type):
    """
    Given a feature type (e.g. "exon") return all the features in the GFF file
    corresponding to that type as a dict of interval->feature.
    """
    gff = HTSeq.GFF_Reader(gff_file)
    features_dict = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    for feature in gff:
        if feature.type == "feature_type":
            features_dict[ feature.iv ] += feature.name
    return features_dict

def annotate_count_aln(aln_file, ann_file, output_dir, feature_type):
    """
    Annotated some alignments, and count that shit
    """
    aln_file_base, _    = os.path.splitext(os.path.basename(aln_file))
    ann_file_base, _    = os.path.splitext(os.path.basename(ann_file))
    ann_sam_file        = os.path.join(output_dir, "{aln_file_base}_annotated_{ann_file_base}-{feature_type}.sam".format(**locals()))
    ann_cnt_file        = os.path.join(output_dir, "{aln_file_base}_annotated_{ann_file_base}-{feature_type}.counts.tsv".format(**locals()))
    print("Getting counts for feature type \"{feature_type}\" in alignment file \"{aln_file}\" "
          "from annotation file \"{ann_file}\"".format( feature_type=feature_type, aln_file=os.path.basename(aln_file),
                                                        ann_file=os.path.basename(ann_file)), file=sys.stderr)
    with open(os.devnull, 'w') as null:
        # counts.counts_reads_in_features prints things to stdout --> redirect to /dev/null
        with redirect_stream(sys.stdout, null):
            try:
                # Usage: count_reads_in_features( sam_filename, gff_filename, stranded, overlap_mode, feature_type, id_attribute, quiet, minaqual, samout )
                return counts.count_reads_in_features(aln_file, ann_file, False, "union", feature_type, 'ID', False, 10, "")
            except ValueError:
            # Skip features that give errors due to e.g. incorrect formatting, missing attributes
                return None


class RedirectStream(object):
    """
    Context manager to swap one stream for another temporarily. Useful for e.g. temporarily redirecting STDOUT to a file.
    """
    def __init__(self, original_stream, tmp_stream):
        self.original_stream    = original_stream
        self.tmp_stream         = tmp_stream

    def __enter__(self):
        original_stream = self.tmp_stream

    def __exit__(self):
        original_stream = self.original_stream


if __name__ == "__main__":
    parser = argparse.ArgumentParser("I got your docstring right here, pal")
    parser.add_argument("-a", "--annotation-file", dest="ann_file", required=True, help="The GTF/GFF annotation file to use (must match reference genome).")
    parser.add_argument("-r", "--reference-genome", dest="ref_file", help="The reference genome matching the alignment and the annotation file.")
    parser.add_argument("-s", "--alignment-file", dest="aln_file", help="The input alignment file.")
    parser.add_argument("-o", "--output-dir", dest="output_dir", help="The directory to use for writing output.")

    kwargs = vars(parser.parse_args())
    main(**kwargs)
