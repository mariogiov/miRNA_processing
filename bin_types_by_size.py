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
        feature_types = sorted(get_feature_types(ann_file).keys())
    else:
        feature_types = [feature_type]

    feature_counts = {}

    # Annotate
    for feature_type in feature_types:
        # TODO here it probably makes more sense to pull the gff file into memory to keep from reloading it each time
        # or to use the get_features_given_type to simulate a file only containing those features
        length_counts = annotate_count_aln(aln_file, ann_file, output_dir, feature_type)
        feature_counts[feature_type] = length_counts

    in_base, _ = os.path.splitext(os.path.basename(aln_file))
    # quick hack to dump output for later use
    tmp_output_file = os.path.join(output_dir, "feature_counts-{}.tsv".format(in_base))
    with open(tmp_output_file, 'w') as f:
        f.write( yaml.dump(feature_counts) )


def get_feature_types(annotation_file):
    """
    Given an annotation file, return a list of the feature types present
    and their frequencies.
    """
    gff = HTSeq.GFF_Reader(annotation_file)
    return collections.Counter( ( f.type for f in gff ) )


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

    # counts_reads_in_features prints to stdout
    with open(os.devnull, 'w') as f:
        with redirect_stream(sys.stdout, f):
            length_counts = None
            try:
                # Usage: count_reads_in_features( sam_filename, gff_filename, stranded, overlap_mode, feature_type, id_attribute, quiet, minaqual, samout )
                return counts.count_reads_in_features(aln_file, ann_file, False, "union", feature_type, 'ID', False, 10, "")
            except ValueError:
            # Skip features that give errors due to e.g. incorrect formatting, missing attributes
                return None

            #return length_counts, ann_sam_file, ann_cnt_file


def count_features_in_aln_manual(aln_file, ann_file, feature_type, output_dir=None, reference_file=None):
    """
    Given an alignment file, an annotation file, and an RNA type,
    returns the number of counts of that feature type present in the alignment.
    See http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html#counting-reads-by-genes
    """

    aln_basename, aln_ext   = os.path.splitext(os.path.basename(bam_file))
    ann_basename, _         = os.path.splitext(os.path.basename(annotation_file))

    # See http://www-huber.embl.de/users/anders/HTSeq/doc/alignments.html
    aln_reader_fn = {   '.bam': HTSeq.BAM_Reader,
                        '.sam': HTSeq.SAM_Reader,}[aln_ext]

    # Read SAM or BAM
    aln_reader = aln_reader_fn(aln_file)

    # Not keeping these files at the moment
    #if not output_dir:
    #    output_dir = os.path.dirname(bam_file)
    #    print("Output directory name not specified; using directory containing SAM file, \"{}\".".format(output_dir), file=sys.stderr)
    #aln_filename = os.path.join(output_dir, "{feature_type}_{sam_basename}.sam".format(**locals()))
    #cnt_filename = os.path.join(output_Dir, "{feature_type}_{sam_basename}.counts.csv".format(**locals()))

    feature_array = HTSeq.GenomicArrayOfSets( "auto", stranded=False )


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


def bam_to_sam(bam_file, ref_file, return_fh=True):
    """
    Given a BAM file, returns a SAM file with header
    or a filehandle.
    """
    # TODO
    pass


def sam_to_bam(sam_file, ref_file):
    """
    Given a SAM file and the correct reference genome,
    returns a sorted BAM file.
    """
    sam_file = os.path.realpath(sam_file)
    sam_basename, _ = os.path.splitext(os.path.basename(sam_file))
    # TODO

def annotate_alignment_no_htseq(sam_file, annotation_file, feature_type, idattr='ID', output_dir=None, reference_file=None):
    """
    Given an alignment file, an annotation file, and an RNA type,
    produces an annotated SAM file using HTSeq.
    Implemented as shell commands.
    """

    sam_file = os.path.realpath(sam_file)
    sam_basename, _ = os.path.splitext(os.path.basename(sam_file))

    annotation_file = os.path.realpath(annotation_file)
    annotation_basename, _ = os.path.splitext(os.path.basename(annotation_file))

    if not output_dir:
        output_dir = os.path.dirname(sam_file)
        print("Output directory name not specified; using directory containing SAM file, \"{}\".".format(output_dir))

    output_ann_aln_filename = os.path.join(output_dir, "{feature_type}_{sam_basename}.sam".format(**locals()))
    output_counts_filename = os.path.join(output_dir, "{feature_type}_{sam_basename}.counts.csv".format(**locals()))

    gff = HTSeq.GFF_Reader(annotation_file)

    #conversion_cmd = shlex.split("samtools view {}".format(sam_file))
    #htseqc_cmd = shlex.split("htseq-count -o {output_ann_aln_filename} -t {feature_type} -s no -q -i {idattr} - {annotation_file}".format(**locals()))
    #sort_cmd = shlex.split("sort -n -k 2 -r")

    #with open(output_counts_filename, 'w') as output_counts_fh:
    #    conversion_fh   = Popen(conversion_cmd, stdout=PIPE)
    #    #conversion_rc   = conversion_fh.wait()
    #    htseqc_fh        = Popen(htseqc_cmd, stdin=conversion_fh.stdout, stdout=PIPE)
    #    #htseq_rc        = htseq_fh.wait()
    #    sort_fh         = Popen(sort_cmd, stdin=htseqc_fh.stdout, stdout=output_counts_fh)
    #    #sort_rc         = sort_fh.wait()

    return (output_ann_aln_filename, output_counts_filename)

#def get_feature_types_no_htseq(annotation_file):
#    """
#    Given an annotation file, return a sorted list of the feature types present.
#    """
#
#    annotation_file = os.path.realpath(annotation_file)
#
#    # Anyone can use numpy, let's use Popen
#    cmd_cut = shlex.split('cut -f 3 {}'.format(annotation_file) )
#    cmd_grep = shlex.split('grep -v "^#"')
#    cmd_sort = shlex.split('sort')
#    cmd_uniq = shlex.split('uniq')
#
#    # Somewhere here I should be closing pipes, i.e. grep_out.stdout.close()
#    # See http://docs.python.org/2/library/subprocess.html#replacing-shell-pipeline
#    cut_out = Popen(cmd_cut, stdout=PIPE)
#    grep_out = Popen(cmd_grep, stdin=cut_out.stdout, stdout=PIPE)
#    sort_out = Popen(cmd_sort, stdin=grep_out.stdout, stdout=PIPE)
#    cmd_uniq = Popen(cmd_uniq, stdin=sort_out.stdout, stdout=PIPE)
#
#    return cmd_uniq.communicate()[0].split()


if __name__ == "__main__":
    parser = argparse.ArgumentParser("I got your docstring right here, pal")
    parser.add_argument("-a", "--annotation-file", dest="ann_file", required=True, help="The GTF/GFF annotation file to use (must match reference genome).")
    parser.add_argument("-r", "--reference-genome", dest="ref_file", help="The reference genome matching the alignment and the annotation file.")
    parser.add_argument("-s", "--alignment-file", dest="aln_file", help="The input alignment file.")
    parser.add_argument("-o", "--output-dir", dest="output_dir", help="The directory to use for writing output.")

    kwargs = vars(parser.parse_args())

    main(**kwargs)
