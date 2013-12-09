"""
A module to produce that plot I've been talking about.
"""
import argparse
#import htseq
import os
import shlex

from subprocess import Popen, PIPE

def main(sam_file, annotation_file, reference_file, output_dir, feature_type=None):

    if not feature_type:
        feature_types = get_feature_types(annotation_file)
    else:
        feature_types = [feature_type]

    for feature_type in feature_types:
        print("Annotating alignment file \"{sam_file}\" with feature \"{feature_type}\" from annotation file \"{annotation_file}\"...".format(**locals()))
        annotated_file, counts_file = annotate_alignment(sam_file, annotation_file, feature_type, reference_file=reference_file)
        print ("Annotated file: {}".format(annotated_file))
        print ("Counts file: {}".format(counts_file))

def get_feature_types(annotation_file):
    """
    Given an annotation file, return a sorted list of the feature types present.
    """

    annotation_file = os.path.realpath(annotation_file)

    # Anyone can use numpy, let's use Popen
    cmd_cut = shlex.split('cut -f 3 {}'.format(annotation_file) )
    cmd_grep = shlex.split('grep -v "^#"')
    cmd_sort = shlex.split('sort')
    cmd_uniq = shlex.split('uniq')

    # Somewhere here I should be closing pipes, i.e. grep_out.stdout.close()
    # See http://docs.python.org/2/library/subprocess.html#replacing-shell-pipeline
    cut_out = Popen(cmd_cut, stdout=PIPE)
    grep_out = Popen(cmd_grep, stdin=cut_out.stdout, stdout=PIPE)
    sort_out = Popen(cmd_sort, stdin=grep_out.stdout, stdout=PIPE)
    cmd_uniq = Popen(cmd_uniq, stdin=sort_out.stdout, stdout=PIPE)

    return cmd_uniq.communicate()[0].split()

def annotate_alignment(sam_file, annotation_file, feature_type, idattr='ID', output_dir=None, reference_file=None):
    """
    Given an alignment file, an annotation file, and an RNA type,
    produces an annotated SAM file using HTSeq.
    """
    # I promise I'll implement this using the python modules themselves
    # eventually
    # probably

    sam_file = os.path.realpath(sam_file)
    sam_basename, _ = os.path.splitext(os.path.basename(sam_file))

    annotation_file = os.path.realpath(annotation_file)
    annotation_basename, _ = os.path.splitext(os.path.basename(annotation_file))

    if not output_dir:
        output_dir = os.path.dirname(sam_file)
        print("Output directory name not specified; using directory containing SAM file, \"{}\".".format(output_dir))

    output_ann_aln_filename = os.path.join(output_dir, "{feature_type}_{sam_basename}.sam".format(**locals()))
    output_counts_filename = os.path.join(output_dir, "{feature_type}_{sam_basename}.counts.csv".format(**locals()))

    conversion_cmd = shlex.split("samtools view {}".format(sam_file))
    htseqc_cmd = shlex.split("htseq-count -o {output_ann_aln_filename} -t {feature_type} -s no -q -i {idattr} - {annotation_file}".format(**locals()))
    sort_cmd = shlex.split("sort -n -k 2 -r")

    with open(output_counts_filename, 'w') as output_counts_fh:
        conversion_fh   = Popen(conversion_cmd, stdout=PIPE)
        #conversion_rc   = conversion_fh.wait()
        htseqc_fh        = Popen(htseqc_cmd, stdin=conversion_fh.stdout, stdout=PIPE)
        #htseq_rc        = htseq_fh.wait()
        sort_fh         = Popen(sort_cmd, stdin=htseqc_fh.stdout, stdout=output_counts_fh)
        #sort_rc         = sort_fh.wait()

    return (output_ann_aln_filename, output_counts_filename)


def sam_to_bam():
    """
    Given a SAM file and the correct reference genome,
    produces a BAM file.
    """
    pass

def sam_or_bam():
    """
    Given a file, determines if the file is in SAM or BAM format.
    """
    pass

def rna_counts_by_type():
    """
    Given an annotated SAM file, tallies the number of genes by
    their RNA type (e.g. tRNA, rRNA, exon)
    """
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser("I got your docstring right here, pal")
    parser.add_argument("-a", "--annotation-file", dest="annotation_file", required=True, help="The GTF/GFF annotation file to use (must match reference genome).")
    parser.add_argument("-r", "--reference-genome", dest="reference_file", help="The reference genome matching the alignment and the annotation file.")
    parser.add_argument("-s", "--sam-file", dest="sam_file", help="The input sam file.")
    parser.add_argument("-o", "--output-dir", dest="output_dir", help="The directory to use for writing output.")

    kwargs = vars(parser.parse_args())

    main(**kwargs)
