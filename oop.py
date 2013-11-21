#! /usr/bin/env python2.7

import HTSeq
import sys
import pysam
import gzip
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib.pyplot as plt
import pylab
import csv


class Arguments:
	def __init__(self):
		self.command = sys.argv[1]
		self.in_handle = sys.argv[2]
		#self.in_handle = open(sys.argv[2], "rU")
		self.out_handle = open(sys.argv[3], "w")
		
	def inHandle(self):
		return self.inHandle

class SAMFile:
	def __init__(self):
		self.args = Arguments()
	
	def change_sam(self):
		pass

class FastqFile:
	def __init__(self):
		self.lengths = []
		self.args = Arguments()

	def move_barcode(self):
		pass
		
	def MoveBarcode(self):
		pass
	
	def GetLengths(self):
		self.__in_handle = open(self.args.in_handle, "rU")
		for rec in SeqIO.parse(self.__in_handle, "fastq"):
			self.lengths.append(len(rec.seq))
           
	def PlotHist(self):
		if len(self.lengths) == 0: 
			self.GetLengths()
		self.Plots(self.lengths, title, xlabel, ylabel, out_handle)
	
	def CountBarcodes(self):
		__list = {}
		__tag = ''
		for record in HTSeq.FastqReader(self.args.in_handle):
			__tag =  record.name.split(' ')[1][6:]
			if __tag in __list:
				__list[__tag] = __list[__tag] + 1
			else:
				__list[__tag] = 1
		writer = csv.writer(self.args.out_handle, delimiter = '\t')
		for barcode in __list:
			writer.writerow([barcode, __list[barcode]])
		#writer.writerow(['klart', 'klart'])


class BarcodedFq(FastqFile):
	def __init__(self):
		self.lengths = []
		self.args = Arguments()
		
	#Comment
	def MoveBarcode(self):
		#print self.args.in_handle
		for record in HTSeq.FastqReader(self.args.in_handle):
			tag =  record.seq[:8]
			record = record[8:]
			record.name = record.name.split(' ')[0] + tag
			record.write_to_fastq_file(self.args.out_handle)


class BarcodedSAM(SAMFile):
	def __init__(self):
		self.args = Arguments()


	def MoveBarcode(self):
		__header = open(self.args.in_handle, 'rU')
		__first = '@'
		__line = __header.readline()
		while __first == '@':
			self.args.out_handle.write(__line)
			__line = __header.readline()
			__first = __line[0]
			
		alignments = HTSeq.SAM_Reader(self.args.in_handle)
		for aln in alignments:
			barcode = aln.read.name[-8:]
			aln.read.name = aln.read.name[:-8]
			aln.optional_fields = aln.optional_fields + [('BZ', barcode)]
			self.args.out_handle.write(aln.get_sam_line() + '\n')
	
	def CountBarcodes(self):
		__attr_list = {}
		__list = {}
		__barcodes = {}
		__gff = HTSeq.GFF_Reader(sys.argv[4], end_included=True)
		for __line in __gff:
			if __line.type == 'miRNA_primary_transcript':
				__attr = str(__line.attr['Name']).lower()
				#print __attr
				if __attr not in __list:
					__list[__attr] = {'no_barcode':0}
		__alignments = HTSeq.SAM_Reader(self.args.in_handle)
		for __aln in __alignments:
			__attr = __aln.optional_field('XF')
			__barcode = __aln.optional_field('BZ')
			if __attr in __list:
				if  __barcode not in __list[__attr]:
					__list[__attr][__barcode] = 1
				else:
					__list[__attr][__barcode] += 1
		
		writer = csv.writer(self.args.out_handle, delimiter = '\t')
		for attrkey in __list:
			for bckey in __list[attrkey]:
				writer.writerow([attrkey, bckey, __list[attrkey][bckey]])

class Plots:
	def __init__(self, name):
		self.name = name
		self.sets = []
		self.ylabel = 'y'
		self.xlabel = 'x'
		self.title = 'x'
	    
	#plots a histogram of 'sets' and saves to 'out_handle'
	def histogram(self):
		pylab.hist(self.sets, bins=1)
		xpylab.xlabel(self.xlabel)
		pylab.ylabel(self.ylabel)
		pylab.title(self.title)
		pylab.savefig(self.out_handle)

class __main__:
	def __init__(self):
		#self.test = sys.argv[1]
		#self.in_handle = open(sys.argv[1], "rU")
		#self.out_handle = open(sys.argv[2], "w")
		return 'main'

	if sys.argv[1] == 'fq_barcoding':
		file = BarcodedFq()
		file.MoveBarcode()
		file.GetLengths()
	
	elif sys.argv[1] == 'sam_barcoding':
		file = BarcodedSAM()
		file.MoveBarcode()
	
	elif sys.argv[1] == 'sam_annotation':
		file = BarcodedSAM()
		file.CountBarcodes()
		
	elif sys.argv[1] == 'fq_barcodes':
		file = FastqFile()
		file.CountBarcodes()
	
if __name__ == __main__:
	print 'sdsfasdf'	

#move_barcode(in_handle, out_handle)

#in_handle = open(sys.argv[1], "rU")
#out_handle = open(sys.argv[2], "w")