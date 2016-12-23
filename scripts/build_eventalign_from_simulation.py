#!/usr/bin/env python

# from Rob Egan's fork of Deepnano
# https://bitbucket.org/robegan21/deepnano/commits/661986be77ac3d6b5178e947a3bab965e5a56f2b
# TODO: add options for 1D, 1D template only, 2D

from __future__ import print_function
from optparse import OptionParser
import Bio
from Bio import SeqIO
import os
import sys
import numpy
import pysam
import h5py
import netcdf_helpers
import itertools
import csv
from math import sqrt
from multiprocessing import Pool, cpu_count

import warnings

warnings.simplefilter("error")


def parseArgs(argv):
	
	#command line options
	parser = OptionParser()
	parser.add_option("-g", "--genome", dest="refFileName", default="data/ecoli_k12.fasta")
	parser.add_option("-m", "--model", dest="modelFileName", default="data/r7.3_e6_70bps_6mer_template_median68pA.model")
	parser.add_option("-r", "--reads", dest="readsFileName")
	parser.add_option("-k", "--kmer", dest="kmer", default=6)
	parser.add_option("-e", "--eventalign", dest="eventalignFileName")
	parser.add_option("-t", "--threads", type=int, default=cpu_count(), dest="numThreads")
	parser.add_option("-n", "--num-reads", default=10, dest="numReads")
	parser.add_option("--methyl", default=False, action="store_true",dest="methyl")
	
	#parse command line options
	(options, args) = parser.parse_args()
	
	# Making sure all mandatory options appeared.
	mandatories = ['refFileName', "eventalignFileName", "readsFileName"]
	for m in mandatories:
		if not options.__dict__[m]:
			print("mandatory option %s is missing" % m)
			parser.print_help()
			exit(-1)
	
	return (options, args)

def loadRef(fasta):
	refs = dict()
	handle = open(fasta, "rU")
	for record in SeqIO.parse(handle, "fasta"):
		#print "%s with %d bases\n" % (record.id, len(record.seq))
		refs[record.id] = str(record.seq)
	handle.close()
	return refs
	
def loadModel(modelFile):
	model = dict()
	with open(modelFile, "r") as tsv:
		init = False
		for line in tsv:
			line = line.strip().split("\t")
			if not init:
				if line[0] == "kmer":
					# header line
					init = True
			else:
				# kmer	level_mean	level_stdv	sd_mean	sd_stdv	weight
				model[line[0]] = [float(line[i]) for i in range(1,len(line))]
	
	return model

def restrictMethylToCpG(labels):
	# manually remove impossible labels - M only in MG sites
	removeLabels=[]
	i = 0
	for label in labels:
		methyl = False
		for c in label:
			if methyl and c != 'G':
				removeLabels.append(i)
				break
			elif c == 'M':
				methyl = True
			else:
				methyl = False
		i += 1
	return list(numpy.delete(labels, removeLabels))

def methylateSeq(seq):
	# simply replacing CG with MG leaves out sequences ending in C that could be methylated
	# as a result, we can't use the last base of the sequence effectively	
	return seq.replace("CG","MG")
	
def reverseComplement(seq):
	# reverse
	seq = seq[::-1]	
	# complement
	# from http://stackoverflow.com/questions/6116978/python-replace-multiple-strings
	# in order to avoid double changes, use symbols first
	repls = ('A', '1'), ('T', '2'), ('C', '3'), ('G', '4')
	seq = reduce(lambda a, kv: a.replace(*kv), repls, seq)
	repls = ('1', 'T'), ('2', 'A'), ('3', 'G'), ('4', 'C')
	return reduce(lambda a, kv: a.replace(*kv), repls, seq)

def normalize(inp, normalizeMeans, normalizeStdvs):
	assert(len(inp) == len(normalizeMeans))
	inp = [(inp[i] - normalizeMeans[i])/normalizeStdvs[i] for i in range(len(inp))]
	return inp

def createFast5(options, i, chromosomeLength, ref, chromosome, model):
	# print progress
	if i % (int(options.numReads) / 78) == 0:
		print('=',end='')
		sys.stdout.flush()
	# create a read
	read_name = options.readsFileName + "_read" +  str(i) + ".fast5"
	strand = "t"
	read_start = numpy.random.randint(chromosomeLength)
	read_length = int(numpy.random.gamma(1,2000))
	read_end = min(read_start + read_length, chromosomeLength - options.kmer)
	raw = []
	events = []
	read_squiggle = []
	for j in range(read_end - read_start):
		# create an event
		ref_kmer = ref[chromosome][read_start+j:read_start+j+options.kmer]
		model_line = model[ref_kmer]
		# level_mean	level_stdv	sd_mean	sd_stdv	weight
		level = numpy.random.normal(model_line[0], model_line[1])
		# inverse Gaussian
		# var = mu^3/lambda => lambda = mu^3/sigma^2
		stdv = numpy.random.wald(model_line[2],model_line[2]**3/(model_line[3]**2))
		num_samples = int(numpy.random.gamma(shape=1.7, scale=55))
		if num_samples == 0:
			# no samples - skip!
			continue
		
		samples = numpy.random.normal(level,stdv,num_samples)
		sample_mean = numpy.mean(samples)
		sample_stdv = sqrt(numpy.var(samples))
		sample_length = num_samples / 3012.0
		squiggle_start = len(raw)
		squiggle_end = squiggle_start + len(samples) - 1
		# contig	position	reference_kmer	read_name	strand	event_index	event_level_mean	event_stdv	event_length	model_kmer	model_mean	model_stdv
		read_squiggle.append([chromosome, read_start + j, ref_kmer, read_name, strand, j, sample_mean, sample_stdv, sample_length, ref_kmer, level, stdv])
		raw.extend(samples)
		events.append((sample_mean, sample_stdv, squiggle_start, num_samples))
	
	# output .fast5 read
	with h5py.File(read_name, "w") as fast5:
		raw_group = fast5.create_group("Raw/Reads/Read_" + str(i))
		raw_group.create_dataset("Signal", data=raw)
		event_group = fast5.create_group("Analyses/EventDetection_000/Reads/Read_" + str(i))
		events = numpy.array(events, dtype=[('mean','float64'),('stdv','float64'),('start','int64'),('length','int32')])
		event_group.create_dataset("Events", data=events)
	return read_squiggle

def multi_run_wrapper(args):
   return createFast5(*args)

def extractReads(options):
	print('loading...')
	numReads = int(options.numReads)
	ref = loadRef(options.refFileName)
	model = loadModel(options.modelFileName)
	for chromosome in ref:
		# create a set of numReads reads on each chromosome
		print('building chromosome ' + str(chromosome) + '...')
		chromosomeLength = len(ref[chromosome])
		# could we do this in parallel?
    	pool = Pool(options.numThreads)
    	print('[',end='')
    	sys.stdout.flush()
    	squiggle = pool.map(multi_run_wrapper, [[options, i, chromosomeLength, ref, chromosome, model] for i in range(int(options.numReads))])
    	print(']')
	print('writing...')
	# output eventalign
	with open(options.eventalignFileName, "w") as eventalign:
		writer = csv.writer(eventalign, delimiter='\t')
		writer.writerow(["contig","position","reference_kmer","read_name","strand","event_index","event_level_mean","event_stdv","event_length","model_kmer","model_mean","model_stdv"])
		for read in squiggle:
			for line in read:
				writer.writerow(line)

def run(argv):
	(options, args) = parseArgs(argv)
	extractReads(options)

if __name__ == "__main__":
	run(sys.argv)
