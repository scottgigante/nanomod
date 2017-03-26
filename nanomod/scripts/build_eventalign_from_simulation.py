################################################################################
#                                                                              #
# This file is part of Nanomod.                                                #
#                                                                              #
# Nanomod is free software: you can redistribute it and/or modify              #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
# Nanomod is distributed in the hope that it will be useful,                   #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with Nanomod.  If not, see <http://www.gnu.org/licenses/>.             #
#                                                                              #
# build_eventalign_from_simulation.py: run a simulation of ONT reads and write #
# the corresponding event and sequence data into an eventalign file.           #
#                                                                              #
# TODO: add options for 1D, 1D template only, 2D                               #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

"""python build_eventalign_from_simulation.py -g ../../data/ecoli_k12.fasta -m ../../models/template_median68pA.model -r ../../data/simulation/canonical/run1 -o ../../data/simulation.canonical --model-seed 1 --num-reads 30000; python build_eventalign_from_simulation.py -g ../../data/ecoli_k12.fasta -m ../../models/template_median68pA.model -r ../../data/simulation/modified/run1 -o ../../data/simulation.modified --model-seed 1 --methyl 1 --num-reads 30000; python build_eventalign_from_simulation.py -g ../../data/ecoli_k12.fasta -m ../../models/template_median68pA.model -r ../../data/simulation/test/run1 -o ../../data/simulation.test --model-seed 1 --methyl 0.3 --num-reads 30000
 """

from __future__ import print_function
import argparse
import Bio
from Bio import SeqIO, Seq, SeqRecord
import os
import sys
import numpy
import pysam
import h5py
import itertools
import csv
from math import sqrt
from multiprocessing import Pool, cpu_count
from difflib import ndiff
from copy import deepcopy
import traceback
from shutil import copyfile

import warnings

_HERTZ = 4000.0
_TEMPLATE_FILE = "../../models/empty_template.fast5"

warnings.simplefilter("error")

# replace a pattern only if at the start of a string
# http://code.activestate.com/recipes/577252-lreplace-and-rreplace-replace-the-beginning-and-en/
def lreplace(pattern, sub, string):
    return sub + string[len(pattern):] if string.startswith(pattern) else string

# replace a pattern only if at the end of a string
def rreplace(pattern, sub, string):
    return string[:-len(pattern)] + sub if string.endswith(pattern) else string

# undo base modification to a sequence
# assume positive site recognition at either end if possible
#
# @args seq The sequence to be unmodified
# @return The unmodified sequence
def unmodifySeq(seq, sequenceMotif):
	if type(seq) == Seq.Seq:
		returnBioSeq = True
		seq = str(seq)
	else:
		returnBioSeq = False

	# we have to check partial matches at either end
	pattern = sequenceMotif[1]
	sub = sequenceMotif[0]
	for i in range(1, len(pattern)):
		seq = lreplace(pattern[i:], sub[i:], seq)
		seq = rreplace(pattern[:i], sub[:i], seq)
	# now replace in the middle
	seq = seq.replace(pattern,sub)
	return Seq.Seq(seq) if returnBioSeq else seq

def parseArgs(argv):

	#command line options
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--genome", dest="refFileName",
			default="data/ecoli_k12.fasta", required=True)
	parser.add_argument("-m", "--model", dest="modelFileName",
			default="models/template_median68pA.model")
	parser.add_argument("-r", "--reads", dest="readsFileName", required=True)
	parser.add_argument("-k", "--kmer", default=6)
	parser.add_argument("-o", "--out-prefix", dest="outPrefix",
			required=True)
	parser.add_argument("-t", "--threads", type=int, default=cpu_count(),
			dest="numThreads")
	parser.add_argument("-n", "--num-reads", default=100, dest="numReads")
	parser.add_argument("--methyl", default=0, type=float)
	parser.add_argument("--bases", default='ACGTM')
	parser.add_argument("-s", "--sequence-motif", dest="sequenceMotif",
			nargs=2, default=["CG","MG"], metavar="CANONICAL MODIFIED",
			help=("Motif of canonical and modified site (e.g., -s CG MG) "
			"(required)"))
	parser.add_argument("--model-seed", default=None, type=int, dest="modelSeed")
	parser.add_argument("--data-seed", default=None, type=int, dest="dataSeed")

	#parse command line options
	options = parser.parse_args()

	return options

def loadRef(fasta):
	refs = dict()
	handle = open(fasta, "rU")
	for record in SeqIO.parse(handle, "fasta"):
		#print "%s with %d bases\n" % (record.id, len(record.seq))
		refs[record.id] = str(record.seq)
	handle.close()
	return refs

def loadModel(modelFile, options):
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

	# create methylated model
	levelMeanOffsetMean = numpy.random.normal(0, 2, 6)
	levelMeanOffsetSd = numpy.random.wald(0.8, 1, 6)
	levelSdOffsetMean = numpy.random.normal(0.1, 0.05, 6)
	levelSdOffsetSd = numpy.random.wald(0.01, 0.1, 6)
	sdMeanOffsetMean = numpy.random.normal(0, 0.05, 6)
	sdMeanOffsetSd = numpy.random.wald(0.01, 0.1, 6)
	sdSdOffsetMean = numpy.random.normal(0.1, 0.05, 6)
	sdSdOffsetSd = numpy.random.wald(0.01, 0.1, 6)
	labels = restrictMethylToCpG(itertools.product(*[options.bases]*options.kmer))
	for l in labels:
		if l not in model:
			unmodified = unmodifySeq(l, options.sequenceMotif)
			model[l] = deepcopy(model[unmodified])
			diff = [i for i in list(ndiff(l, unmodified)) if not i[0] == '+' ]
			locs = [i for i in range(len(diff)) if diff[i][0] == '-']
			for i in locs:
				model[l][0] += numpy.random.normal(levelMeanOffsetMean[i],
						levelMeanOffsetSd[i])
				model[l][1] += numpy.random.normal(levelSdOffsetMean[i],
						levelSdOffsetSd[i])
				model[l][2] += numpy.random.normal(sdMeanOffsetMean[i],
						sdMeanOffsetSd[i])
				model[l][3] += numpy.random.normal(sdSdOffsetMean[i],
						sdSdOffsetSd[i])
	return model

def restrictMethylToCpG(labels):
	# manually remove impossible labels - M only in MG sites
	labels = numpy.array(list(labels))
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
	return [''.join(i) for i in numpy.delete(labels, removeLabels, axis=0).tolist()]

def methylateSeq(seq):
	# simply replacing CG with MG leaves out sequences ending in C that could be methylated
	# as a result, we can't use the last base of the sequence effectively
	return seq.replace("CG","MG")

def methylateChromosome(ref, options):
	ref = ref.split(options.sequenceMotif[0])
	newRef = []
	for seq in ref:
		newRef.append(seq)
		if numpy.random.random() < options.methyl:
			# TODO: include CGIs?
			newRef.append(options.sequenceMotif[1])
		else:
			newRef.append(options.sequenceMotif[0])
	return ''.join(newRef)

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
	if i % (max(int(options.numReads) / 78, 1)) == 0:
		print('=',end='')
		sys.stdout.flush()
	# create a read
	read_name = os.path.abspath(options.readsFileName + "_read" +  str(i) + ".fast5")
	copyfile(_TEMPLATE_FILE, read_name)
	strand = "t"
	read_start = numpy.random.randint(chromosomeLength)
	read_length = int(numpy.random.gamma(1,2000))
	read_end = min(read_start + read_length, chromosomeLength - options.kmer)
	raw = []
	events = []
	basecall = []
	read_squiggle = []
	skips = 0
	for j in range(read_end - read_start):
		# create an event
		model_kmer = ref[chromosome][read_start+j:read_start+j+options.kmer]
		model_line = model[model_kmer]
		ref_kmer = unmodifySeq(model_kmer, options.sequenceMotif)
		# level_mean	level_stdv	sd_mean	sd_stdv	weight
		level = numpy.random.normal(model_line[0], model_line[1])
		# inverse Gaussian
		# var = mu^3/lambda => lambda = mu^3/sigma^2
		stdv = numpy.random.wald(model_line[2],model_line[2]**3/(model_line[3]**2))
		num_samples = int(numpy.random.gamma(shape=1.7, scale=55))
		if num_samples == 0:
			# no samples - skip!
			skips += 1
			continue

		samples = numpy.random.normal(level,stdv,num_samples)
		sample_mean = numpy.mean(samples)
		sample_stdv = sqrt(numpy.var(samples))
		sample_length = num_samples / _HERTZ
		squiggle_start = len(raw)
		squiggle_end = squiggle_start + len(samples) - 1
		# contig	position	reference_kmer	read_name	strand	event_index	event_level_mean	event_stdv	event_length	model_kmer	model_mean	model_stdv
		read_squiggle.append([chromosome, read_start + j, ref_kmer, read_name, strand, j - skips, sample_mean, sample_stdv, sample_length, model_kmer, level, stdv])
		raw.extend(samples)
		events.append((sample_mean, sample_stdv, squiggle_start, num_samples))
		basecall.append((sample_mean, squiggle_start/_HERTZ, sample_stdv, sample_length, ref_kmer, 0, 1, 1, ref_kmer, 0, 0, 0, 0, 0))
		assert(len(basecall) == j-skips+1)
		assert(len(basecall) == len(read_squiggle))

	# output .fast5 read
	with h5py.File(read_name, "r+") as fast5:
		raw_group = fast5.get("Raw/Reads")
		raw_group.move("Read_0","Read_{}".format(i))
		raw_dset = raw_group.get("Read_{}/Signal".format(i))
		raw = numpy.array(raw, dtype=raw_dset.dtype)
		raw_dset.resize((raw.shape[0],))
		raw_dset[()] = raw

		event_group = fast5.get("Analyses/EventDetection_000/Reads")
		event_group.move("Read_0","Read_{}".format(i))
		event_dset = event_group.get("Read_{}/Events".format(i))
		events = numpy.array(events, dtype=event_dset.dtype)
		event_dset.resize((events.shape[0],))
		event_dset[()] = events

		basecall_dset = fast5.get("Analyses/Basecall_1D_000/BaseCalled_template/Events")
		basecall = numpy.array(basecall, dtype=basecall_dset.value.dtype)
		basecall_dset.resize((basecall.shape[0],))
		basecall_dset[()] = basecall

	return read_squiggle, ref[chromosome][read_start:read_end+options.kmer]

def multi_run_wrapper(args):
	try:
		return createFast5(*args)
	except Exception as e:
		print('Caught exception in worker thread:')
		traceback.print_exc()
		print()
		raise e

def extractReads(options):
	numpy.random.seed(options.modelSeed)
	print('loading...')
	numReads = int(options.numReads)
	options.readsFileName = os.path.abspath(options.readsFileName)
	ref = loadRef(options.refFileName)
	model = loadModel(options.modelFileName, options)
	for chromosome in ref:
		# create a set of numReads reads on each chromosome
		print('building chromosome ' + str(chromosome) + '...')
		chromosomeLength = len(ref[chromosome])
		ref[chromosome] = methylateChromosome(ref[chromosome], options)
		# TODO: what about reverse reads?
		numpy.random.seed(options.dataSeed)
    	pool = Pool(options.numThreads)
    	print('[',end='')
    	sys.stdout.flush()
    	squiggle = pool.map(multi_run_wrapper, [[options, i, chromosomeLength, ref, chromosome, model] for i in range(int(options.numReads))])
    	print(']')
	print('writing...')
	# output eventalign
	with open("{}.eventalign".format(options.outPrefix), "w") as eventalign, open("{}.fasta".format(options.outPrefix), 'w') as fasta:
		writer = csv.writer(eventalign, delimiter='\t')
		writer.writerow(["contig", "position", "reference_kmer", "read_name",
				"strand", "event_index", "event_level_mean", "event_stdv",
				"event_length", "model_kmer", "model_mean", "model_stdv"])
		for read, seq in squiggle:
			try:
				record = SeqRecord.SeqRecord(Seq.Seq(seq), id=read[0][3], description="{} {}".format(options.readsFileName, read[0][3]))
				SeqIO.write(record, fasta, "fasta")
				for line in read:
					writer.writerow(line)
			except IndexError:
				# read empty
				pass

def run(argv):
	options = parseArgs(argv)
	extractReads(options)

if __name__ == "__main__":
	run(sys.argv)
