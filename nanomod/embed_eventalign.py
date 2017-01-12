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
# embed_eventalign.py: use eventalign data to label individual fast5 reads     #
# in order to feed labelled fast5 files into neural network training.          #
#                                                                              #
# TODO: remove Numpy FutureWarning from selecting fields in structured array.  #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

from Bio import SeqIO
import os
import pysam
import h5py
import numpy as np
import csv
from numpy.lib.recfunctions import append_fields, drop_fields
from multiprocessing import Pool, cpu_count
from shutil import copyfile

from seq_tools import *
from utils import log, makeDir, multiprocessWrapper
from check_skip_stay_prob import getSkipStayConstraints

__outdir_name__ = "nanomod" # subdirectory to be created in the output dir

# get the path of the output file corresponding to a fast5 file
#
# @args fast5 basename of the fast5 file
# @args options Namespace object from argparse
# @return path to output fast5 file
def getOutfile(fast5, options):
	return os.path.join(options.outPrefix, __outdir_name__, fast5)

# write a labelled fast5 file
#
# @args options Namespace object from argparse
# @args events Structured numpy array including event labels
# @args fast5Path path to original fast5 file
# @args initial_ref_index starting position of read on reference genome
# @args last_ref_index ending position of read on reference genome
# @args chromosome Name of the chromosome to which fast5 is aligned
# @args forward Boolean value representing forward or reverse read along genome
# @args numSkips Number of skipped bases on genome
# @args numStays Number of repeated events without moving along genome
# @args kmer Length of kmer labels
# @return String representing whether fast5 file should be train or val data
def writeFast5(options, events, fast5Path, initial_ref_index, last_ref_index, 
		chromosome, forward, numSkips, numStays, kmer):
	
	# all events are considered good emissions - bad idea?
	events = append_fields(events, 'good_emission', events["kmer"] != 'X'*kmer)
	# remove move field - nanonet adds it
	events = drop_fields(events, "move")
	
	# remove unlabelled events -  TODO: is this helpful?
	events = events[events["kmer"] != 'X'*kmer]
	
	# TODO: check if this read passes or fails
	numEvents = events.shape[0]
	skipProb = float(numSkips)/numEvents
	stayProb = float(numStays)/numEvents
	stepProb = 1-skipProb-stayProb
	
	filename = fast5Path.split('/')[-1]
	newPath = getOutfile(filename, options)
	copyfile(fast5Path, newPath)
	with h5py.File(newPath,'r+') as fast5:
		analysesGroup = fast5.get("Analyses")
		alignToRefGroup = analysesGroup.create_group("AlignToRef")
		
		# create events
		eventsGroup = alignToRefGroup.create_group(("CurrentSpaceMapped"
				"_template"))
		eventsGroup.create_dataset("Events", data=events)
		
		# create attrs
		summaryGroup = alignToRefGroup.create_group("Summary")
		attrsGroup = summaryGroup.create_group("current_space_map_template")
		attrs = attrsGroup.attrs
		attrs.create("genome_start", initial_ref_index)
		attrs.create("genome_end", last_ref_index)
		attrs.create("sequence_length", last_ref_index - initial_ref_index)
		attrs.create("ref_name", chromosome)
		attrs.create("direction", "+" if forward else "-")
		attrs.create("num_skips", numSkips)
		attrs.create("num_stays", numStays)
		attrs.create("num_events", numEvents)
		attrs.create("skip_prob", skipProb)
		attrs.create("stay_prob", stayProb)
		attrs.create("step_prob", stepProb)
		
		# create alignment group
		alignGroup = analysesGroup.create_group("Alignment")
		fastaGroup = alignGroup.create_group("Aligned_template")
		# our current eventalign files don't have the same chromosome name!???
		# TODO: don't do this pls
		fasta = options.genome["Chromosome"].format("fasta")
		fastaGroup.create_dataset("Fasta", 
				data=options.genome["Chromosome"].format("fasta"))
		
		# create attrs
		summaryGroup = alignGroup.create_group("Summary")
		attrsGroup = summaryGroup.create_group("genome_mapping_template")
		attrs = attrsGroup.attrs
		attrs.create("genome", chromosome)
		
	return (skipProb < options.constraints['maxSkips'] and 
			stayProb < options.constraints['maxStays'] and 
			stepProb > options.constraints['minSteps'])

# step through eventalign data corresponding to a fast5 file and label events
#
# TODO: can we determine kmer length earlier?
#
# @args options Namespace object from argparse
# @args idx headers index of eventalign file
# @args fast5Path path to original fast5 file
# @return two-element array containing basename of fast5 file and string reference to either training or validation dataset
def processRead(options, idx, fast5Path):
	
	fast5File = fast5Path.split('/')[-1]
	eventalign = np.load(os.path.join(options.tempDir, fast5File + ".npy"))
	skip = False
	kmer=None # could initialise this earlier
	
	# open fast5 file
	try:
		fast5 = h5py.File(fast5Path, 'r')
	except IOError:
		log("Failed to read {}".format(fast5Path), 0, options)
		return None
		
	# r9 for now
	events = np.array(fast5.get("Analyses/Basecall_1D_000/BaseCalled_template/Events"))
	if kmer is None:
		kmer = len(eventalign[0][idx['ref_kmer']])
	events = append_fields(events, 'kmer', ['X'*kmer] * events.shape[0])
	events = append_fields(events, 'seq_pos', [-1] * events.shape[0]) # is there a better option here? keep counting?
	initial_event_index = int(eventalign[0][idx['event_index']])
	initial_ref_index = int(eventalign[0][idx['ref_pos']])
	last_ref_index = initial_ref_index-1
	chromosome = eventalign[0][idx['contig']]
	start=False
	numSkips = 0
	numStays = 0
	
	for line in eventalign:
		
		if not line[idx['strand']] == 't':
			# finished reading the template
			# now reading complement of the same read - skip
			break
		
		seq = line[idx['ref_kmer']]
		current_ref_index = int(line[idx['ref_pos']])
		
		if not start:
			if int(line[idx['event_index']]) > initial_event_index:
				# forward strand
				forward = True
			else:
				# reverse strand
				forward = False
			start = True
			seq_pos = 0
		else:
			seq_pos_diff = abs(current_ref_index - last_ref_index)
			seq_pos = last_seq_pos + seq_pos_diff
		
		if not forward:
			seq = reverseComplement(seq)
		
		if options.methyl:
			unmethyl_seq = seq
			seq = methylateSeq(seq)
		
		if last_ref_index == current_ref_index:
			numStays += 1
		elif abs(last_ref_index - current_ref_index) > 1:
			# does skipping two bases at once count as 2 skips or one???
			# TODO: count size of skip, take negative into account
			numSkips += abs(last_ref_index - current_ref_index) - 1
		last_ref_index = current_ref_index
		last_seq_pos = seq_pos
		
		events["seq_pos"][int(line[idx['event_index']])] = seq_pos
		events["kmer"][int(line[idx['event_index']])] = seq
	
	fast5.close()
	pass_quality = writeFast5(options, events, fast5Path, initial_ref_index, last_ref_index, chromosome, forward, numSkips, numStays, kmer)
	return [fast5File, pass_quality]

# write temporary .npy files to split eventalign tsv file into chunks 
# corresponding to one fast5 file each - this allows multiprocessing later on
#
# @args options Namespace object from argparse
# @args eventalign filename of eventalign tsv
# @args refs dictionary linking read names and fast5 file paths
# @return filenames a list of paths to fast5 files which have been processed
# @return idx a dictionary linking headers with positions in the tsv
# @return premadeFilenames a list of basenames of fast5 files which had already been processed prior to this call to Nanomod
def writeTempFiles(options, eventalign, refs):
	if not os.path.exists(options.tempDir):
		os.makedirs(options.tempDir)
	
	with open(eventalign) as tsv:
		
		# open file
		reader = csv.reader(tsv, delimiter="\t")
		headers = reader.next()
		idx={}
		idx['contig'] = headers.index("contig")
		idx['ref_pos'] = headers.index("position")
		idx['ref_kmer'] = headers.index("reference_kmer")
		idx['model_kmer'] = headers.index("model_kmer")
		idx['read_name'] = headers.index("read_name")
		idx['strand'] = headers.index("strand")
		idx['event_index'] = headers.index("event_index")
		idx['mean'] = headers.index("event_level_mean")
		idx['stdv'] = headers.index("event_stdv")
		idx['length'] = headers.index("event_length")
		
		current_read_name = ""
		skip=True
		tmp=None
		filenames = []
		premadeFilenames = []
		n=0
		
		for line in reader:
			
			if line[idx['read_name']] == current_read_name and skip:
				continue
			elif line[idx['read_name']] != current_read_name:
				if tmp is not None and tmp != []:
					try:
						np.save(filename,tmp)
						n += 1
					except IOError as e:
						log("Failed to write {}.npy: {}".format(filename, e), 0, 
								options)
				if options.numReads > 0 and n >= options.numReads:
					# we're done!
					break
				
				# new read, initialise
				tmp = []
				current_read_name = line[idx['read_name']]
				try:
					fast5Path = refs[current_read_name]
				except KeyError:
					skip = True
					log(current_read_name + " not found", 0, options)
					continue
					
				fast5Name = fast5Path.split('/')[-1]
				if line[idx['strand']] != 't':
					# not template, skip
					skip=True
					continue
				
				outfile = getOutfile(fast5Name, options)
				filename = os.path.join(options.tempDir, fast5Name)
				if (not options.force) and os.path.isfile(outfile):
					log(outfile + " already exists. Use --force to recompute.", 
							2, options)
					skip=True
					n += 1
					premadeFilenames.append(fast5Name)
					continue
				elif (not options.force) and os.path.isfile(filename + ".npy"):
					log(("{}.npy already exists. Use --force to " 
							"recompute.").format(filename), 2, options)
					skip=True
					filenames.append(fast5Path)
					assert(tmp == [])
					continue
				
				filenames.append(fast5Path)
				skip=False
			
			assert(skip==False)
			tmp.append(line)
		
		# last one gets missed
		if tmp is not None and tmp != []:
			np.save(filename,tmp)
		
	return filenames, idx, premadeFilenames

# check quality of premade fast5 files
def checkPremade(options, filename):
	with h5py.File(getOutfile(filename, options), 'r') as fh:
		attrs = fh.get(("Analyses/Basecall_1D_000/Summary/basecall_1d"
				"_template")).attrs
		numEvents = float(attrs["called_events"])
		skipProb = attrs["num_skips"]/numEvents
		stayProb = attrs["num_stays"]/numEvents
		stepProb = 1-skipProb-stayProb
	return [filename, (skipProb < options.constraints['maxSkips'] and 
			stayProb < options.constraints['maxStays'] and 
			stepProb > options.constraints['minSteps'])]

# write text files for training and validation consisting of fast5 basenames
#
# @args options Namespace object from argparse
# @args trainData Dataset consisting or two-element arrays of filename and dataset to which file should be added (either "train" or "val")
# @return None
def writeTrainfiles(options, trainData):
	trainFilename = options.outPrefix + ".train.txt"
	valFilename = options.outPrefix + ".val.txt"
	
	with open(trainFilename, 'w') as trainFile, open(valFilename, 'w') as valFile, open(trainFilename + ".small", 'w') as smallTrainFile,	open(valFilename + ".small", 'w') as smallValFile:

		# write trainFile / valFile headers
		trainFile.write("#filename\n")
		valFile.write("#filename\n")
		smallTrainFile.write("#filename\n")
		smallValFile.write("#filename\n")
		
		for filename, pass_quality in trainData:
			if os.path.isfile(getOutfile(filename, options)):
				if np.random.rand() > options.valFraction:
					log("Training set: {}".format(filename), 2, options)
					trainFile.write(filename + "\n")
					if pass_quality:
						smallTrainFile.write(filename + "\n")
				else:
					log("Validation set: {}".format(filename), 2, options)
					valFile.write(filename + "\n")
					if pass_quality:
						smallValFile.write(filename + "\n")

def processReadWrapper(args):
	return multiprocessWrapper(processRead, args)

def checkPremadeWrapper(args):
	return multiprocessWrapper(checkPremade, args)

# main script: embed labels to reference genome into fast5 files
#
# @args options Namespace object from argparse
# @args fasta Fasta filename corresponding to reads (MUST be produced by poretools)
# @args eventalign Filename of eventalign tsv output
# @return None
def embedEventalign(options, fasta, eventalign):
	makeDir(os.path.join(options.outPrefix, __outdir_name__))
	log("Loading fast5 names and references...", 1, options)
	refs = loadRef(fasta)
	options.genome = loadGenome(options.genome)
	log("Calculating skip/stay count constraints...", 1, options)
	options.constraints = getSkipStayConstraints(options.reads, 
			options.dataFraction)
	log("Splitting eventalign into separate files...", 1, options)
	filenames, idx, premadeFilenames = writeTempFiles(options, eventalign, refs)
	pool = Pool(options.threads)
	log("Embedding labels into {} fast5 files...".format(len(filenames)), 1,
			options)
	trainData = pool.map(processReadWrapper, 
			[[options, idx, i] for i in filenames])
	log(("Adding data for {} premade fast5" 
			"files...").format(len(premadeFilenames)), 1, options)
	trainData.extend(pool.map(checkPremadeWrapper, 
			[[options, i] for i in premadeFilenames]))
	log(("Saving datasets for {} total fast5" 
			"files...").format(len(trainData)), 1, options)
	writeTrainfiles(options, trainData)
