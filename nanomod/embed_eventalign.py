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
# TODO: deal with secondary / auxiliary reads better                           #
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
import logging
from numpy.lib.recfunctions import append_fields, drop_fields
from multiprocessing import Pool, cpu_count
from shutil import copyfile

from seq_tools import *
from utils import callSubProcess, makeDir, multiprocessWrapper, preventOverwrite
from check_skip_stay_prob import getSkipStayConstraints

# get the path of the output file corresponding to a fast5 file
#
# @args fast5 basename of the fast5 file
# @args options Namespace object from argparse
# @return path to output fast5 file
def getOutfile(fast5, options):
	return os.path.join(options.outPrefix, fast5)

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
		chromosome, forward, numSkips, numStays, readLength, kmer, genome):
	
	# all events are considered good emissions - bad idea?
	events = append_fields(events, 'good_emission', events["kmer"] != 'X'*kmer)
	# remove move field - nanonet adds it
	events = drop_fields(events, "move")
	
	# median normalise
	# TODO: we do this twice - make this into a routine.
	# TODO: maybe swap to nanoraw instead of eventalign
	if not options.noNormalize:
		med = np.median(np.concatenate([np.repeat(events["mean"][i], int(events["length"][i] * 4000)) for i in range(events.shape[0])])) # TODO: pull sample rate from uniqueglobalkey
		median_deviation = [events["mean"][i] - med for i in range(events.shape[0])]
		MAD = 1/np.sum(events["length"]) * sum(median_deviation[i] * events["length"][i] for i in range(events.shape[0]))
		events["mean"] = np.array([median_deviation[i] / MAD for i in range(events.shape[0])])
	
	# remove unlabelled events 
	# long stretches of unlabelled events severely decrease accuracy
	events = events[events["kmer"] != 'X'*kmer]
	
	# TODO: if the read fails QC, should we just move on?
	numEvents = events.shape[0]
	skipProb = float(numSkips)/numEvents
	stayProb = float(numStays)/numEvents
	stepProb = 1-skipProb-stayProb
	passQuality = (skipProb < options.constraints['maxSkips'] and 
			stayProb < options.constraints['maxStays'] and 
			stepProb > options.constraints['minSteps'] and 
			readLength >= options.readLength)
	
	filename = fast5Path.split('/')[-1]
	newPath = getOutfile(filename, options)
	copyfile(fast5Path, newPath)
	logging.debug("Saving {}".format(filename))
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
		fasta = genome[chromosome]['record'].format("fasta")
		fastaGroup.create_dataset("Fasta", 
				data=np.array(fasta, dtype='|S{}'.format(len(fasta))))
		
		# create attrs
		summaryGroup = alignGroup.create_group("Summary")
		attrsGroup = summaryGroup.create_group("genome_mapping_template")
		attrs = attrsGroup.attrs
		attrs.create("genome", chromosome)
		
	return passQuality

# step through eventalign data corresponding to a fast5 file and label events
#
# TODO: can we determine kmer length earlier?
#
# @args options Namespace object from argparse
# @args idx headers index of eventalign file
# @args fast5Path path to original fast5 file
# @return two-element array containing basename of fast5 file and string reference to either training or validation dataset
def processRead(options, idx, fast5Path, genome, modified, kmer):
	
	fast5File = fast5Path.split('/')[-1]
	try:
		eventalign = np.load(os.path.join(options.tempDir, fast5File + ".npy"))
	except IOError:
		loggoing.warning("Failed to load {}.npy".format(os.path.join(options.tempDir, fast5File)))
		# why on earth?!
		return [fast5File, False]
	skip = False
	
	# open fast5 file
	try:
		fast5 = h5py.File(fast5Path, 'r')
		# r9 for now
		events = fast5.get("Analyses/Basecall_1D_000/BaseCalled_template/Events").value
		events = append_fields(events, 'kmer', ['X'*kmer] * events.shape[0])
		events = append_fields(events, 'seq_pos', [-1] * events.shape[0]) # is there a better option here? keep counting?
	except IOError:
		logging.warning("Failed to read {}".format(fast5Path))
		return [fast5File, False]
	except TypeError:
		logging.warning("Data corrupted in {}".format(fast5Path))
		return [fast5File, False]

	try:
		readLength = fast5.get("Analyses/Basecall_1D_000/Summary/basecall_1d_template").attrs["sequence_length"]
	except AttributeError:
		# attributes not found
		readLength = options.readLength
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
		
		seq = getKmer(genome, chromosome, current_ref_index, kmer, forward)
		
		if last_ref_index == current_ref_index:
			numStays += 1
		elif abs(last_ref_index - current_ref_index) > 1:
			# count the size of the skip 
			# skipping two bases at once counts as two skips
			numSkips += abs(last_ref_index - current_ref_index) - 1
		last_ref_index = current_ref_index
		last_seq_pos = seq_pos
		
		events["seq_pos"][int(line[idx['event_index']])] = seq_pos
		events["kmer"][int(line[idx['event_index']])] = seq
	
	fast5.close()
	try:
		pass_quality = writeFast5(options, events, fast5Path, initial_ref_index, 
				last_ref_index, chromosome, forward, numSkips, numStays, 
				readLength, kmer, genome)
	except Exception as e:
		logging.warning("Failed to save {}.".format(fast5File))
		print e.message
		return [fast5File, False]
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
		
		kmer=None
		current_read_name = ""
		skip=True
		tmp=None
		filenames = set()
		premadeFilenames = set()
		n=0
		
		for line in reader:
			
			if line[idx['read_name']] == current_read_name and skip:
				continue
			elif line[idx['read_name']] != current_read_name:
				if kmer is None:
					# need to initialise kmer length
					kmer = len(line[idx['ref_kmer']])
				if tmp is not None and len(tmp) > 0:
					try:
						np.save(filename,tmp)
						logging.debug("Saving {}.npy".format(filename))
						n += 1
					except IOError as e:
						logging.warning("Failed to write {}.npy: {}".format(filename, e))
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
					logging.warning(current_read_name + " not found.")
					continue
					
				fast5Name = fast5Path.split('/')[-1]
				if line[idx['strand']] != 't':
					# not template, skip
					skip=True
					continue
				
				outfile = getOutfile(fast5Name, options)
				filename = os.path.join(options.tempDir, fast5Name)
				if (not options.force) and os.path.isfile(outfile):
					logging.debug(outfile + " already exists. Use --force to recompute.")
					skip=True
					n += 1
					premadeFilenames.add(fast5Name)
					continue
				elif (not options.force) and os.path.isfile(filename + ".npy"):
					logging.debug(("{}.npy already exists. Use --force to " 
							"recompute.").format(filename))
					skip=True
					filenames.add(fast5Path)
					assert(tmp == [])
					continue
				
				filenames.add(fast5Path)
				skip=False
			
			assert(skip==False)
			tmp.append(line)
		
		# last one gets missed
		# TODO: should we exclude fail reads earlier to save time?
		if tmp is not None and len(tmp) > 0:
			np.save(filename,tmp)
			logging.debug("Saving {}.npy".format(filename))
		
	return filenames, idx, premadeFilenames, kmer

# check quality of premade fast5 files
# @args options Namespace object from argparse
# @args filename Name of premade fast5 file
def checkPremade(options, filename):
	try:
		with h5py.File(getOutfile(filename, options), 'r') as fh:
			attrs = fh.get(("Analyses/Basecall_1D_000/Summary/basecall_1d"
					"_template")).attrs
			numEvents = float(attrs["called_events"])
			skipProb = attrs["num_skips"]/numEvents
			stayProb = attrs["num_stays"]/numEvents
			stepProb = 1-skipProb-stayProb
			readLength = attrs["sequence_length"]
			passQuality = (skipProb < options.constraints['maxSkips'] and 
					stayProb < options.constraints['maxStays'] and 
					stepProb > options.constraints['minSteps'] and
					readLength >= options.readLength)
	except IOError:
		logging.warning("Failed to read {}".format(filename))
		passQuality = False
	except AttributeError:
		# can't find attributes
		logging.info("Attributes missing on {}".format(filename))
		passQuality = True
	return [filename, passQuality]

# write text files for training and validation consisting of fast5 basenames
#
# @args options Namespace object from argparse
# @args trainData Dataset consisting or two-element arrays of filename and dataset to which file should be added (either "train" or "val")
# @return None
def writeTrainfiles(options, trainData, outPrefix):
	trainFilename = outPrefix + ".train.txt"
	valFilename = outPrefix + ".val.txt"
	
	with open(trainFilename, 'w') as trainFile, open(valFilename, 'w') as valFile, open(trainFilename + ".small", 'w') as smallTrainFile,	open(valFilename + ".small", 'w') as smallValFile:

		# write trainFile / valFile headers
		trainFile.write("#filename\n")
		valFile.write("#filename\n")
		smallTrainFile.write("#filename\n")
		smallValFile.write("#filename\n")
		
		for filename, pass_quality in trainData:
			if os.path.isfile(getOutfile(filename, options)):
				if np.random.rand() > options.valFraction:
					logging.debug("Training set: {}".format(filename))
					trainFile.write(filename + "\n")
					if pass_quality:
						smallTrainFile.write(filename + "\n")
				else:
					logging.debug("Validation set: {}".format(filename))
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
def embedEventalign(options, fasta, eventalign, reads, outPrefix, modified):
	output = "{}.train.txt.small".format(outPrefix)
	if preventOverwrite(output, options):
		return 1
	makeDir(options.outPrefix)
	
	logging.info("Loading fast5 names and references...")
	refs = loadRef(fasta)
	genome = loadGenome(options, modified)
	
	if "random" in options.selectMode:
		# no point populating constrained files, we will overwrite later
		options.constraints = {'maxSkips' : 0,
			'maxStays' : 0, 
			'minSteps' : 1 }
	else:
		logging.info("Calculating skip/stay count constraints...")
		options.constraints = getSkipStayConstraints(reads, 
				options.dataFraction, options.selectMode, options.readLength)
		logging.debug(str(options.constraints))
	
	logging.info("Splitting eventalign into separate files...")
	filenames, idx, premadeFilenames, kmer = writeTempFiles(options, eventalign, refs)
	pool = Pool(options.threads)
	
	logging.info("Embedding labels into {} fast5 files...".format(len(filenames)))
	trainData = pool.map(processReadWrapper, 
			[[options, idx, i, genome, modified, kmer] for i in filenames])
	
	logging.info("Adding data for {} premade fast5 files...".format(len(premadeFilenames)))
	trainData.extend(pool.map(checkPremadeWrapper, 
			[[options, i] for i in premadeFilenames]))
	
	logging.info("Saving datasets for {} total fast5 files...".format(len(trainData)))
	writeTrainfiles(options, trainData, outPrefix)
