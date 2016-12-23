from Bio import SeqIO
import os
import pysam
import h5py
import numpy as np
import csv
from numpy.lib.recfunctions import append_fields
from multiprocessing import Pool, cpu_count
from shutil import copyfile

from seq_tools import *
from utils import log, makeDir

__outdir_name__ = "nanomod"

def getOutfile(fast5, options):
	return os.path.join(options.outPrefix, __outdir_name__, fast5)

def writeFast5(options, events, fast5Path, initial_ref_index, last_ref_index, chromosome, forward, numSkips, numStays, kmer):

	# all events are considered good emissions - bad idea?
	events = append_fields(events, 'good_emission', events["kmer"] != 'X'*kmer)
	# remove move field - nanonet adds it
	names = list(events.dtype.names)
	names.remove("move")
	newEvents = events[names].copy()
	
	filename = fast5Path.split('/')[-1]
	newPath = getOutfile(filename, options)
	copyfile(fast5Path, newPath)
	with h5py.File(newPath,'r+') as fast5:
		analysesGroup = fast5.get("Analyses")
		alignToRefGroup = analysesGroup.create_group("AlignToRef")
		
		# create events
		eventsGroup = alignToRefGroup.create_group("CurrentSpaceMapped_template")
		eventsGroup.create_dataset("Events", data= newEvents)
		
		# create attrs
		summaryGroup = alignToRefGroup.create_group("Summary")
		attrsGroup = summaryGroup.create_group("current_space_map_template")
		attrs = attrsGroup.attrs
		attrs.create("genome_start", initial_ref_index)
		attrs.create("genome_end", last_ref_index)
		attrs.create("ref_name", chromosome)
		attrs.create("direction", "+" if forward else "-")
		attrs.create("num_skips", numSkips)
		attrs.create("num_stays", numStays)
		
		# create alignment group
		alignGroup = analysesGroup.create_group("Alignment")
		fastaGroup = alignGroup.create_group("Aligned_template")
		# our current eventalign files don't have the same chromosome name!???
		fasta = options.genome["Chromosome"].format("fasta")
		fastaGroup.create_dataset("Fasta", data=options.genome["Chromosome"].format("fasta"))
		
		# create attrs
		summaryGroup = alignGroup.create_group("Summary")
		attrsGroup = summaryGroup.create_group("genome_mapping_template")
		attrs = attrsGroup.attrs
		attrs.create("genome", chromosome)
		
	if np.random.rand() < options.valFraction:
		return "val"
	else:
		return "train"

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
		return ["",""]
		
	# r7.3 for now
	events = np.array(fast5.get("Analyses/Basecall_1D_000/BaseCalled_template/Events"))
	if kmer is None:
		kmer = len(eventalign[0][idx['ref_kmer']])
	newEvents = append_fields(events, 'kmer', ['X'*kmer] * events.shape[0])
	newEvents = append_fields(newEvents, 'seq_pos', [0] * events.shape[0])
	initial_event_index = int(eventalign[0][idx['event_index']])
	initial_ref_index = int(eventalign[0][idx['ref_pos']])
	last_ref_index = initial_ref_index-1
	chromosome = eventalign[0][idx['contig']]
	start=False
	numSkips = 0
	numStays = 0
	
	# note that eventalign must be run with --scale-events of shift/drift will be an issue
	for line in eventalign:
		
		if not line[idx['strand']] == 't':
			# finished reading the template, now reading complement of the same read - skip
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
		elif last_ref_index < current_ref_index-1:
			# does skipping two bases at once count as 2 skips or one???
			numSkips += 1
		last_ref_index = current_ref_index
		last_seq_pos = seq_pos
		
		newEvents["seq_pos"][int(line[idx['event_index']])] = seq_pos
		newEvents["kmer"][int(line[idx['event_index']])] = seq
	
	fast5.close()
	#try:
	trainType = writeFast5(options, newEvents, fast5Path, initial_ref_index, last_ref_index, chromosome, forward, numSkips, numStays, kmer)
	return [fast5File, trainType]
	#except Exception as e:
#		log("Failed to write " + fast5File, 0, options)
#		log(str(e), 0, options)
#		return ["",""]

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
						log("Failed to write {}.npy: {}".format(filename, e), 0, options)
				if options.numReads != -1 and n >= options.numReads:
					# we're done!
					break
				
				# new read, initialise
				tmp = []
				current_read_name = line[idx['read_name']]
				try:
					fast5Path = refs[current_read_name]
				except KeyError:
					skip = True
					log(fast5Path + " not found", 0, options)
					continue
				fast5Name = fast5Path.split('/')[-1]
				if line[idx['strand']] != 't':
					# not template, skip
					skip=True
					continue
				
				outfile = getOutfile(fast5Name, options)
				filename = os.path.join(options.tempDir, fast5Name)
				if (not options.force) and os.path.isfile(outfile):
					log(outfile + " already exists. Use --force to recompute.", 0, options)
					skip=True
					premadeFilenames.append(fast5Name)
					continue
				elif (not options.force) and os.path.isfile(filename + ".npy"):
					log("{}.npy already exists. Use --force to recompute.".format(filename), 2, options)
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

def processReadWrapper(args):
   return processRead(*args)

def appendPremade(options, filenames, trainData):
	for fn in filenames:
		if np.random.rand() < options.valFraction:
			trainType="val"
		else:
			trainType="train"
		trainData.append([trainType, fn])
	return trainData

def writeTrainfiles(options, trainData):
	trainFilename = options.outPrefix + ".train.txt"
	valFilename = options.outPrefix + ".val.txt"
	
	with open(trainFilename, 'w') as trainFile, open(valFilename, 'w') as valFile:

		# write trainFile / valFile headers
		trainFile.write("#filename\n")
		valFile.write("#filename\n")
		
		for line in trainData:
			if line[1] == 'train':
				trainFile.write(line[0] + "\n")
			elif line[1] == 'val':
				valFile.write(line[0] + "\n")

def embedEventalign(options, fasta, eventalign):
	makeDir(os.path.join(options.outPrefix, __outdir_name__))
	log("Loading fast5 names and references...", 1, options)
	refs = loadRef(fasta)
	options.genome = loadGenome(options.genome)
	log("Splitting eventalign into separate files...", 1, options)
	filenames, idx, premadeFilenames = writeTempFiles(options, eventalign, refs)
	pool = Pool(options.threads)
	log("Embedding labels into fast5 files...", 1, options)
	trainData = pool.map(processReadWrapper, [[options, idx, i] for i in filenames])
	trainData = appendPremade(options, premadeFilenames, trainData)
	writeTrainfiles(options, trainData)
