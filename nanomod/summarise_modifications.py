################################################################################
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
# index_modifications.py: Index where modified bases are in each read          #
# Input: fast5 reads and nanomod trained network                               #
# Output: fasta, bam, and file indexing probability of methylation per-base.   #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 19 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

from multiprocessing import Pool, Lock
from Bio import SeqIO
import pysam
import os
import json
import numpy as np

from utils import callSubProcess, multiprocessWrapper, log
from seq_tools import loadGenome, getSeqDiff

def getBayesPercentage(numMods, numUnmods, numReads, alpha):
	# TODO: can we use coverage vs. numMods + numUnmods to get a clever value for alpha?
	return float(numMods + alpha)/(numMods + numUnmods + alpha)

def checkModifiedPositions(options):
	modPositions = dict()
	canonGenome = loadGenome(options)
	modGenome = loadGenome(options, modified=True)
	for contig in canonGenome:
		modPositions[contig] = dict()
		modPositions[contig]['+'] = getSeqDiff(
				canonGenome[contig]['seq'], 
				modGenome[contig]['seq'])
		modPositions[contig]['-'] = getSeqDiff(
				canonGenome[contig]['reverseSeq'], 
				modGenome[contig]['reverseSeq'])	
	return modPositions, canonGenome

def checkReadModification(queryName, readPos, readLength, cigarTuples, reverse, 
		modDir):
	# TODO: check that this checks the correct position
	if reverse:
		readPos = -readPos-1 
		readPos += readLength
	startHardClip = cigarTuples[0][1] if cigarTuples[0][0] == 5 else 0
	readPos += startHardClip
	
	modFile = os.path.join(modDir, queryName)
	with open(modFile, 'rb') as handle:
		modIndex = json.load(handle)
	for base in modIndex:
		if readPos in modIndex[base]:
			return base
	return False

def processBase(pileupColumn, alpha, modPositions, genome, contig, modDir, 
		dtype):
	refPos = pileupColumn.reference_pos
	
	# check if there are any modifiable bases in this position on reference
	checkFwd = False
	checkReverse = False
	refName = pileupColumn.reference_name	
	if refPos in modPositions[refName]['+']:
		checkFwd = True
	if (-refPos-1+len(genome)) in modPositions[refName]['-']:
		# TODO: we have 345000 sites when we should have double that - is reverse working?
		checkReverse = True
	if not (checkFwd or checkReverse):
		# ref genome is not modified here
		return
	
	# count modifications
	numMods = 0
	numUnmods = 0
	numNoMatch = 0
	for read in pileupColumn.pileups:
		# check direction
		reverse = read.alignment.is_reverse
		if ((reverse and not checkReverse) or
				(not reverse and not checkFwd)):
			# wrong direction strand, nothing to check
			continue
		
		# check if we have a match
		readPos = read.query_position
		if readPos is None:
			# deletion
			numNoMatch += 1
			continue
		if not (read.alignment.query_sequence[readPos] == genome[refPos]):
			# mismatch
			numNoMatch += 1
			continue
		
		if checkReadModification(read.alignment.query_name, readPos, 
				read.alignment.query_length, read.alignment.cigartuples, 
				reverse, modDir) is not False:
			numMods += 1
		else:
			numUnmods += 1
	numReads = numMods + numUnmods + numNoMatch
	modProb = getBayesPercentage(numMods, numUnmods, numReads, alpha)
	return np.array([(contig, refPos, modProb, numMods, numUnmods, 
			numReads)], dtype=dtype)

def processWindow(bamFile, contig, startPos, endPos, alpha, modPositions, 
		genome, modDir, dtype):
	bam = pysam.AlignmentFile(bamFile, "rb")
	output = np.empty(shape=(0,), dtype=dtype)
	for base in bam.pileup(contig, startPos, endPos):
		refPos = base.reference_pos
		if refPos >= startPos and refPos < endPos:
			# window is expanded for some reason
			baseOutput = processBase(base, alpha, modPositions, genome, contig, 
					modDir, dtype)
			if baseOutput is not None:
				output = np.append(output, baseOutput)
	bam.close()
	return output
	
def processWindowWrapper(args):
	return multiprocessWrapper(processWindow, args)

def countModifications(bamFile, modDir, options):
	# TODO: allow specifying a window to reduce time cost
	log("Indexing modification sites on reference genome...", 1, options)
	modPositions, genome = checkModifiedPositions(options)
	
	# S32 string, <u4 uint32, <f4 float32, u1 uint8 <u2 uint16
	dtype = np.dtype("S32,<u4,<f4,u1,u1,u1")
	aggDtype = np.dtype("S32,<u4,<u4,<f4,<u2,<u2,<u2")
	if options.window is None or options.window < 2:
		fmt = "\t".join(["%s","%i","%lf","%i","%i","%i"])
		header = "\t".join(["Chromosome", "Position", "Modified Percentage", 
				"Count Modified", "Count Unmodified", "Coverage"])
		modCounts = np.empty(shape=(0,), dtype=dtype)
	else:
		fmt = "\t".join(["%s","%i","%i","%lf","%i","%i","%i"])
		header = "\t".join(["Chromosome", "Start", "Stop",
				"Modified Percentage", "Count Modified", "Count Unmodified", 
				"Coverage"])
		modCounts = np.empty(shape=(0,), dtype=aggDtype)
	
	log("Analyzing read modifications...", 1, options)
	basesPerCall = 2000 
	# pileup seems to take a minimum of 2000 even if you specify less
	for contig in genome:
		contigModCounts = np.empty(shape=(0,), dtype=dtype)
		refStart = 1
		refEnd = len(genome[contig]['seq'])
		if options.region is not None:
			if contig != options.region[0]:
				continue
			else:
				refStart = max(refStart, options.region[1])
				if options.region[2] is not None:
					refEnd = min(refEnd, options.region[2])
		
		windowArray = [i for i in xrange(refStart, refEnd, basesPerCall)]
		p = Pool(options.threads)
		argsList = [[bamFile, contig, i, min(i+basesPerCall, refEnd), 
				options.alpha, modPositions, genome[contig]['seq'], 
				modDir, dtype] for i in windowArray]
		windows = p.map(processWindowWrapper, argsList)
		for w in windows:
			contigModCounts = np.append(contigModCounts, w)
		
		if options.window is None or options.window < 2:
			modCounts = np.append(modCounts, contigModCounts)
		else:
			# aggregate calls
			backStep = (options.window - 1) / 2
			forwardStep = (options.window - 1) / 2 + (options.window - 1) % 2
			for i in range(contigModCounts.shape[0]):
				startPos = contigModCounts['f1'][i] - backStep
				endPos = contigModCounts['f1'][i] + forwardStep + 1
				# coarse window
				countsWindow = contigModCounts[range(max(0,i-backStep), 
						min(contigModCounts.shape[0],i+forwardStep+1))]
				afterStart = countsWindow['f1'] >= startPos
				beforeEnd = countsWindow['f1'] < endPos
				# fine window
				countsWindow = countsWindow[np.logical_and(afterStart, 
						beforeEnd)]
				numMods = sum(countsWindow['f3'])
				numUnmods = sum(countsWindow['f4'])
				numReads = sum(countsWindow['f5'])
				modProb = getBayesPercentage(numMods, numUnmods, numReads, 
							options.alpha)
				modCounts = np.append(modCounts, np.array([(contig, 
						startPos, endPos, modProb, numMods, numUnmods, 
						numReads)], dtype=aggDtype))
					
	return fmt, header, modCounts
		
