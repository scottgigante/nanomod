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

from utils import callSubProcess, multiprocessWrapper
from seq_tools import loadGenome, getSeqDiff

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

def processBase(pileupColumn, modPositions, genome, modDir):
	refPos = pileupColumn.reference_pos
	
	# check if there are any modifiable bases in this position on reference
	checkFwd = False
	checkReverse = False
	refName = pileupColumn.reference_name	
	if refPos in modPositions[refName]['+']:
		checkFwd = True
	if (-refPos-1+len(genome)) in modPositions[refName]['-']:
		checkReverse = True
	if not checkFwd or checkReverse:
		# ref genome is not modified here
		return
	
	# count modifications
	numMods = 0
	numReads = 0
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
			continue
		if not (read.alignment.query_sequence[readPos] == genome[refPos]):
			# mismatch
			continue
		
		if checkReadModification(read.alignment.query_name, readPos, 
				read.alignment.query_length, read.alignment.cigartuples, 
				reverse, modDir) is not False:
			numMods += 1
		numReads += 1
	return [refPos, numMods, numReads] # TODO: are these returned in order? If so, remove reference pos

def processWindow(bamFile, contig, startPos, endPos, modPositions, genome, 
		modDir):
	bam = pysam.AlignmentFile(bamFile, "rb")
	output = np.ndarray(shape=(0,3))
	for base in bam.pileup(contig, startPos, endPos):
		refPos = base.reference_pos
		if refPos >= startPos and refPos < endPos:
			# window is expanded for some reason
			baseOutput = processBase(base, modPositions, genome, modDir)
			if baseOutput is not None:
				output = np.append(output, [baseOutput], axis=0)
	bam.close()
	return output
	
def processWindowWrapper(args):
	return multiprocessWrapper(processWindow, args)

def countModifications(bamFile, modDir, options):
	# TODO: allow specifying a window to reduce time cost
	modPositions, genome = checkModifiedPositions(options)
	modCounts = dict()
	
	basesPerCall = 2000 
	# pileup seems to take a minimum of 2000 even if you specify less
	for contig in genome:
		refLength = len(genome[contig]['seq'])
		windowArray = [i for i in xrange(1, refLength, basesPerCall)]
		p = Pool(options.threads)
		argsList = [[bamFile, contig, i, i+basesPerCall, modPositions, 
				genome[contig]['seq'], modDir] for i in windowArray]
		# TODO: need to run pileup within single thread. Should be okay to read from multiple locations.
		windows = p.map(processWindowWrapper, argsList)
		if len(windows) > 1:
			modCounts[contig] = np.array(windows[0])
			for i in range(1,len(windows)):
				modCounts[contig] = np.append(modCounts[contig], windows[i],
						axis=0)
			modCounts[contig] = modCounts[contig].tolist() # TODO: should we avoid using numpy? Keep numpy and remove JSON (probably.)
	
	return modCounts
	

	
