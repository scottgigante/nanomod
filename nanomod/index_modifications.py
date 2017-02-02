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

import os
import json
from multiprocessing import Pool

from seq_tools import unmodifyFasta, __canonical__
from utils import makeDir, multiprocessWrapper, preventOverwrite, log

def indexRead(record, outputDir, options):
	outputFile = os.path.join(outputDir, record.id)
	if preventOverwrite(outputFile, options):
		return
	modIndex = dict()
	for i in range(len(record.seq)):
		if record.seq[i] not in __canonical__: 
			try:
				modIndex[record.seq[i]].append(i)
			except KeyError:
				# not yet initialised
				modIndex[record.seq[i]] = [i]
	with open(outputFile, 'wb') as handle:
		json.dump(modIndex, fp=handle)

def indexReadWrapper(args):
	return multiprocessWrapper(indexRead, args)

def indexAndCleanModifications(fastaFile, options):
	outputDir = "{}.modIndex".format(options.outPrefix)
	makeDir(outputDir)
	
	log("Cleaning fasta of modifications...", 1, options)
	# TODO: should we avoid overwrite here? would have to reload fasta
	unmodifiedFastaFile = "{}.unmodified.fasta".format(options.outPrefix)
	fasta = unmodifyFasta(fastaFile, unmodifiedFastaFile, options.sequenceMotif)
	
	log("Indexing modifications on reads...", 1, options)
	p = Pool(options.threads)
	p.map(indexReadWrapper, [[i, outputDir, options] for i in fasta])
	
	return unmodifiedFastaFile, outputDir
	

	
