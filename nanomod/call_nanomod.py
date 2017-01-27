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
# call_nanomod.py: Call base modifications using a pre-trained neural network. #
# Input: fast5 reads and nanomod trained network                               #
# Output: fasta, bam, and file indexing probability of methylation per-base.   #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 19 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

import argparse
import json
from multiprocessing import cpu_count

from utils import callSubProcess
from index_modifications import indexAndCleanModifications
from summarise_modifications import countModifications
from seq_tools import loadGenome
from . import __exe__

def parseArgs(argv):
	"""parse command line args
	@args argv sys.argv
	@return options Namespace object from argparse
	"""
	
	#command line options
	parser = argparse.ArgumentParser(prog="nanomod",
			description=("Call DNA base modifications, using a pre-trained "
			"nanomod network. Requires: nanonetcall, bwa."), 
			epilog=("Example usage: nanomod call --reads " 
			"sample_data/r9/modified --genome sampla_data/ecoli_k12.fasta "
			"--model models/5mc.nanomod.npy --output-prefix data/test --threads"
			" 16"))
	parser.add_argument("-r","--reads", dest="reads", 
			required=True, 
			help=("Directory in which fast5 reads are stored (required)"))
	parser.add_argument("-g", "--genome", required=True, 
			dest="genome", help="Reference genome in fasta format (required)")
	parser.add_argument("-o", "--output-prefix", dest="outPrefix", 
			required=True, help="Prefix for nanomod output files")
	parser.add_argument("-m", "--model", 
			default="models/5mc.nanomod.npy", #required=True,
			help="Nanomod pre-trained network (required)")
	parser.add_argument("-s", "--sequence-motif", dest="sequenceMotif",
			nargs=2, default=["CG","MG"], 
			help=("Motif of canonical and modified site (e.g., -s CG MG) "
			"(required)")) # TODO: include in model?
	parser.add_argument("-t", "--threads", type=int, default=cpu_count(), 
			dest="threads", 
			help="Number of threads to be used in multiprocessing")
	parser.add_argument("-v","--verbose", default=0, action="count", 
			dest="verbosity", help=("Can be used multiple times (e.g., -vv) "
			"for greater levels of verbosity."))
	parser.add_argument("--chemistry", default="r9") # TODO: can we infer this?
	parser.add_argument("--num-reads", type=int, default=-1, dest="numReads", 
			help="Limit the number of reads to be analysed")
	parser.add_argument("--force", default=False, action="store_true", 
			dest="force", help="Force recreation of extant files")
	
	#parse command line options
	options = parser.parse_args(argv)
	
	return options

# run main script
# 
# @args argv sys.argv
# @return None
def callNanomod(argv):
	options = parseArgs(argv)
	
	fastaFile = "{}.fasta".format(options.outPrefix)
	callSubProcess([__exe__['nanonetcall'], "--chemistry", options.chemistry, 
			"--jobs", options.threads, "--model", options.model, "--output",
			fastaFile, "--section", "template", options.reads], options, 
			shell=False, newFile=fastaFile)
	
	unmodifiedFastaFile, modDir = indexAndCleanModifications(fastaFile, options)
	
	sortedBamFile = "{}.sorted.bam".format(options.outPrefix)
	# TODO: we do this twice - separate into a new file?
	callSubProcess(('{} mem -x ont2d -t {} {} {} | samtools view -Sb - '
			'| samtools sort -f - {}').format(__exe__['bwa'], options.threads,
			options.genome, unmodifiedFastaFile, sortedBamFile), options, 
			outputFile=sortedBamFile, newFile=sortedBamFile)
	callSubProcess([__exe__['samtools'], 'index', sortedBamFile], options, 
			shell=False)
	
	modificationCounts = countModifications(sortedBamFile, modDir, options)
	
	outFile = "{}.methylation.txt".format(options.outPrefix)
	with open(outFile, 'w') as handle:
		json.dump(modificationCounts, fp=handle, indent=4)
