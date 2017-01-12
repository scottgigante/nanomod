################################################################################
#                                                                              #
# __main__.py: Runs Scott Gigante's Nanomod.                                   #
# Input: Reference genome and fast5 reads                                      #
# Output: fasta, bam, eventalign, labelled fast5 reads and (eventually) a      #
# trained network for nanonetcall                                              #
#                                                                              #
# This file is part of Nanomod. Say something about a GNU licence?             #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

# python -m nanomod -r data/Nott_R9_run2/downloads/pass/ -g data/ecoli_k12.fasta -o data/Nott_R9_run2 -v --temp-dir temp; python -m nanomod -r data/2016_09_11_E_COLI_MTHYLTD_R9/reads/downloads/pass/ -g data/ecoli_k12.fasta -o data/2016_09_11_E_COLI_MTHYLTD_R9 -v --temp-dir temp --methyl --force --data-fraction 0.1; python -m nanomod -r data/2016_09_18_E_COLI_NON_MTHYLTD_R9/reads/downloads/pass/ -g data/ecoli_k12.fasta -o data/2016_09_18_E_COLI_NON_MTHYLTD_R9 -v --temp-dir temp --force --data-fraction 0.1 

import argparse
import subprocess
import os
import sys
from multiprocessing import cpu_count
import tempfile
import shutil

from utils import log, makeDir, callSubProcess
from build_eventalign import buildEventalign
from embed_eventalign import embedEventalign
from . import init

# create directories, move things to the right place, etc
#
# @args options Namespace object from argparse
# @return None
def initialiseArgs(options):
	# check executables
	init()
	
	# make the temp dir
	makeDir(options.tempDir)
	options.tempDir = os.path.join(options.tempDir, 
			os.path.basename(options.outPrefix))
	makeDir(options.tempDir)
	makeDir(options.outPrefix)
	
	# move models files to current working directory for eventalign
	modelsFile = os.path.basename(options.nanopolishModels)
	cwd = os.getcwd()
	newModels = os.path.join(cwd, modelsFile)
	if (not options.force and newModels != options.nanopolishModels and 
			os.path.isfile(newModels)):
		log(('{} already exists in current working directory. Use --force to '
				'overwrite.').format(modelsFile), 0, options)
		options.nanopolishModels = newModels
	else:
		log('Copying {} to current working directory.'.format(
				options.nanopolishModels), 1, options)
		modelsPath = os.path.dirname(options.nanopolishModels)
		with open(options.nanopolishModels, 'r') as models:
			for line in models.readlines():
				model = line.strip()
				shutil.copy(os.path.join(modelsPath, model), 
						os.path.join(cwd, model))
		shutil.copy(options.nanopolishModels, newModels)

# move temporary files out of current working directory, delete temp directory
#
# @args options Namespace object from argparse
# @return None
def clean(options):
	# delete temp files
	shutil.rmtree(options.tempDir)
	
	# delete models from cwd
	modelsFile = os.path.basename(options.nanopolishModels)
	cwd = os.getcwd()
	newModels = os.path.join(cwd, modelsFile)
	if newModels != options.nanopolishModels:
		with open(options.nanopolishModels, 'r') as models:
			for line in models.readlines():
				model = line.strip()
				os.remove(os.path.join(cwd, model))
		os.remove(newModels)

# parse command line args
#
# @args argv sys.argv
# @return options Namespace object from argparse
def parseArgs(argv):
	
	#command line options
	parser = argparse.ArgumentParser(prog="nanomod")
	parser.add_argument("-r","--reads",dest="reads", required=True, 
			help="Directory in which fast5 reads are stored")
	parser.add_argument("-g", "--genome", required=True, 
			dest="genome", help="reference genome in fasta format")
	parser.add_argument("-o", "--output-prefix", dest="outPrefix", 
			required=True, help="prefix for nanomod output files")
	parser.add_argument("-k", "--kmer", type=int, default=3,
			help="Length of kmer for network training.")
	parser.add_argument("-t", "--threads", type=int, default=cpu_count(), 
			dest="threads", 
			help="number of threads to be used in multiprocessing")
	parser.add_argument("-p", "--parallel-sequences", default=125, 
			type=int, dest="parallelSequences", 
			help="Number of parallel threads of execution in GPU training.")
	parser.add_argument("-v","--verbose", default=0, action="count", 
			dest="verbosity", help=("can be used multiple times (e.g., -vv) "
			"for greater levels of verbosity."))
	parser.add_argument("--temp-dir", 
			default="{}/nanomod".format(tempfile.gettempdir()), 
			dest="tempDir", help="directory for nanomod temporary files")
	parser.add_argument("--nanopolish-models", 
			default="models/nanopolish_models.fofn", dest="nanopolishModels", 
			help=("nanopolish models for eventalign." 
			" Note: will be copied to current working directory"))
	parser.add_argument("--currennt-template", dest="currenntTemplate",
			default="models/default_template.json", 
			help="Template file for CURRENNT network architecture.")
	parser.add_argument("--force", default=False, action="store_true", 
			dest="force", help="force recreation of extant files")
	parser.add_argument("--num-reads", type=int, default=-1, dest="numReads", 
			help="limit the number of reads to be computed")
	parser.add_argument("--val-fraction", type=float, default=0.05, 
			dest="valFraction", help="fraction of data used for validation set")
	parser.add_argument("--data-fraction", type=float, default=1,
			dest="dataFraction", help="fraction of data to be sent to nanonet")
	parser.add_argument("--methyl", default=False, action="store_true", 
			dest="methyl")
	
	#parse command line options
	options = parser.parse_args()
	
	return options

# run main script
# 
# @args argv sys.argv
# @return None
def run(argv):
	options = parseArgs(argv)
	initialiseArgs(options)
	try:
		fasta, eventalign = buildEventalign(options)
		embedEventalign(options, fasta, eventalign)
		#callSubProcess(("./scripts/select_data_fraction.sh {0} {1}.train.txt"
		#		" {1}.val.txt").format(options.dataFraction, options.outPrefix), 
		#		options)
		#trainNanonet(options)
	finally:
		# we'd better clean up after ourselves, even if it crashes
		clean(options)

if __name__ == "__main__":
	run(sys.argv)
