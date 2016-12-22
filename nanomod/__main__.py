#!/usr/bin/env python

# Runs Scott Gigante's Nanomod.

import argparse
import subprocess
import os
import sys
from multiprocessing import cpu_count
import tempfile
import shutil

from utils import log
from build_eventalign import buildEventalign

# create directories, move things to the right place, etc
def initialiseArgs(options):
	# make the temp dir
	try:
		os.mkdir(options.tempDir)
	except OSError:
		# already exists
		pass
	
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

def parseArgs(argv):
	
	#command line options
	parser = argparse.ArgumentParser(prog="nanomod")
	parser.add_argument("-r","--reads",dest="reads",required=True, 
			help="Directory in which fast5 reads are stored")
	parser.add_argument("-g", "--genome", default="data/ecoli_k12.fasta", 
			dest="genome", help="reference genome in fasta format") # ??? required
	parser.add_argument("-o", "--output-prefix", default="data/nanomod", 	
			dest="outPrefix", help="prefix for nanomod output files")
	parser.add_argument("-t", "--threads", type=int, default=cpu_count(), 
			dest="threads", 
			help="number of threads to be used in multiprocessing")
	parser.add_argument("--temp-dir", default="/tmp/embed_eventalign", 
			dest="tempDir", help="directory for nanomod temporary files")
	parser.add_argument("--nanopolish-models", 
			default="models/nanopolish_models.fofn", dest="nanopolishModels", 
			help=("nanopolish models for eventalign." 
			" Note: will be copied to current working directory"))
	parser.add_argument("--cache_path", default=tempfile.gettempdir(), 
			help="Path for temporary files.")
	parser.add_argument("-v","--verbose", default=0, action="count", 
			dest="verbosity", help=("can be used multiple times (e.g., -vv) "
			"for greater levels of verbosity."))
	parser.add_argument("--force", default=False, action="store_true", 
			dest="force", help="force recreation of extant files")
	
	#parse command line options
	options = parser.parse_args()
	
	return options

def run(argv):
	options = parseArgs(argv)
	initialiseArgs(options)
	buildEventalign(options)
	clean(options)

if __name__ == "__main__":
	run(sys.argv)
