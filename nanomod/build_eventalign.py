################################################################################
#                                                                              #
# build_eventalign.py: a simple script to run the eventalign pipeline, from    #
# Metrichor labelled reads to the final eventalign tsv file. Pipeline as       #
# described in nanopolish eventalign's README.                                 #
#                                                                              #
# This file is part of Nanomod. Say something about a GNU licence?             #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

import os
import glob
from multiprocessing import Pool
from functools import partial
import subprocess

from utils import callSubProcess, log
from . import __exe__

# custom use of subprocess for poretools to allow special return value for multiprocessing
#
# @args call The system call to be run
# @args options Namespace object from argparse
# @args idx The index of fasta file to be written
# @return Name of the fasta file that was written
def callPoretools(call, options, idx):
	outfile = os.path.join(options.tempDir, "{}.fasta".format(idx))
	log("Poretools: {} files, output to {}.".format(len(call)-2, outfile), 2, 
			options)
	f = open(outfile, 'w')
	subprocess.call(call, stdout=f, close_fds=True, shell=False)
	return outfile

# wrapper script to multiprocess callPoretools
#
# @args args Array version of args to callPoretools
# @return Name of the fasta file that was written
def callPoretoolsWrapper(args):
	return callPoretools(*args)

# runs poretools multithreaded on small chunks of fasta files
# 
# Could be deprecated by find . -name "*.fast5" | parallel "poretools fasta" 
#
# TODO: why can't this handle more than eight cores? Arbitrary.
#
# @args poretools Path to poretools
# @args options Namespace object from argparse
# @args output Path to final output fasta file
# @return None
def multithreadPoretools(poretools, options, output):
	# check first that we actually want to edit
	if callSubProcess("touch {}".format(output), options, newFile=output) == 1:
		return 1
		
	cwd = os.getcwd()
	os.chdir(options.reads)
	files = [os.path.join(cwd, options.reads, f) for f in glob.glob("*.fast5")]
	os.chdir(cwd)
	prefix = [poretools, "fasta"]
	
	filesPerCall = 1000 # in practice, os probably supports more
	calls = [prefix + files[i:i + filesPerCall] for i in xrange(0, len(files),
			filesPerCall)]
	args = [[call, options] for call in calls]
	idx=0
	for arg in args:
		arg.append(idx)
		idx += 1
	
	pool = Pool(min(options.threads, 8)) # can't handle more than eight cores
	tempFastaList = pool.map(callPoretoolsWrapper, args)
	
	# write results to a single fasta
	with open(output, 'w') as outfile:
		for fa in tempFastaList:
			with open(fa) as infile:
				for line in infile:
					outfile.write(line)
	
# main script. Build eventalign file according to nanopolish's pipeline.
#
# @args options Namespace object from argparse
# @return fastaFile Filename of fasta corresponding to reads
# @return eventalignFile Filename of eventalign tsv
def buildEventalign(options):
	
	fastaFile = '{}.fasta'.format(options.outPrefix)
	multithreadPoretools(__exe__['poretools'], options, fastaFile)
	#poretoolsMaxFiles = 1000
	#callSubProcess(('find {} -name "*.fast5" | parallel -l {} "{} ' 
	#		'fasta"').format(options.reads, poretoolsMaxFiles, 
	#		__exe__['poretools']), options, newFile=fastaFile, 
	#		outputFile=fastaFile)
	
	callSubProcess('{} index {}'.format(__exe__['bwa'], options.genome), 
			options, newFile="{}.fai".format(options.genome))
	
	sortedBamFile = "{}.sorted.bam".format(options.outPrefix)
	callSubProcess(('{} mem -x ont2d -t {} {} {} | samtools view -Sb - '
			'| samtools sort -f - {}').format(__exe__['bwa'], options.threads,
			options.genome, fastaFile, sortedBamFile), options, 
			newFile=sortedBamFile)
	
	callSubProcess('{} index {}'.format(__exe__['samtools'], sortedBamFile), 
			options, newFile="{}.bai".format(sortedBamFile))
	
	eventalignFile = "{}.eventalign".format(options.outPrefix)
	callSubProcess(('{} eventalign -t {} --print-read-names -r {} -b {} -g {}' 
			' --models {}').format(__exe__['nanopolish'], options.threads, 
			fastaFile, sortedBamFile, options.genome, options.nanopolishModels), 
			options, newFile=eventalignFile, outputFile=eventalignFile)
	
	return fastaFile, eventalignFile