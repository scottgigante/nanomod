import os
import glob
from multiprocessing import Pool
from functools import partial
import subprocess

from utils import callSubProcess, log
from . import __exe__

def callPoretools(call, options, idx):
	outfile = os.path.join(options.tempDir, "{}.fasta".format(idx))
	log("Poretools: {} files, output to {}.".format(len(call)-2, outfile), 2, 
			options)
	f = open(outfile, 'w')
	subprocess.call(call, stdout=f, close_fds=True, shell=False)
	return outfile

def callPoretoolsWrapper(args):
	return callPoretools(*args)

def multithreadPoretools(poretools, options, output):
	# check first that we actually want to edit
	if callSubProcess("touch {}".format(output), options, newFile=output) == 1:
		return 1
		
	cwd = os.getcwd()
	os.chdir(options.reads)
	files = [os.path.join(cwd, options.reads, f) for f in glob.glob("*.fast5")]
	os.chdir(cwd)
	prefix = [poretools, "fasta"]
	
	filesPerCall = 1000 # in practice, os probably supports ~200
	calls = [prefix + files[i:i + filesPerCall] for i in xrange(0, len(files),
			filesPerCall)]
	args = [[call, options] for call in calls]
	idx=0
	for arg in args:
		arg.append(idx)
		idx += 1
	
	pool = Pool(min(options.threads, 8))
	tempFastaList = pool.map(callPoretoolsWrapper, args)
	
	# write results to a single fasta
	with open(output, 'w') as outfile:
		for fa in tempFastaList:
			with open(fa) as infile:
				for line in infile:
					outfile.write(line)
	

def buildEventalign(options):
	
	fastaFile = '{}.fasta'.format(options.outPrefix)
	multithreadPoretools(__exe__['poretools'], options, fastaFile)
	
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