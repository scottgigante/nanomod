from utils import callSubProcess
from . import __exe__
import os

def buildEventalign(options):
	
	fastaFile = '{}.fasta'.format(options.outPrefix)
	# poretools needs to be run from directory where fasta file will be output
	fastaDir = os.path.dirname(fastaFile)
	fastaBasename = os.path.basename(fastaFile)
	cwd = os.getcwd()
	try:
		os.chdir(fastaDir)
	except OSError:
		# already in cwd
		pass
	callSubProcess('{} fasta {}'.format(__exe__['poretools'], os.path.join(cwd, 
			options.reads)), options, newFile=fastaBasename, 
			outputFile=fastaBasename)
	os.chdir(cwd)
	
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