from utils import callSubProcess
from . import __exe__

def buildEventalign(options):
	fastaFile = '{}.fasta'.format(options.outPrefix)
	callSubProcess('{} fasta {}'.format(__exe__['poretools'], options.reads), 
			options, newFile=fastaFile, outputFile=fastaFile)
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