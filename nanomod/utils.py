from __future__ import print_function
import os
import sys
import subprocess

__log_levels__ = ['[warning] ','[debug] ','[log] ']

def log(message, level, options):
	if level > options.verbosity:
		return
	if level == 0:
		output = sys.stderr
	else:
		output = sys.stdout
	print(__log_levels__[level] + message, file=output)

def callSubProcess(call, options, newFile=None, outputFile=None):
	# we should use shell=False to make this more portable
	# print call to debug
	log(call, 1, options)
	
	# check if we are accidentally overwriting something
	if newFile is not None and os.path.isfile(newFile):
		if options.force:
			log("{} already exists. Forcing recreation.".format(newFile), 
					1, options)
		else:
			log("{} already exists. Use --force to recompute.".format(newFile), 
					0, options)
			return
	if outputFile is not None:
		out=open(outputFile, 'w')
	elif options.verbosity < 2:
		# don't print stdout from subprocess
		out=open(os.devnull, 'w')
	else:
		out=sys.stdout
	subprocess.call(call, stdout=out, close_fds=True, shell=True)
	
	
	# check the file was created, if given
	assert(newFile is None or os.path.isfile(newFile))
