################################################################################
#                                                                              #
# This file is part of Nanomod.                                                #
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
# utils.py: helper scripts for Nanomod.                                        #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

from __future__ import print_function
import os
import sys
import subprocess

__log_levels__ = ['[warning] ','[debug] ','[log] ']

# print a message to log depending on verbosity
#
# @args message The message to be logged
# @args level The level of urgency of the message:
#		0: warning
#		1: debug
# 		2: log
# @args options Namespace object from argparse
def log(message, level, options):
	if level > options.verbosity:
		return
	if level == 0:
		output = sys.stderr
	else:
		output = sys.stdout
	print(__log_levels__[level] + message, file=output)

# call subprocess to create a new file using an external script if the file 
# doesn't already exist
#
# TODO: can we move away from shell calls?
#
# @args call The system call to be run
# @args options Namespace object from argparse
# @args newFile File to be created
# @args outputFile File to which we shall pipe stdout
# @args close_fds Boolean to close outputFile after writing
# @args shell Use local shell
# @args mode Writing mode: choose from w, wb, a
def callSubProcess(call, options, newFile=None, outputFile=None, close_fds=True, 
		shell=True, mode='w'):
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
			return 1
	if outputFile is not None:
		out=open(outputFile, mode)
	elif options.verbosity < 2:
		# don't print stdout from subprocess
		out=open(os.devnull, 'w')
	else:
		out=sys.stdout
	subprocess.call(call, stdout=out, close_fds=close_fds, shell=shell)
	
	# check the file was created, if given
	assert(newFile is None or os.path.isfile(newFile))
	return 0

# create a directory, if it doesn't already exist
#
# @args dir The directory to be created
# @return None	
def makeDir(dir):
	try:
		os.mkdir(dir)
	except OSError:
		# already exists
		pass
		
# multiprocessing.Pool.map() wrapper for functions taking more than one argument
# @args args an array of arguments, the first of which is the name of the function to be called
# @return return value of called function
def multiprocessWrapper(args):
	func = args[0]
	args = args[1:]
	try:
		return func(*args)
	except Exception as e:
		print('Caught exception in worker thread: {}({})'.format(func, 
				", ".join(args)))
		traceback.print_exc()
		print()
		raise e