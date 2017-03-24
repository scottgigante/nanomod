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
import traceback
import json
import logging

__log_levels__ = [logging.WARNING, logging.INFO, logging.DEBUG]

# set up logging
def configureLog(level=0):
    if level >= len(__log_levels__):
        level = len(__log_levels__) - 1
    elif level < 0:
        level = 0
    logging.basicConfig(format='[%(levelname)s] %(message)s', level=__log_levels__[level])

def preventOverwrite(file, options):
    """ Check if we are accidentally overwriting something
    @param file Path to file to be written
    @param options Namespace from argparse """
    if file is not None and os.path.isfile(file):
        if options.force:
            logging.info("{} already exists. Forcing recreation.".format(file))
        else:
            logging.warning("{} already exists. Use --force to recompute.".format(file))
            return True
    return False

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
    try:
        # we should use shell=False to make this more portable
        # print call to debug
        logging.info(call if shell else " ".join(call))
    
        if preventOverwrite(newFile, options):
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
    except OSError as e:
        logging.warning("OSError: has a dependency changed path? Try deleting nanomod/.*.config.json and try again.")
        # TODO: fix this automatically - maybe a --check-dependencies option?
        raise e
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

def loadJson(filename):
    with open(filename, 'r') as inFh:
        data = json.load(inFh)
    return data

def saveJson(filename, data):
    with open(filename, 'w') as outFh:
        json.dump(data, fp=outFh, indent=4)
        
# multiprocessing.Pool.map() wrapper for functions taking more than one argument
# @args func the function to be called
# @args args an array of arguments for the function
# @return return value of called function
def multiprocessWrapper(func, args):
    try:
        return func(*args)
    except Exception as e:
        logging.error('Caught exception in worker thread:')
        traceback.print_exc()
        print()
        raise e