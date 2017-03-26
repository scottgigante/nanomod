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

def configureLog(level=0):
    """
    Set up logging based on user input to verbosity

    :param level: Positive integer, capped at 2
    """
    if level >= len(__log_levels__):
        # too high, not available
        level = len(__log_levels__) - 1
    elif level < 0:
        # too low, not available
        level = 0
    logging.basicConfig(format='[%(levelname)s] %(message)s', level=__log_levels__[level])

def preventOverwrite(file, force):
    """
    Check if we are accidentally overwriting something

    :param file: Path to file to be written
    :param force: Boolean value, force overwrite or not

    :returns: Boolean value, prevent overwrite or not
    """
    if file is not None and os.path.isfile(file):
        if force:
            logging.info("{} already exists. Forcing recreation.".format(file))
        else:
            logging.warning("{} already exists. Use --force to recompute.".format(file))
            return True
    return False

def callSubProcess(call, force, newFile=None, outputFile=None, close_fds=True,
        shell=True, mode='w'):
    """
    Call subprocess to create a new file using an external script if the file
    doesn't already exist
    TODO: can we move away from shell calls?

    :param call: The system call to be run as a string or array_like
    :param: Boolean value, True if we want to overwrite extant files, False otherwise
    :param newFile: File that should be created after run
    :param outputFile: File to which we shall pipe stdout
    :param close_fds: Boolean value, close outputFile after writing or not
    :param shell: Use local shell call - true for string call, false for array_like call
    :param mode: Writing mode: choose from w, wb, a
    """
    try:
        # we should use shell=False to make this more portable
        # print call to debug
        logging.info(call if shell else " ".join(call))

        if preventOverwrite(newFile, force):
            return 1

        if outputFile is not None:
            out=open(outputFile, mode)
        elif logging.getLogger().isEnabledFor(logging.DEBUG):
            out=sys.stdout # TODO: can we send this to logging?
        else:
            # don't print stdout from subprocess
            out=open(os.devnull, 'w')
        subprocess.call(call, stdout=out, close_fds=close_fds, shell=shell)

        # check the file was created, if given
        assert(newFile is None or os.path.isfile(newFile))
    except OSError as e:
        logging.warning("OSError: has a dependency changed path? Try deleting nanomod/.*.config.json and try again.")
        # TODO: fix this automatically - maybe a --check-dependencies option?
        raise e
    return 0

def makeDir(dir):
    """
    Create a directory, if it doesn't already exist

    :param dir: The directory to be created
    """
    try:
        os.mkdir(dir)
    except OSError:
        # already exists
        pass

def loadJson(filename):
    """
    Load a JSON file

    :param filename: string File to be loaded

    :returns: dictionary or array_like JSON output
    """
    with open(filename, 'r') as inFh:
        data = json.load(inFh)
    return data

def saveJson(filename, data):
    """
    Save a JSON file

    :param filename: string File to be loaded
    :param data: dictionary or array_like JSON input
    """
    with open(filename, 'w') as outFh:
        json.dump(data, fp=outFh, indent=4)

def multiprocessWrapper(func, args):
    """
    multiprocessing.Pool.map() wrapper for functions taking more than one argument

    :param func: the function to be called
    :param args: an array of arguments for the function

    :returns: return value of called function
    """
    try:
        return func(*args)
    except Exception as e:
        logging.error('Caught exception in worker thread:')
        traceback.print_exc()
        print()
        raise e