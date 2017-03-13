################################################################################
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
# __main__.py: Runs Scott Gigante's Nanomod.                                   #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 19 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

import argparse
import sys

from train_nanomod import trainNanomod
from call_nanomod import callNanomod
from . import init

def parseArgs():
    """parse command line args
    @args argv sys.argv
    @return options Namespace object from argparse
    """
    
    #command line options
    parser = argparse.ArgumentParser(prog="nanomod",
            description=("Nanopore base modification caller."), 
            epilog=("Commands:\n"
            "train\tTrain a neural network for a new type of modification\n"
            "call\tUse an existing trained network to call modifications on a"
            " set of fast5 files\n\n"
            "Example usage: nanomod call -r sample_data/r9/modified -m "
            "models/5mc.nanomod -o data/test.fa"),
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("command", choices=("train","call"))
    parser.add_argument("options", nargs=argparse.REMAINDER)
    
    #parse command line options
    options = parser.parse_args()
    
    return options

# run main script
# 
# @args argv sys.argv
# @return None
def run():
    options = parseArgs()
    # check executables
    init(options.command)
    if options.command == "train":
        trainNanomod(options.options)
    elif options.command == "call":
        callNanomod(options.options)
    else:
        # panic
        raise argparse.ArgumentError(command, "Command {} not found".format(options.command))
        

if __name__ == "__main__":
    run()
