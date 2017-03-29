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

from train_nanomod import trainNanomod
from call_nanomod import callNanomod
from parse_args import parseArgs
from . import init

def run():
    """
    Run Nanomod.
    """
    command, options = parseArgs()
    # check executables
    init(command)
    if command == "train":
        trainNanomod(options)
    elif command == "call":
        callNanomod(options)
    else:
        # panic
        raise argparse.ArgumentError(command, "Command {} not found".format(command))

if __name__ == "__main__":
    """
    Run from command line.
    """
    run()
