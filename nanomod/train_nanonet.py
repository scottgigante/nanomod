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
# train_nanonet.py: run Nanonettrain.                                          #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

from utils import callSubProcess
import logging
from . import __exe__
from seq_tools import expandAlphabet

def trainNanonet(options):
    """
    Train nanonet from nanomod output.

    At present, we just print the command for the user to run, since we assume nanomod is not run on a GPU-enabled machine.

    :param options: Namespace object from argparse

    TODO: what should we return? can we try running nanonet?
    """
    logging.info(("{0} --train {1} --train_list {1}.train.txt --val {1}"
                   " --val_list {1}.val.txt --output {1} --kmer_length "
                   "{2} --parallel_sequences {3} --workspace {4} --cache_path "
                   "{4} --model {5} --bases {6} --cuda").format(__exe__['nanonettrain'],
                   options.outPrefix, options.kmer, options.parallelSequences,
                   options.tempDir, options.expandedTemplate, "".join(expandAlphabet(options.sequenceMotif))))
    return None
