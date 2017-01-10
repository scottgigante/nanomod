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
from . import __exe__

# main script. Run nanonettrain
#
# @args options Namespace object from argparse
# @return None
def trainNanonet(options):
	
	callSubProcess(("{0} --train {1} --train_list {1}.train.txt --val {1} "
				   "--val_list {1}.val.txt --output {1} --kmer_length {2}"
				   "--parallel_sequences {3} --workspace {4} --cache_path {4}"
				   "--model {5} --cuda").format(__exe__['nanonettrain'], 
				   options.outPrefix, options.kmer, options.parallelSequences, 
				   options.tempDir, options.currenntTemplate)
	# return ???