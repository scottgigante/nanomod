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
# check_skip_stay_prob.py:                                                     #
#                                                                              #
# TODO: add method to check bounds given desired proportion                    #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 11 Jan 2017                                                            #
#                                                                              #
################################################################################

import h5py
import sys
import os
import numpy as np
from multiprocessing import Pool
import json

def getStats(ary):
	return np.mean(ary), np.std(ary)

def checkProbs(filename):
	with h5py.File(filename, 'r') as fh:
		attrs = fh.get("Analyses/Basecall_1D_000/Summary/basecall_1d_template").attrs
		numEvents = float(attrs["called_events"])
		skipProb = attrs["num_skips"]/numEvents
		stayProb = attrs["num_stays"]/numEvents
		stepProb = 1-skipProb-stayProb
	return [skipProb, stayProb, stepProb]

# get the empirical skip and stay probabilities for a set of reads
# @param dir The directory to search for fast5 files
# @param proportion The desired proportion of reads to be retained
# @return maxSkips
# @return maxStays
# @return minSteps
def getSkipStayConstraints(dir, proportion):
	files = [os.path.join(dir,file) for file in os.listdir(dir) if file.endswith(".fast5")]
	
	p = Pool()
	probs = np.array(p.map(checkProbs, files)).transpose()	
	
	skipMean, skipStdv = getStats(probs[0])
	stayMean, stayStdv = getStats(probs[1])
	stepMean, stepStdv = getStats(probs[2])
	
	cutoffs = []
	for read in probs.transpose():
		cutoffs.append(max((read[0] - skipMean)/skipStdv, (read[1] - stayMean)/stayStdv, (stepMean - read[2])/stepStdv))
	cutoffs.sort()
	numStdvs = cutoffs[int(len(cutoffs)*proportion)]
	
	maxSkips = skipMean + numStdvs * skipStdv
	maxStays = stayMean + numStdvs * stayStdv
	minSteps = stepMean - numStdvs * stepStdv
	
	return {'maxSkips' : maxSkips,
			'maxStays' : maxStays, 
			'minSteps' : minSteps }

if __name__ == "__main__":
	print getSkipStayConstraints(sys.argv[1], sys.argv[2])