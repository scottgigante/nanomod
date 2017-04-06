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
# TODO: include numReads as an alternative to dataFraction
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 11 Jan 2017                                                            #
#                                                                              #
################################################################################

import h5py
import sys
import os
import fnmatch
import numpy as np
from multiprocessing import Pool
import json
import logging
from operator import itemgetter

from . import __modes__

def getStats(ary):
    return np.mean(ary), np.std(ary)

def checkProbs(filename):
    try:
        with h5py.File(filename, 'r') as fh:
            attrs = fh.get("Analyses/Basecall_1D_000/Summary/basecall_1d_template").attrs
            numEvents = float(attrs["called_events"])
            skipProb = float(attrs["num_skips"])/numEvents
            stayProb = float(attrs["num_stays"])/numEvents
            stepProb = 1-skipProb-stayProb
            qscore = float(attrs["mean_qscore"])
            readLength = int(attrs["sequence_length"])
    except IOError:
        logging.warning("Failed to open {}".format(filename)) # TODO: use logging here
        skipProb, stayProb, stepProb, qscore = 1, 1, 0, 0, 0
    return [skipProb, stayProb, stepProb, qscore, readLength]

# get the empirical skip and stay probabilities for a set of reads
# @param dir The directory to search for fast5 files
# @param proportion The desired proportion of reads to be retained
# @param mode Choose any or all of 'skip', 'stay', and 'step'
# @return maxSkips
# @return maxStays
# @return minSteps
def selectBestReads(dir, proportion, mode=__modes__, threads=1, readLength=2000):
    for m in mode:
        if m not in __modes__:
            logging.warning("Mode {} not recognised")
    files = []
    for root, _, filenames in os.walk(dir):
        for filename in fnmatch.filter(filenames, '*.fast5'):
            files.append(os.path.join(root, filename))
    selectNum = max(int(len(files) * proportion),1)

    if "random" in mode:
        return files[np.random.choice(len(files), selectNum, replace=False)]

    p = Pool(threads)
    probs = np.array(p.map(checkProbs, files), dtype=np.float32).transpose()

    skipMean, skipStdv = getStats(probs[0])
    stayMean, stayStdv = getStats(probs[1])
    stepMean, stepStdv = getStats(probs[2])
    qscoreMean, qscoreStdv = getStats(probs[3])

    cutoffs = []
    for i in xrange(len(probs.transpose())):
        read = probs.transpose()[i]
        filename = files[i]
        if read[4] >= readLength:
            deviations = []
            if "skip" in mode:
                deviations.append((read[0] - skipMean)/skipStdv)
            if "stay" in mode:
                deviations.append((read[1] - stayMean)/stayStdv)
            if "step" in mode:
                deviations.append((stepMean - read[2])/stepStdv)
            if "qscore" in mode:
                deviations.append((qscoreMean - read[3])/qscoreStdv)
            cutoffs.append([filename, np.mean(deviations)])
    cutoffs = sorted(cutoffs, key=itemgetter(1), reverse=True)

    bestFiles = list(map(itemgetter(0), cutoffs))

    if not len(bestFiles) >= selectNum:
        logging.warning("Insufficient reads of adequate read length. Reducing data fraction to {0:.2f} ({} reads).".format(float(len(bestFiles)) / len(files), selectNum))
    else:
        logging.info("Selecting {} reads.".format(selectNum))
        bestFiles = bestFiles[:selectNum]
    return bestFiles
