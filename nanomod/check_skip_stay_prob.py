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
# TODO: rename file to something more meaningful
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
import logging
from operator import itemgetter

from . import __modes__

def getStats(ary):
    """
    Get the mean and standard deviation of an array

    :param ary: array-like Array to be checked

    :returns mean: float Mean of ary
    :returns stdv: float Standard deviation of ary
    """
    return np.mean(ary), np.std(ary)

def checkProbs(filename):
    """
    Check probabilities and related stats for a fast5 file

    :param filename: string Path to fast5 file

    :returns: List of stats: skip probability, stay probability, step probability, mean quality score, read length

    TODO: rename this to something more meaningful
    """
    try:
        with h5py.File(filename, 'r') as fh:
            attrs = fh.get("Analyses/Basecall_1D_000/Summary/basecall_1d_template").attrs
            numEvents = float(attrs["called_events"])
            skipProb = float(attrs["skip_prob"])
            stayProb = float(attrs["stay_prob"])
            stepProb = 1-skipProb-stayProb
            qscore = float(attrs["mean_qscore"])
            readLength = int(attrs["sequence_length"])
    except IOError:
        logging.warning("Failed to open {}".format(filename))
        skipProb, stayProb, stepProb, qscore, readLength = 1, 1, 0, 0, 0
    except AttributeError:
        logging.warning("{}: summary not found".format(filename))
        skipProb, stayProb, stepProb, qscore, readLength = 1, 1, 0, 0, 0
    return [skipProb, stayProb, stepProb, qscore, readLength]

def selectBestReads(files, proportion, numReads, mode=__modes__, threads=1, readLength=5000):
    """
    Choose the best subset of reads from a directory

    :param files: list of strings Paths to files to be analysed
    :param proportion: float The desired proportion of reads to be retained (overwritten by numReads)
    :param numReads: int The number of reads to be retained
    :param mode: list of strings Any or all of 'skip', 'stay', 'step', 'qscore', and 'random'
    :param threads: int Number of threads to run in parallel
    :param readLength: int Minimum read length to retain

    :returns bestFiles: list of strings Paths to files to be retained
    """
    # check modes
    for m in mode:
        if m not in __modes__:
            logging.warning("Mode {} not recognised")

    # check number of reads to retain
    if numReads > 0:
        selectNum = min(len(files), numReads)
    else:
        selectNum = max(int(len(files) * proportion),1)

    if "random" in mode:
        # not much to do
        files = np.array(files)
        logging.info("Selecting {} reads.".format(selectNum))
        return files[np.random.choice(len(files), selectNum, replace=False)]

    # get read stats
    p = Pool(threads)
    probs = np.array(p.map(checkProbs, files), dtype=np.float32).transpose()

    # find stats on each read
    skipMean, skipStdv = getStats(probs[0])
    stayMean, stayStdv = getStats(probs[1])
    stepMean, stepStdv = getStats(probs[2])
    qscoreMean, qscoreStdv = getStats(probs[3])

    # sort files by stats
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
    cutoffs = sorted(cutoffs, key=itemgetter(1))

    bestFiles = list(map(itemgetter(0), cutoffs))

    # select best reads
    if not len(bestFiles) >= selectNum:
        logging.warning("Insufficient reads of adequate read length. Reducing data fraction to {0:.2f} ({} reads).".format(float(len(bestFiles)) / len(files), selectNum))
    else:
        logging.info("Selecting {} reads.".format(selectNum))
        bestFiles = bestFiles[:selectNum]
    return bestFiles
