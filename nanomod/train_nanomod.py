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
# train_nanomod.py: Train a CURRENNT network for calling base modifications.   #
# Input: Reference genome and fast5 reads                                      #
# Output: fasta, bam, eventalign, labelled fast5 reads and (eventually) a      #
# trained network for nanonetcall                                              #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

# python -m nanomod -c data/Nott_R9_run2/downloads/pass/ -g data/ecoli_k12.fasta -o data/Nott_R9_run2/nanomod --select-mode skip step stay --data-fraction 0.1 -v --temp-dir temp --force

import os
import sys
import shutil
import logging
import numpy as np

from utils import callSubProcess, recursiveFindFast5, preventOverwrite
from build_eventalign import buildEventalign
from embed_eventalign import embedEventalign
from train_nanonet import trainNanonet

def clean(options):
    """
    Clean up after ourselves.
    Move temporary files out of current working directory, delete temp directory.

    :param options: Namespace object from argparse
    """
    # delete temp files
    shutil.rmtree(options.tempDir)

def buildTrainingSet(options, reads, outPrefix, modified=False):
    """
    Build training set for a single MinION run

    :param options: Namespace from argparse
    :param reads: Directory for fast5 reads
    :param outPrefix: Prefix for output files
    """
    fasta, eventalign = buildEventalign(options, reads, outPrefix)
    embedEventalign(options, fasta, eventalign, reads, outPrefix, modified)

def combineTrainingSets(options, t1, t2, outPrefix):
    """
    Combine training sets from multiple nanomod runs.
    Primary function is to combine canonical and modified training sets

    :param options: Namespace from argparse
    :param t1: Prefix for training set 1
    :param t2: Prefix for training set 2
    :param outPrefix: Prefix for output files
    """
    trainFiles = ("{}.train.txt".format(t1), "{}.train.txt".format(t2), "{}.train.txt".format(outPrefix))
    valFiles = ("{}.val.txt".format(t1), "{}.val.txt".format(t2), "{}.val.txt".format(outPrefix))
    for file1, file2, combinedFile in [trainFiles, valFiles]:
        # grab header
        callSubProcess("head -n 1 {}".format(file1), options.force,
                stdout=combinedFile, mode='w')
        # remove headers
        callSubProcess("sed -i '1d' {}".format(file1), options.force)
        callSubProcess("sed -i '1d' {}".format(file2), options.force)
        callSubProcess("paste -d '\n' {} {}".format(file1, file2), options.force,
                stdout=combinedFile, mode='a')

def equaliseTrainingSets(f1, f2):
    """
    Equalise size of training sets by selecting randomly from the larger list.

    :param f1: list of strings Paths to files from training set 1
    :param f2: list of strings Paths to files from training set 2

    :returns: f1, f2 Tuple of reduced lists
    """

    l1, l2 = len(f1), len(f2)
    if l1 > l2:
        f1 = np.array(f1)[np.random.choice(l1, l2, replace=False)]
    elif l2 > l1:
        f2 = np.array(f2)[np.random.choice(l2, l1, replace=False)]
    logging.info("Randomly restricting both datasets to {} files.".format(max(l1, l2)))
    return f1, f2

def trainNanomod(options):
    """
    Run nanomod training

    :param options: Namespace from argparse
    """
    try:
        if options.modifiedReads is not None:
            canonicalPrefix = "{}.canonical".format(options.outPrefix)
            modifiedPrefix = "{}.modified".format(options.outPrefix)
            canonicalFiles = recursiveFindFast5(options.canonicalReads)
            modifiedFiles = recursiveFindFast5(options.modifiedReads)
            canonicalFiles, modifiedFiles = equaliseTrainingSets(canonicalFiles,
                    modifiedFiles)
            buildTrainingSet(options, canonicalFiles, canonicalPrefix)
            buildTrainingSet(options, modifiedFiles, modifiedPrefix,
                    modified=True)
            combineTrainingSets(options, canonicalPrefix, modifiedPrefix,
                    options.outPrefix)
            trainNanonet(options)
        else:
            canonicalFiles = recursiveFindFast5(options.canonicalReads)
            buildTrainingSet(options, canonicalFiles, options.outPrefix)
            trainNanonet(options)
    finally:
        # we'd better clean up after ourselves, even if it crashes
        clean(options)
