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

from utils import callSubProcess
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

    # delete models from cwd
    modelsFile = os.path.basename(options.nanopolishModels)
    cwd = os.getcwd()
    newModels = os.path.join(cwd, modelsFile)
    if newModels != options.nanopolishModels:
        with open(options.nanopolishModels, 'r') as models:
            for line in models.readlines():
                model = line.strip()
                os.remove(os.path.join(cwd, model))
        os.remove(newModels)

def buildTrainingSet(options, reads, outPrefix, modified=False):
    """
    Build training set for a single MinION run

    :param options: Namespace from argparse
    :param reads: Directory for fast5 reads
    :param outPrefix: Prefix for output files
    """
    fasta, eventalign = buildEventalign(options, reads, outPrefix)
    embedEventalign(options, fasta, eventalign, reads, outPrefix, modified)
    if "random" in options.selectMode:
        # overwrite small training sets with random selection
        callSubProcess(("./nanomod/scripts/select_data_fraction.sh {0} "
                "{1}.train.txt {1}.val.txt").format(options.dataFraction,
                outPrefix), options.force)

def combineTrainingSets(options, t1, t2, outPrefix):
    """
    Combine training sets from multiple nanomod runs.
    Primary function is to combine canonical and modified training sets

    :param options: Namespace from argparse
    :param t1: Prefix for training set 1
    :param t2: Prefix for training set 2
    :param outPrefix: Prefix for output files
    """
    trainFile1 = "{}.train.txt.small".format(t1)
    trainFile2 = "{}.train.txt.small".format(t2)
    valFile1 = "{}.val.txt.small".format(t1)
    valFile2 = "{}.val.txt.small".format(t2)
    trainFileCombined = "{}.train.txt.small".format(outPrefix)
    valFileCombined = "{}.val.txt.small".format(outPrefix)
    callSubProcess("cat {}".format(trainFile1), options.force,
            outputFile=trainFileCombined, mode='w')
    # skip header second time
    callSubProcess("tail -n +2 {}".format(trainFile2), options.force,
            outputFile=trainFileCombined, mode='a')
    callSubProcess("cat {}".format(valFile1), options.force,
            outputFile=valFileCombined, mode='w')
    # skip header second time
    callSubProcess("tail -n +2 {}".format(valFile2), options.force,
            outputFile=valFileCombined, mode='a')

def trainNanomod(options):
    """
    Run nanomod training

    :param options: Namespace from argparse
    """
    try:
        if options.modifiedReads is not None:
            canonicalPrefix = "{}.canonical".format(options.outPrefix)
            modifiedPrefix = "{}.modified".format(options.outPrefix)
            # TODO: should we modify datafraction in case we have drastically different amounts of data for modified and canonical?
            buildTrainingSet(options, options.canonicalReads, canonicalPrefix)
            buildTrainingSet(options, options.modifiedReads, modifiedPrefix,
                    modified=True)
            combineTrainingSets(options, canonicalPrefix, modifiedPrefix,
                    options.outPrefix)
            trainNanonet(options)
        else:
            buildTrainingSet(options, options.canonicalReads, options.outPrefix)
            trainNanonet(options)
    finally:
        # we'd better clean up after ourselves, even if it crashes
        clean(options)
