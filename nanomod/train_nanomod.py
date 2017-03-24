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

import argparse
import subprocess
import os
import sys
from multiprocessing import cpu_count
import tempfile
import shutil
import logging

from utils import makeDir, callSubProcess, configureLog
from build_eventalign import buildEventalign
from embed_eventalign import embedEventalign
from pickle_to_currennt import convertPickle
from expand_model_alphabet import expandModelAlphabet
from . import __modes__

# create directories, move things to the right place, etc
#
# @args options Namespace object from argparse
# @return None
def initialiseArgs(options):
    
    # make the temp dir
    makeDir(options.tempDir)
    options.tempDir = os.path.join(options.tempDir, 
            os.path.basename(options.outPrefix))
    makeDir(options.tempDir)
    makeDir(options.outPrefix)
    
    # make sure sequenceMotif is uppercase
    for i in range(len(options.sequenceMotif)):
        options.sequenceMotif[i] = options.sequenceMotif[i].upper()
    
    # move models files to current working directory for eventalign
    modelsFile = os.path.basename(options.nanopolishModels)
    cwd = os.getcwd()
    newModels = os.path.join(cwd, modelsFile)
    if (not options.force and newModels != options.nanopolishModels and 
            os.path.isfile(newModels)):
        logging.warning(('{} already exists in current working directory. Use --force to '
                'overwrite.').format(modelsFile))
        options.nanopolishModels = newModels
    else:
        logging.info('Copying {} to current working directory.'.format(
                options.nanopolishModels))
        modelsPath = os.path.dirname(options.nanopolishModels)
        with open(options.nanopolishModels, 'r') as models:
            for line in models.readlines():
                model = line.strip()
                shutil.copy(os.path.join(modelsPath, model), 
                        os.path.join(cwd, model))
        shutil.copy(options.nanopolishModels, newModels)
    
    # convert ONT model to json
    options.currenntTemplate = "{}.json".format(options.nanonetTemplate)
    convertPickle(options)
    
    # expand model to include modifications
    options.expandedTemplate = "{}.mod.json".format(options.nanonetTemplate)
    expandModelAlphabet(options)

# move temporary files out of current working directory, delete temp directory
#
# @args options Namespace object from argparse
# @return None
def clean(options):
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

# 
def parseArgs(argv):
    """parse command line args
    @args argv sys.argv
    @return options Namespace object from argparse
    """
    
    #command line options
    parser = argparse.ArgumentParser(prog="nanomod",
            description=("Generates a neural network model to call DNA base "
            "modifications, using Oxford Nanopore Technologies' Nanonettrain."
            " Requires: poretools, bwa, samtools, nanopolish, nanonettrain"), 
            epilog=("Example usage: nanomod --canonical-reads " 
            "sample_data/r9/canonical --modified-reads sample_data/r9/modified "
            "--genome sample_data/ecoli_k12_mg1655.fa --output-prefix data/test"
            " --sequence-motif CG MG --threads 16 --verbose --data-fraction 0.5"
            " --select-mode random"))
    parser.add_argument("-c","--canonical-reads", dest="canonicalReads", 
            required=True, 
            help=("Directory in which canonical-base fast5 reads are stored " 
            "(required)"))
    parser.add_argument("-m", "--modified-reads", dest="modifiedReads",
            help="Directory in which modified-base fast5 reads are stored")
    parser.add_argument("-g", "--genome", required=True, 
            dest="genome", help="Reference genome in fasta format (required)")
    parser.add_argument("-o", "--output-prefix", dest="outPrefix", 
            required=True, help="Prefix for nanomod output files")
    parser.add_argument("-s", "--sequence-motif", dest="sequenceMotif",
            nargs=2, default=["CG","MG"], 
            help=("Motif of canonical and modified site (e.g., -s CG MG) "
            "(required if modified reads are provided)"))
    parser.add_argument("-k", "--kmer", type=int, default=5,
            help="Length of kmer for network training")
    parser.add_argument("-t", "--threads", type=int, default=cpu_count(), 
            dest="threads", 
            help="Number of threads to be used in multiprocessing")
    parser.add_argument("-p", "--parallel-sequences", default=125, 
            type=int, dest="parallelSequences", 
            help="Number of parallel threads of execution in GPU training")
    parser.add_argument("-v","--verbose", default=0, action="count", 
            dest="verbosity", help=("Can be used multiple times (e.g., -vv) "
            "for greater levels of verbosity."))
    parser.add_argument("--temp-dir", 
            default="{}/nanomod".format(tempfile.gettempdir()), 
            dest="tempDir", help="Directory for nanomod temporary files")
    parser.add_argument("--nanopolish-models", 
            default="models/nanopolish_models.fofn", dest="nanopolishModels", 
            help=("Nanopolish models for eventalign. " 
            "Note: will be copied to current working directory"))
    parser.add_argument("--nanonet-template", dest="nanonetTemplate",
            default="models/r9_template.npy", 
            help="Nanonet model file for network initialisation")
    parser.add_argument("--error-rate", default=0.001, type=float, 
    		dest="rate", help="Randomly introduce errors to avoid overfitting.")
    parser.add_argument("--force", default=False, action="store_true", 
            dest="force", help="Force recreation of extant files")
            # TODO: make this ranked? eg force 1, force 8?
    parser.add_argument("--val-fraction", type=float, default=0.05, 
            dest="valFraction", help="Fraction of data used for validation set")
    parser.add_argument("--data-fraction", type=float, default=1,
            dest="dataFraction", help="Fraction of data to be sent to nanonet")
    parser.add_argument("--select-mode", default=__modes__, nargs="*",
            help=("Method for choosing reads to send to nanonet; choose any or" 
            " all from random, {}").format(", ".join(__modes__)), 
            dest="selectMode")
    parser.add_argument("--read-length", default=2000, type=int, 
            dest="readLength", help="Minimum read length for training reads.")
    parser.add_argument("--no-normalize", default=False, action="store_true", 
            dest="noNormalize", 
            help="do not apply median normalization before run")
    parser.add_argument("--num-reads", type=int, default=-1, dest="numReads", 
            help="Limit the number of reads to be analysed")
    
    #parse command line options
    options = parser.parse_args(argv)
    configureLog(options.verbosity)
    
    return options

# build training set for a single minION run
#
# @args options Namespace from argparse
# @args reads Directory for fast5 reads
# @args outPrefix Prefix for output files
def buildTrainingSet(options, reads, outPrefix, modified=False):
    fasta, eventalign = buildEventalign(options, reads, outPrefix)
    embedEventalign(options, fasta, eventalign, reads, outPrefix, modified)
    if "random" in options.selectMode:
        # overwrite small training sets with random selection
        callSubProcess(("./nanomod/scripts/select_data_fraction.sh {0} "
                "{1}.train.txt {1}.val.txt").format(options.dataFraction,
                outPrefix), options)

def combineTrainingSets(options, t1, t2, outPrefix):
    trainFile1 = "{}.train.txt.small".format(t1)
    trainFile2 = "{}.train.txt.small".format(t2)
    valFile1 = "{}.val.txt.small".format(t1)
    valFile2 = "{}.val.txt.small".format(t2)
    trainFileCombined = "{}.train.txt.small".format(outPrefix)
    valFileCombined = "{}.val.txt.small".format(outPrefix)
    callSubProcess("cat {}".format(trainFile1), options, 
            outputFile=trainFileCombined, mode='w')
    # skip header second time
    callSubProcess("tail -n +2 {}".format(trainFile2), options, 
            outputFile=trainFileCombined, mode='a')
    callSubProcess("cat {}".format(valFile1), options, 
            outputFile=valFileCombined, mode='w')
    # skip header second time
    callSubProcess("tail -n +2 {}".format(valFile2), options, 
            outputFile=valFileCombined, mode='a')

# run main script
# 
# @args argv sys.argv
# @return None
def trainNanomod(argv):
    options = parseArgs(argv)
    initialiseArgs(options)
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
            #trainNanonet(options)
        else:
            buildTrainingSet(options, options.canonicalReads, options.outPrefix)
    finally:
        # we'd better clean up after ourselves, even if it crashes
        clean(options)