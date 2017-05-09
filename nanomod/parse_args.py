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
import os
from multiprocessing import cpu_count
import tempfile
import shutil
import logging
import numpy as np
import re

from utils import makeDir, configureLog
from pickle_to_currennt import convertPickle
from expand_model_alphabet import expandModelAlphabet
from . import __modes__

__description__ = {
    'train' : """Generates a neural network model to call DNA base
                modifications, using Oxford Nanopore Technologies' Nanonettrain.
                Requires: poretools, bwa, samtools, nanopolish, nanonettrain""",
    'call' : """Call DNA base modifications, using a pre-trained nanomod network.
                Requires: nanonetcall, bwa."""
    }
__epilog__ = {
    'train' : """Example usage: nanomod --canonical-reads sample_data/r9/canonical
                --modified-reads sample_data/r9/modified --genome
                sample_data/ecoli_k12_mg1655.fa --output-prefix data/test
                --sequence-motif CG MG --threads 16 --verbose --data-fraction 0.5
                --select-mode random""",
    'call' : """Example usage: nanomod call --reads
                sample_data/r9/modified --genome sampla_data/ecoli_k12.fasta
                --model models/5mc.nanomod.npy --output-prefix data/test --threads
                16"""
    }

def parseRegion(region):
    """
    Parse a region specification.

    :param region: String region specification

    :raises argparse.ArgumentTypeError: raises an error if format not recognised

    :returns contig: String contig / chromosome name
    :returns start: integer start position (0-based)
    :returns end: integer end position (1-based)

    >>> parseRegion("chr1:1000-2000")
    ("chr1", 1000, 2000)

    """
    region = ''.join(region.split()) # remove whitespace
    region = re.split(':|-', region)
    start = 0
    end = None
    if len(region) < 1:
        raise argparse.ArgumentTypeError("Region must specify a reference name")
    elif len(region) > 3:
        raise argparse.ArgumentTypeError("Region format not recognized")
    else:
        contig = region[0]
        try:
            start = int(re.sub(",|\.", "", region[1]))
        except IndexError:
            pass
        try:
            end = int(re.sub(",|\.", "", region[2]))
        except IndexError:
            pass
        finally:
            return contig, start, end

def initialiseTrainArgs(options):
    """
    Perform various initalisation tasks for training.

    :param options: Namespace object from argparse

    :returns: Namespace object from argparse
    """
    # make the temp dir
    makeDir(options.tempDir)
    options.tempDir = os.path.join(options.tempDir,
            os.path.basename(options.outPrefix))
    makeDir(options.tempDir)
    makeDir(options.outPrefix)

    # restrict dataFraction to [0,1]
    options.dataFraction = min(max(options.dataFraction, 0),1)

    # make sure sequenceMotif is uppercase
    for i in range(len(options.sequenceMotif)):
        options.sequenceMotif[i] = options.sequenceMotif[i].upper()

    # convert ONT model to json
    options.currenntTemplate = "{}.json".format(options.nanonetTemplate)
    convertPickle(options)

    # expand model to include modifications
    options.expandedTemplate = "{}.mod.json".format(options.nanonetTemplate)
    expandModelAlphabet(options.currenntTemplate, options.expandedTemplate, options.kmer,
            options.sequenceMotif, options.force)

    # verify select modes
    if "random" in options.selectMode and options.selectMode != ["random"]:
        logging.warning("--select-mode random cannot be used in conjunction with other modes. Using random only.")

def initialiseArgs(command, options):
    """
    Perform various initalisation tasks.

    :param command: String The nanomod command to be executed
    :param options: Namespace object from argparse

    :returns: Namespace object from argparse
    """
    # sequence motif should all be uppercase
    options.sequenceMotif = [s.upper() for s in options.sequenceMotif]

    np.random.seed(options.seed)

    options.threads = max(min(options.threads, cpu_count()), 1)

    if command == "train":
        initialiseTrainArgs(options)
    elif command == "call":
        # currently no call specific initialisation
        pass
    else:
        # this shouldn't ever happen
        raise NameError("{} command not defined".format(command))

    return options

def addCommonArgs(parser):
    """
    Add arguments common to all running modes

    :param parser: The argument parser to be modified

    :returns: The modified argument parser
    """
    parser.add_argument("-g", "--genome", required=True,
            dest="genome", help="Reference genome in fasta format (required)")
    parser.add_argument("-o", "--output-prefix", dest="outPrefix",
            required=True, help="Prefix for nanomod output files")
    parser.add_argument("-s", "--sequence-motif", dest="sequenceMotif",
            nargs=2, default=["CG","MG"], metavar="CANONICAL MODIFIED",
            help=("Motif of canonical and modified site (e.g., -s CG MG) "
            "(required)")) # TODO: include in model?
    parser.add_argument("-v","--verbose", default=0, action="count",
            dest="verbosity", help=("Can be used multiple times (e.g., -vv) "
            "for greater levels of verbosity."))
    parser.add_argument("--num-reads", type=int, default=-1, dest="numReads",
            help="Limit the number of reads to be analysed")
    parser.add_argument("--force", default=False, action="store_true",
            dest="force", help="Force recreation of extant files") # TODO: make this ranked? eg force 1, force 8?
    parser.add_argument("--no-normalize", default=False, action="store_true",
            dest="noNormalize",
            help="do not apply median normalization before run") # TODO: can we infer this?
    parser.add_argument("-t", "--threads", type=int, default=cpu_count(),
            dest="threads",
            help="Number of threads to be used in multiprocessing")
    parser.add_argument("--seed", type=int, default=None,
            help="Random seed for read selection.")
    return parser

def addCallArgs(parser):
    """
    Add arguments unique to nanomod call

    :param parser: The argument parser to be modified

    :returns: The modified argument parser
    """
    #command line options
    parser.add_argument("-r","--reads", dest="reads",
            required=True,
            help=("Directory in which fast5 reads are stored (required)"))
    parser.add_argument("-m", "--model",
            default="models/5mc.nanomod.npy", #required=True,
            help="Nanomod pre-trained network (required)")
    parser.add_argument("--region", type=parseRegion, metavar="CHR:START-END",
            help="Name of contig to analyze")
    parser.add_argument("--window", type=int, metavar="SIZE",
            help="Size of window over which to aggregate calls")
    parser.add_argument("--chemistry", default="r9") # TODO: can we infer this?
    parser.add_argument("--alpha", default=0.5, type=float,
            help="Parameter for uninformative Beta prior")

    return parser

def addTrainArgs(parser):
    """
    Add arguments unique to nanomod train

    :param parser: The argument parser to be modified

    :returns: The modified argument parser
    """
    parser.add_argument("-c","--canonical-reads", dest="canonicalReads",
            required=True,
            help=("Directory in which canonical-base fast5 reads are stored "
            "(required)"))
    parser.add_argument("-m", "--modified-reads", dest="modifiedReads",
            help="Directory in which modified-base fast5 reads are stored")
    parser.add_argument("-k", "--kmer", type=int, default=5,
            help="Length of kmer for network training")
    parser.add_argument("-p", "--parallel-sequences", default=80,
            type=int, dest="parallelSequences",
            help="Number of parallel threads of execution in GPU training")
    parser.add_argument("--temp-dir",
            default="{}/nanomod".format(tempfile.gettempdir()),
            dest="tempDir", help="Directory for nanomod temporary files")
    parser.add_argument("--nanonet-template", dest="nanonetTemplate",
            default="models/r9_template.npy",
            help="Nanonet model file for network initialisation")
    parser.add_argument("--error-rate", default=0.001, type=float,
            dest="rate", help="Randomly introduce errors to avoid overfitting.")
    parser.add_argument("--val-fraction", type=float, default=0.05,
            dest="valFraction", help="Fraction of data used for validation set")
    parser.add_argument("--data-fraction", type=float, default=0.1,
            dest="dataFraction", help="Fraction of data to be sent to nanonet (overwritten by --num-reads)")
    parser.add_argument("--select-mode", choices=tuple(__modes__), default=["skip", "stay", "qscore"], nargs="*",
            help=("Method for choosing reads to send to nanonet; choose any or"
            " all from {}").format(", ".join(__modes__)),
            dest="selectMode")
    parser.add_argument("--read-length", default=5000, type=int,
            dest="readLength", help="Minimum read length for training reads.")
    parser.add_argument("--mapping-quality", default=60, type=int,
            dest="mappingQuality", help="Minimum mapping quality for training reads.")

    return parser

def parseCommandArgs(command, argv):
    """
    Common actions for parsing command line arguments
    :param command: Nanomod command to be called
    :param argv: The remainder of the command line arguments after command is removed
    :raises NameError: raises error if a nonexistent command is called
    :returns: Namespace object from argparse
    """
    parser = argparse.ArgumentParser(prog="nanomod",
            description=(__description__[command]),
            epilog=__epilog__[command])
    parser = addCommonArgs(parser)
    if command == "call":
        parser = addCallArgs(parser)
    elif command == "train":
        parser = addTrainArgs(parser)
    else:
        raise NameError("{} command not defined".format(command))

    options = parser.parse_args(argv)
    configureLog(options.verbosity)
    logging.debug("nanomod {} {}".format(command, " ".join(argv)))

    initialiseArgs(command, options)

    return options

def parseArgs():
    """
    Parse Nanomod command line arguments

    :raises argparse.ArgumentError: raises error if a nonexistent command is called

    :returns command: Nanomod command to be called
    :returns options: argparse namespace for command line arguments
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
    command = parser.add_argument("command", choices=("train","call"))
    parser.add_argument("options", nargs=argparse.REMAINDER)

    #parse command line options
    initialOptions = parser.parse_args()

    try:
        options = parseCommandArgs(initialOptions.command, initialOptions.options)
    except NameError as e:
        if str(e).endswith(" command not defined"):
            raise argparse.ArgumentError(command, "Command {} not found".format(initialOptions.command))
        else:
            raise

    return initialOptions.command, options
