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
# call_nanomod.py: Call base modifications using a pre-trained neural network. #
# Input: fast5 reads and nanomod trained network                               #
# Output: fasta, bam, and file indexing probability of methylation per-base.   #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 19 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

import numpy as np
import h5py
import os
import multiprocessing
import functools
import logging

from utils import callSubProcess, configureLog, preventOverwrite, makeDir, recursiveFindFast5
from index_modifications import indexAndCleanModifications
from summarise_modifications import countModifications
from shutil import copyfile
from build_eventalign import buildSortedBam
from . import __exe__

def normalizeRead(outPrefix, input):
    read, subfolder = input
    outDir = os.path.join(outPrefix, "nanomod", str(subfolder))
    makeDir(outDir)
    newRead = os.path.join(outDir, os.path.basename(read))
    copyfile(read, newRead)
    try:
        with h5py.File(newRead, 'r+') as outFile:
            raw_group = outFile["Raw/Reads"].values()[0]
            signal = raw_group["Signal"][()]
            med = np.median(signal)
            median_deviation = signal - med
            MAD = np.median(abs(median_deviation))
            signal = median_deviation / MAD
            try:
                del outFile["Analyses"]
            except KeyError:
                # no analyses, no worries
                pass
            del raw_group["Signal"]
            raw_group["Signal"] = signal
    except Exception as e:
        logging.warning("Failed to normalize file {}: {}".format(os.path.basename(read), str(e)))
        os.remove(newRead)

def normalizeReads(outPrefix, reads, threads):
    """
    Normalise raw signal using median normalization

    :param reads: String path to fast5 files

    :returns: String path to normalized fast5 files

    TODO: does nanonetcall use raw signal or event data?
    """
    # TODO: we do this twice - make this into a routine.
    newReads = os.path.join(outPrefix, "nanomod")
    makeDir(outPrefix)
    makeDir(newReads)
    readList = recursiveFindFast5(reads)
    argList = ((readList[i], i//4000) for i in range(len(readList)))
    pool = multiprocessing.Pool(threads)
    pool.map(functools.partial(normalizeRead, outPrefix), argList)
    pool.close()
    pool.join()
    return newReads

def callNanomod(options):
    """
    Call base modifications on reads based on a previously built nanomod model

    :param options: Namespace object from argparse
    """

    fastaFile = "{}.fasta".format(options.outPrefix)
    if not preventOverwrite(fastaFile, options.force):
        if not options.noNormalize:
            # median normalise raw signal
            options.reads = normalizeReads(options.outPrefix, options.reads, options.threads)
        fastaFile = "{}.fasta".format(options.outPrefix)
        args = [__exe__['nanonetcall'], "--chemistry", options.chemistry,
                "--jobs", str(options.threads), "--model", options.model,
                "--output", fastaFile, "--section", "template", options.reads]
        if options.numReads > 0:
            args.extend(["--limit", str(options.numReads)])
        callSubProcess(args, options.force, shell=False, newFile=fastaFile)

    unmodifiedFastaFile, modDir = indexAndCleanModifications(fastaFile, options)

    sortedBamFile = buildSortedBam(options.threads, options.genome, unmodifiedFastaFile, options.outPrefix, options.force)

    fmt, header, modCounts = countModifications(sortedBamFile, modDir, options)
    #wig = getWiggleTrack(modificationCounts)

    outFile = "{}.methylation.txt".format(options.outPrefix)
    with open(outFile, 'wb') as handle:
        np.savetxt(handle, modCounts, fmt=fmt, header=header)
# TODO: specify dtype
