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
import re
import h5py
import os

from utils import callSubProcess, configureLog
from index_modifications import indexAndCleanModifications
from summarise_modifications import countModifications
from seq_tools import loadGenome
from shutil import copyfile
from . import __exe__

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

def normaliseReads(reads):
    """
    Normalise raw signal using median normalization
    
    :param reads: String path to fast5 files
    
    :returns: String path to normalized fast5 files
    
    TODO: does nanonetcall use raw signal or event data?
    """
    # TODO: we do this twice - make this into a routine.
    newReads = os.path.join(options.reads, "nanomod")
    os.mkdir(newReads)
    for f in os.listdir(options.reads):
        if f.endswith(".fast5"):
            path = os.path.join(options.reads, f)
            newPath = os.path.join(options.reads, "nanomod", f)
            copyfile(path, newPath)
            with h5py.File(newPath, 'r+') as f:
                raw_group = f["Raw/Reads"].values()[0]
                signal = raw_group["Signal"][()]
                med = np.median(signal)
                median_deviation = signal - med
                MAD = np.median(abs(median_deviation))
                signal = median_deviation / MAD
                del raw_group["Signal"]
                del f["Analyses"]
                raw_group["Signal"] = signal
    return newReads

def callNanomod(options):
    """
    Call base modifications on reads based on a previously built nanomod model
     
    :param options: Namespace object from argparse
    """
    
    if not options.noNormalise:
        # median normalise raw signal
        options.reads = normaliseReads(options.reads)
    
    fastaFile = "{}.fasta".format(options.outPrefix)
    args = [__exe__['nanonetcall'], "--chemistry", options.chemistry, 
            "--jobs", str(options.threads), "--model", options.model, 
            "--output", fastaFile, "--section", "template", options.reads]
    if args.numReads > 0:
        args.extend(["--limit", str(options.numReads)])
    callSubProcess(args, options.force, shell=False, newFile=fastaFile)
    
    unmodifiedFastaFile, modDir = indexAndCleanModifications(fastaFile, options)
    
    sortedBamFile = buildSortedBam(options.threads, options.genome, fastaFile, options.force)
    
    fmt, header, modCounts = countModifications(sortedBamFile, modDir, options)
    #wig = getWiggleTrack(modificationCounts)
    
    outFile = "{}.methylation.txt".format(options.outPrefix)
    with open(outFile, 'wb') as handle:
        np.savetxt(handle, modCounts, fmt=fmt, header=header)
# TODO: specify dtype
