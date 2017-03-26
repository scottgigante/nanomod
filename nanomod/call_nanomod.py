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
            end = int(re.sub(",|\.", "", region[2]))
        except IndexError:
            # no worries, we don't have to specify a window
            pass
        finally:
            return contig, start, end

# run main script
# 
# @args argv sys.argv
# @return None
def callNanomod(options):
    
    if not options.noNormalise:
        # median normalise
        # TODO: we do this twice - make this into a routine.
        # TODO: oh no now we're modifying the raw data!
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
                    raw_group["Signal"] = signal
        options.reads = newReads
                
    
    fastaFile = "{}.fasta".format(options.outPrefix)
    callSubProcess([__exe__['nanonetcall'], "--chemistry", options.chemistry, 
            "--jobs", str(options.threads), "--model", options.model, 
            "--output",    fastaFile, "--limit", str(options.numReads), "--section", "template", options.reads], 
            options, shell=False, newFile=fastaFile)
    # TODO: don't inclide --limit if numReads is -1
    # TODO: fix bases
    
    unmodifiedFastaFile, modDir = indexAndCleanModifications(fastaFile, options)
    
    sortedBamFile = "{}.sorted.bam".format(options.outPrefix)
    # TODO: we do this twice - separate into a new file?
    callSubProcess(('{} mem -x ont2d -t {} {} {} | samtools view -Sb - '
            '| samtools sort -f - {}').format(__exe__['bwa'], options.threads,
            options.genome, unmodifiedFastaFile, sortedBamFile), options, 
            outputFile=sortedBamFile, newFile=sortedBamFile)    
    callSubProcess([__exe__['samtools'], 'index', sortedBamFile], options, 
            shell=False)
    
    fmt, header, modCounts = countModifications(sortedBamFile, modDir, options)
    #wig = getWiggleTrack(modificationCounts)
    
    outFile = "{}.methylation.txt".format(options.outPrefix)
    with open(outFile, 'wb') as handle:
        np.savetxt(handle, modCounts, fmt=fmt, header=header)
# TODO: specify dtype
