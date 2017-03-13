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
# index_modifications.py: Index where modified bases are in each read          #
# Input: fast5 reads and nanomod trained network                               #
# Output: fasta, bam, and file indexing probability of methylation per-base.   #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 19 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

from multiprocessing import Pool, Lock
from Bio import SeqIO
import pysam
import os
import json
import numpy as np
import logging

from utils import callSubProcess, multiprocessWrapper
from seq_tools import loadGenome, getSeqDiff

def getBayesPercentage(numMods, numUnmods, alpha):
    # TODO: can we use coverage vs. numMods + numUnmods to get a clever value for alpha?
    return float(numMods + alpha)/(numMods + numUnmods + alpha)

def checkModifiedPositions(options):
    modPositions = dict()
    canonGenome = loadGenome(options)
    modGenome = loadGenome(options, modified=True)
    for contig in canonGenome:
        modPositions[contig] = dict()
        modPositions[contig]['+'] = getSeqDiff(
                canonGenome[contig]['seq'], 
                modGenome[contig]['seq'])
        modPositions[contig]['-'] = getSeqDiff(
                canonGenome[contig]['reverseSeq'], 
                modGenome[contig]['reverseSeq'])    
    return modPositions, canonGenome

def checkReadModification(queryName, readPos, readLength, cigarTuples, reverse, 
        modDir):
    # TODO: check that this checks the correct position
    if reverse:
        readPos = -readPos-1 
        readPos += readLength
    startHardClip = cigarTuples[0][1] if cigarTuples[0][0] == 5 else 0
    readPos += startHardClip
    
    modFile = os.path.join(modDir, queryName)
    with open(modFile, 'rb') as handle:
        modIndex = json.load(handle)
    for base in modIndex:
        if readPos in modIndex[base]:
            return base
    return False
    
def checkMotif(readPos, sequenceMotif, motifModPos, seq, deletion=False):
    # check after modpos
    if deletion:
        readIdx = readPos
    else:
        readIdx = readPos + 1
    motifIdx = motifModPos + 1
    while motifIdx < len(sequenceMotif[0]):
        try:
            if not sequenceMotif[0][motifIdx] == seq[readIdx]:
                return False
        except IndexError:
            # ran out of read
            break
        motifIdx += 1
        readIdx += 1
    
    # check before modpos
    readIdx = readPos - 1
    motifIdx = motifModPos - 1
    while motifIdx >= 0:
        try:
            if not sequenceMotif[0][motifIdx] == seq[readIdx]:
                return False
        except IndexError:
            # ran out of read
            break
        motifIdx -= 1
        readIdx -= 1
    return True

def processBase(pileupColumn, alpha, modPositions, motifModPos, genome, contig, modDir, 
        dtype, options):
    refPos = pileupColumn.reference_pos
    
    # check if there are any modifiable bases in this position on reference
    checkFwd = False
    checkReverse = False
    refName = pileupColumn.reference_name    
    if refPos in modPositions[refName]['+']:
        checkFwd = True
    if (-refPos-1+len(genome)) in modPositions[refName]['-']:
        # TODO: we have 345000 sites when we should have double that - is reverse working?
        checkReverse = True
    if not (checkFwd or checkReverse):
        # ref genome is not modified here
        return
    
    # count modifications
    count = { 'M' : 0, # modified in correct motif
              'U' : 0, # unmodified in correct motif
              'u' : 0, # unmodified in incorrect motif
              'X' : 0, # mismatch in otherwise correct motif
              'x' : 0, # mismatch in incorrect motif
              'D' : 0, # deletion in otherwise correct motif
              'd' : 0  # deletion in incorrect motif
            }
    for read in pileupColumn.pileups:
        # check direction
        reverse = read.alignment.is_reverse
        if ((reverse and not checkReverse) or
                (not reverse and not checkFwd)):
            # wrong direction strand, nothing to check
            continue
        
        # check for mismatch / deletion
        readPos = read.query_position
        if read.is_del == 1:
            # deletion
            if checkMotif(read.query_position_or_next, options.sequenceMotif, motifModPos, read.alignment.query_sequence, deletion=True):
                count['D'] += 1
            else:
                count['d'] += 1
            continue
        if not (read.alignment.query_sequence[readPos] == genome[refPos]):
            # mismatch
            if checkMotif(readPos, options.sequenceMotif, motifModPos, read.alignment.query_sequence):
                count['X'] += 1
            else:
                count['x'] += 1
            continue
        
        try:
            if checkReadModification(read.alignment.query_name, readPos, 
                    read.alignment.query_length, read.alignment.cigartuples, 
                    reverse, modDir) is not False:
                count['M'] += 1
            else:
                if checkMotif(readPos, options.sequenceMotif, motifModPos, read.alignment.query_sequence):
                    count['U'] += 1
                else:
                    count['u'] += 1
        except IOError as e:
            logging.warning(str(e))
    modProb = getBayesPercentage(count['M'], count['U'], alpha)
    return np.array([(contig, refPos, modProb, count['M'], count['U'], count['u'], count['X'], count['x'], count['D'], count['d'])], dtype=dtype)

def processWindow(bamFile, contig, startPos, endPos, alpha, modPositions, motifModPos, 
        genome, modDir, dtype, options):
    bam = pysam.AlignmentFile(bamFile, "rb")
    output = np.empty(shape=(0,), dtype=dtype)
    for base in bam.pileup(contig, startPos, endPos):
        refPos = base.reference_pos
        if refPos >= startPos and refPos < endPos:
            # window is expanded for some reason
            baseOutput = processBase(base, alpha, modPositions, motifModPos, genome, contig, 
                    modDir, dtype, options)
            if baseOutput is not None:
                output = np.append(output, baseOutput)
    bam.close()
    return output
    
def processWindowWrapper(args):
    return multiprocessWrapper(processWindow, args)

def aggregateWindowCounts(modCounts, pos, row, forwardStep, backStep, contig, 
        aggDtype, options):
    """ forwardStep 1-based end position
    backStep 0-based start position """
    # coarse window
    countsWindow = modCounts[range(max(0, row-backStep), 
            min(modCounts.shape[0], row+forwardStep))]
    afterStart = countsWindow['f1'] >= pos - backStep
    beforeEnd = countsWindow['f1'] < pos + forwardStep
    # fine window
    countsWindow = countsWindow[np.logical_and(afterStart, beforeEnd)]
    startPos = min(countsWindow['f1'])
    endPos = max(countsWindow['f1'])
    numMods = sum(countsWindow['f3'])
    numUnmods = sum(countsWindow['f4'])
    numReads = sum(countsWindow['f5'])
    print [startPos, endPos, numReads]
    modProb = getBayesPercentage(numMods, numUnmods, numReads, 
                options.alpha)
    return np.array([(contig, startPos, endPos, modProb, numMods, numUnmods, 
            numReads)], dtype=aggDtype)

def countModifications(bamFile, modDir, options):
    # TODO: allow specifying a window to reduce time cost
    logging.info("Indexing modification sites on reference genome...")
    modPositions, genome = checkModifiedPositions(options)
    motifModPos = [i for i in xrange(len(options.sequenceMotif[0])) if options.sequenceMotif[0][i] != options.sequenceMotif[1][i]][0]
    
    # S32 string, <u4 uint32, <f4 float32, u1 uint8 <u2 uint16
    dtype = np.dtype("S32,<u4,<f4,u1,u1,u1,u1,u1,u1,u1")
    aggDtype = np.dtype("S32,<u4,<u4,<f4,<u2,<u2,<u2,<u2,<u2,<u2,<u2")
    if options.window is None or options.window < 2:
        fmt = "\t".join(["%s","%i","%lf","%i","%i","%i","%i","%i","%i","%i"])
        header = "\t".join(["Chromosome", "Position", "Modified Percentage", 
                "M", "U", "u", "X", "x", "D", "d"])
        modCounts = np.empty(shape=(0,), dtype=dtype)
    else:
        fmt = "\t".join(["%s","%i","%i","%lf","%i","%i","%i","%i","%i","%i","%i"])
        header = "\t".join(["Chromosome", "Start", "Stop",
                "Modified Percentage", "M", "U", "u", "X", "x", "D", "d"])
        modCounts = np.empty(shape=(0,), dtype=aggDtype)
    
    logging.info("Analyzing read modifications...")
    basesPerCall = 2000 
    # pileup seems to take a minimum of 2000 even if you specify less
    for contig in genome:
        contigModCounts = np.empty(shape=(0,), dtype=dtype)
        refStart = 1
        refEnd = len(genome[contig]['seq'])
        if options.region is not None:
            if contig != options.region[0]:
                continue
            else:
                refStart = max(refStart, options.region[1])
                if options.region[2] is not None:
                    refEnd = min(refEnd, options.region[2])
        
        windowArray = [i for i in xrange(refStart, refEnd, basesPerCall)]
        p = Pool(options.threads)
        argsList = [[bamFile, contig, i, min(i+basesPerCall, refEnd), 
                options.alpha, modPositions, motifModPos, genome[contig]['seq'], 
                modDir, dtype, options] for i in windowArray]
        windows = p.map(processWindowWrapper, argsList)
        for w in windows:
            contigModCounts = np.append(contigModCounts, w)
        
        if options.window is None or options.window < 2:
            modCounts = np.append(modCounts, contigModCounts)
        else:
            # aggregate calls
            backStep = (options.window - 1) / 2
            forwardStep = (options.window-1) / 2 + (options.window - 1) % 2 + 1
            for i in range(contigModCounts.shape[0]):
                pos = contigModCounts['f1'][i]
                window = aggregateWindowCounts(contigModCounts, pos, i,
                        forwardStep, backStep, contig, aggDtype, options)
                if modCounts.shape[0] == 0 or not window == modCounts[-1]:
                    # exclude duplicate rows
                    modCounts = np.append(modCounts, window)
                    
    return fmt, header, modCounts
        
