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
# TODO: put modProb back in
# TODO: remove check for c
#!/usr/bin/env python
from __future__ import print_function

from multiprocessing import Pool, Lock
from Bio import SeqIO
import pysam
import os
import json
import numpy as np
import logging
import itertools

from utils import multiprocessWrapper
from seq_tools import loadGenome, getSeqDiff

def getBayesPercentage(numMods, numUnmods, alpha):
    """
    A naive Bayesian estimator of likelihood of modification at this site. Will be deprecated by a GLM.

    :param numMods: int Number of reads observed as modified at this site
    :param numUnmods: int Number of reads observed as unmodified at this site
    :param alpha: Hyperparameter for beta prior
    """
    # TODO: implement glm solution instead
    return float(numMods + alpha)/(numMods + numUnmods + 2 * alpha)

def checkModifiedPositions(genome, sequenceMotif):
    """
    Find all positions on genome where modification motif appears.

    # TODO: this is slow for large genomes and could be easily parallelised

    :param options: Namespace object from argparse

    :returns modPositions: Dictionary of dictionaries: one dict for each contig, each of which has one dict for forward and one for reverse, each of which is a list of integer positions of modified motifs
    :returns canonGenome: Seq.Seq Unmodified genome sequence
    """
    modPositions = dict()
    canonGenome = loadGenome(genome, sequenceMotif)
    modGenome = loadGenome(genome, sequenceMotif, modified=True)
    for contig in canonGenome:
        modPositions[contig] = {'+' : dict(), '-' : dict()}
        for i in getSeqDiff(canonGenome[contig]['seq'],
                modGenome[contig]['seq']):
            modPositions[contig]['+'][i] = True
        for i in getSeqDiff(canonGenome[contig]['reverseSeq'],
                modGenome[contig]['reverseSeq']):
            modPositions[contig]['-'][i] = True
    return modPositions, canonGenome

def checkReadModification(queryName, readPos, readLength, cigarTuples, reverse,
        modDir, modIndices):
    """
    Check if a read is modified at a point.

    :param queryName: string Read unique identifier
    :param readPos: int Genomic position on read
    :param readLength: int Length of read
    :param cigarTuples: list of tuples of ints, converted cigar string from BAM record
    :param reverse: boolean, True if read is reverse read, False otherwise
    :param modDir: string Path to folder of records marking where reads were modified

    :returns: Boolean, True if modified, False otherwise
    """
    # TODO: check that this checks the correct position
    if reverse:
        readPos = readLength-readPos-1
        startHardClip = cigarTuples[-1][1] if cigarTuples[-1][0] == 5 else 0
    else:
        startHardClip = cigarTuples[0][1] if cigarTuples[0][0] == 5 else 0
    readPos += startHardClip
    try:
        modIndex = modIndices[queryName]
    except KeyError:
        try:
            modFile = os.path.join(modDir, queryName)
            with open(modFile, 'rb') as handle:
                modIndex = json.load(handle)
            modIndices[queryName] = modIndex
        except IOError:
            # no file, read has no mods
            modIndices[queryName] = {}
            return False, modIndices
    try:
        return modIndex[str(readPos)], modIndices
    except KeyError:
        # not modified
        return False, modIndices

def checkMotif(readPos, sequenceMotif, motifModPos, seq, deletion=False):
    """
    Check whether the motif is correctly or incorrectly called around a modification.

    :param readPos: int Position on the read
    :param sequenceMotif: two-element list of strings, canonical motif, modified motif
    :param motifModPos: int Position of modified base on motif
    :param seq: string Read sequence
    :param deletion: Boolean value, true if modified base is deleted in read, false otherwise

    :returns: Boolean value, true if motif (excluding modified base) is called correctly, false otherwise

    >>> checkMotif(5, ["CATG", "CMTG"], 1, "AAAACATGGGG", False)
    True

    >>> checkMotif(5, ["CATG", "CMTG"], 1, "AAAACMTGGGG", False)
    True

    >>> checkMotif(5, ["CATG", "CMTG"], 1, "AAAAGATGGGG", False)
    False

    >>> checkMotif(5, ["CATG", "CMTG"], 1, "AAAACTGGGG", True)
    True

    >>> checkMotif(5, ["CATG", "CMTG"], 1, "AAAACGT", False)
    True
    """
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

def processBase(pileupColumn, alpha, modPositions, motifModPos, genome, contig, refPos,
        modDir, modIndices, dtype, sequenceMotif):
    """
    Process modification motifs over a genomic window.

    :param pileupColumn: pysam.PileupColumn Pileup of reads at genomic position
    :param modPositions: list of int List of all genomic locations with modification motif
    :param motifModPos: int Position containing modification within motif
    :param genome: Seq.Seq forward sequence of reference genome
    :param contig: string Contig / chromosome name
    :param modDir: string Path to directory listing all modifications per read
    :param dtype: numpy.dtype Data type for output ndarray
    :param options: Namespace object from argparse

    :returns: numpy.ndarray Position, counts, and predicted modification percentage at this location, or None if no motif exists
    """

    # check if there are any modifiable bases in this position on reference
    checkFwd = False
    checkReverse = False
    try:
        modPositions[contig]['+'][refPos]
        checkFwd = True
    except KeyError:
        pass
    try:
        modPositions[contig]['-'][-refPos-1+len(genome)]
        checkReverse = True
    except KeyError:
        pass
    if not (checkFwd or checkReverse):
        # ref genome is not modified here
        return None, modIndices

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
            if checkMotif(read.query_position_or_next, sequenceMotif, motifModPos, read.alignment.query_sequence, deletion=True):
                count['D'] += 1
            else:
                count['d'] += 1
            continue
        if not (read.alignment.query_sequence[readPos] == genome[refPos]):
            # mismatch
            if checkMotif(readPos, sequenceMotif, motifModPos, read.alignment.query_sequence):
                count['X'] += 1
            else:
                count['x'] += 1
            continue

        mod, modIndices = checkReadModification(read.alignment.query_name, readPos,
                read.alignment.query_length, read.alignment.cigartuples,
                reverse, modDir, modIndices)
        if mod is not False:
            count['M'] += 1
        else:
            if checkMotif(readPos, sequenceMotif, motifModPos, read.alignment.query_sequence):
                count['U'] += 1
            else:
                count['u'] += 1
    #modProb = getBayesPercentage(count['M'], count['U'], alpha)
    #return np.array([(contig, refPos, modProb, count['M'], count['U'], count['u'], count['X'], count['x'], count['D'], count['d'])], dtype=dtype)
    return np.array([(contig, refPos, count['M'], count['U'], count['u'], count['X'], count['x'], count['D'], count['d'])], dtype=dtype), modIndices

def processWindow(bamFile, contig, startPos, endPos, alpha, modPositions, motifModPos,
        genome, modDir, dtype, sequenceMotif):
    """
    Process modification motifs over a genomic window.

    :param bamFile: string Path to sorted indexed bam file
    :param contig: string Contig / chromosome name
    :param startPos: int Window genomic location start (0-based)
    :param endPos: int Window genomic location end (1-based)
    :param alpha: float Hyperparamter for Beta prior
    :param modPositions: list of int List of all genomic locations with modification motif
    :param motifModPos: int Position containing modification within motif
    :param genome: Seq.Seq forward sequence of reference genome
    :param modDir: string Path to directory listing all modifications per read
    :param dtype: numpy.dtype Data type for output ndarray
    :param options: Namespace object from argparse

    :returns: numpy.ndarray Position, counts, and predicted modification percentage at each motif location
    """
    bam = pysam.AlignmentFile(bamFile, "rb")
    output = np.empty(shape=(0,), dtype=dtype)
    modIndices = {}
    for base in bam.pileup(contig, startPos, endPos):
        refPos = base.reference_pos
        if refPos >= startPos and refPos < endPos:
            # window is expanded for some reason
            baseOutput, modIndices = processBase(base, alpha, modPositions, motifModPos, genome, contig, refPos,
                    modDir, modIndices, dtype, sequenceMotif)
            if baseOutput is not None:
                output = np.append(output, baseOutput)
    bam.close()
    return output

def processWindowWrapper(args):
    """
    Multiprocessing wrapper for processWindow

    :param args: List of arguments for processWindow

    :returns: return value from processWindow
    """
    return multiprocessWrapper(processWindow, args)

def aggregateWindowCounts(modCounts, pos, row, forwardStep, backStep, contig,
        aggDtype, alpha):
    """
    Aggregate counts over a window of multiple bases.

    :param modCounts: numpy ndarray Single-base modification count output
    :param pos: int Genomic location at which to start window
    :param row: int Data location (motif site index) at which to start window
    :param backStep: int 0-based start position before pos
    :param forwardStep: 1-based end position after pos
    :param contig: string Contig / chromosome name
    :param addDtype: numpy.dtype Data type for output ndarray
    :param options: Namespace object from argparse

    :returns: numpy.ndarray Modification counts aggregated over windows
    """
    # coarse window
    countsWindow = modCounts[range(max(0, row-backStep),
            min(modCounts.shape[0], row+forwardStep))]
    afterStart = countsWindow['f1'] >= pos - backStep
    beforeEnd = countsWindow['f1'] < pos + forwardStep
    # fine window
    countsWindow = countsWindow[np.logical_and(afterStart, beforeEnd)]
    startPos = min(countsWindow['f1'])
    endPos = max(countsWindow['f1'])
    count = { 'M' : 0, # modified in correct motif
              'U' : 0, # unmodified in correct motif
              'u' : 0, # unmodified in incorrect motif
              'X' : 0, # mismatch in otherwise correct motif
              'x' : 0, # mismatch in incorrect motif
              'D' : 0, # deletion in otherwise correct motif
              'd' : 0  # deletion in incorrect motif
            }
    for base, column in zip(count.keys(), ['f{}'.format(i+3) for i in range(len(count))]):
        count[base] = sum(countsWindow[column])
    modProb = getBayesPercentage(numMods, numUnmods, numReads,
                alpha)
    return np.array([tuple(itertools.chain([contig, startPos, endPos, modProb], count.values()))], dtype=aggDtype)

def countModifications(bamFile, modDir, options):
    """
    Count modifications at recognised motif sites along the genome.

    :param bamFile: string Path to sorted indexed bam file
    :param modDir: string Path to directory listing all modifications per read
    :param options: Namespace object from argparse

    :returns fmt: numpy.savetxt format for output
    :returns header: string Header for output tsv
    :returns modCounts: numpy ndarray Position, counts, and predicted modification percentage at each motif location

    # TODO: use consistent output for window of 1 and window of X
    """
    logging.info("Indexing modification sites on reference genome...")
    modPositions, genome = checkModifiedPositions(options.genome, options.sequenceMotif)
    motifModPos = [i for i in xrange(len(options.sequenceMotif[0])) if options.sequenceMotif[0][i] != options.sequenceMotif[1][i]][0]

    # S32 string, <u4 uint32, <f4 float32, u1 uint8 <u2 uint16
    #dtype = np.dtype("S32,<u4,<f4,u1,u1,u1,u1,u1,u1,u1")
    dtype = np.dtype("S32,<u4,u1,u1,u1,u1,u1,u1,u1")
    #aggDtype = np.dtype("S32,<u4,<u4,<f4,<u2,<u2,<u2,<u2,<u2,<u2,<u2")
    aggDtype = np.dtype("S32,<u4,<u4,<u2,<u2,<u2,<u2,<u2,<u2,<u2")
    if options.window is None or options.window < 2:
        #fmt = "\t".join(["%s","%i","%lf","%i","%i","%i","%i","%i","%i","%i"])
        fmt = "\t".join(["%s","%i","%i","%i","%i","%i","%i","%i","%i"])
        #header = "\t".join(["Chromosome", "Position", "Modified Percentage",
        #        "M", "U", "u", "X", "x", "D", "d"])
        header = "\t".join(["Chromosome", "Position",
                "M", "U", "u", "X", "x", "D", "d"])
        modCounts = np.empty(shape=(0,), dtype=dtype)
    else:
        #fmt = "\t".join(["%s","%i","%i","%lf","%i","%i","%i","%i","%i","%i","%i"])
        fmt = "\t".join(["%s","%i","%i","%i","%i","%i","%i","%i","%i","%i"])
        #header = "\t".join(["Chromosome", "Start", "Stop",
        #        "Modified Percentage", "M", "U", "u", "X", "x", "D", "d"])
        header = "\t".join(["Chromosome", "Start", "Stop",
                "M", "U", "u", "X", "x", "D", "d"])
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

        q, r = divmod(refEnd-refStart, options.threads)
        windowArray = [q*i + min(i, r) + refStart for i in xrange(options.threads + 1)]
        p = Pool(options.threads)
        argsList = [[bamFile, contig, windowArray[i], windowArray[i+1],
                options.alpha, modPositions, motifModPos, genome[contig]['seq'],
                modDir, dtype, options.sequenceMotif] for i in range(len(windowArray)-1)]
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
                        forwardStep, backStep, contig, aggDtype, options.alpha)
                if modCounts.shape[0] == 0 or not window == modCounts[-1]:
                    # exclude duplicate rows
                    modCounts = np.append(modCounts, window)

    return fmt, header, modCounts

