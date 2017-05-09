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
# expandModelAlphabet.py: builds a currennt network json for an enlarged       #
# base alphabet from an existing network: assumes that modified bases start    #
# with initialisation equal to the unmodified bases.                           #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 10 Jan 2017                                                            #
#                                                                              #
################################################################################

import json
import copy
import sys
import numpy as np
from nanonet.util import all_kmers

from seq_tools import unmodifySeq, __canonical__, expandAlphabet
from utils import loadJson, saveJson, preventOverwrite

def generateKmers(inAlpha, outAlpha, inKmerLen, outKmerLen, sequenceMotif):
    """
    Generate a list of kmers for the expanded alphabet, and a map from the modified kmers
    to their corresponding unmodified kmers, when there are any.

    :param inAlpha: List of characters, original alphabet
    :param outAlpha: List of characters, new alphabet
    :param inKmerLen: int length of kmers in input model
    :param outKmerLen: int length of kmers in output model
    :param sequenceMotif: two element array, canonical motif, modified motif

    :returns outKmers: array of strings, all kmers in new model
    :returns reverseMap: dictionary linking new kmer indices to original kmer strings where they exist, and None otherwise

    """

    # generate original kmers
    inKmers = all_kmers("".join(inAlpha), inKmerLen)
    inKmers.append('X'*inKmerLen)
    # create map of kmer to index in inKmers
    inKmersMap = {k:i for i,k in enumerate(inKmers)}
    # generate modified kmers
    outKmers = all_kmers("".join(outAlpha), outKmerLen)
    outKmers.append('X'*outKmerLen)

    # map each methylated kmer to its non-methylated equivalent
    reverseMap = dict()
    if inKmerLen > outKmerLen:
        if inKmerLen % 2 == 0:
            start = (inKmerLen - outKmerLen) / 2
            end = -((inKmerLen - outKmerLen) / 2) - (inKmerLen - outKmerLen) % 2
        else:
            start = (inKmerLen - outKmerLen) / 2 + (inKmerLen - outKmerLen) % 2
            end = -((inKmerLen - outKmerLen) / 2)
        if end == 0:
            end = inKmerLen
        for i in range(len(outKmers)):
            k = unmodifySeq(outKmers[i], sequenceMotif)
            reverseMap[i] = []
            for inKmer in inKmersMap.keys():
                if inKmer[start:end] == k:
                    reverseMap[i].append(inKmersMap[inKmer])
            if reverseMap[i] == []:
                pass
                # TODO: should we randomly initialise (rather than zero) for nonexistent kmers?
    elif inKmerLen < outKmerLen:
        raise ValueError("Cannot build model with longer kmers than input model")
    else:
        for i in range(len(outKmers)):
            k = unmodifySeq(outKmers[i], sequenceMotif)
            reverseMap[i] = []
            if inKmersMap.has_key(k):
                reverseMap[i].append(inKmersMap[k])
            else:
                pass
                # TODO: should we randomly initialise (rather than zero) for nonexistent kmers?
    return outKmers, reverseMap

def createExpandedNetwork(inNetwork, kmers, reverseMap):
    """
    Generate an expanded network dictionary based on a pre-existing network and a list of
    new kmers.

    :param inNetwork: Dictionary Pre-existing network
    :param kmers: list of strings, all possible kmers in new model
    :param reverseMap: dictionary linking new kmer indices to original kmer strings where they exist, and None otherwise

    :returns: Dictionary new network with expanded alphabet
    """
    outNetwork = copy.deepcopy(inNetwork)

    # expand softmax output and multiclass postoutput layer
    outNetwork['layers'][-2]['size'] = len(kmers)
    outNetwork['layers'][-1]['size'] = len(kmers)

    # reallocate output layer weights according to reverse map
    outputLayer = outNetwork['layers'][-2]['name']

    inputSize = len(kmers) * outNetwork['layers'][-3]['size']
    outNetwork['weights'][outputLayer]['input'] = [0] * inputSize
    for i in range(outNetwork['layers'][-3]['size']):
        for j in range(len(kmers)):
            if len(reverseMap[j]) > 0:
                outIdx = i * len(kmers) + j
                for mapping in reverseMap[j]:
                    inIdx = i * inNetwork['layers'][-3]['size'] + mapping
                    outNetwork['weights'][outputLayer]['input'][outIdx] += \
                            inNetwork['weights'][outputLayer]['input'][inIdx]
                outNetwork['weights'][outputLayer]['input'][outIdx] = outNetwork['weights'][outputLayer]['input'][outIdx] / len(reverseMap[j])

    biasSize = len(kmers)
    outNetwork['weights'][outputLayer]['bias'] = [0] * biasSize
    for i in range(biasSize):
        if len(reverseMap[i]) > 0:
            for mapping in reverseMap[i]:
                outNetwork['weights'][outputLayer]['bias'][i] += \
                    inNetwork['weights'][outputLayer]['bias'][mapping]
            outNetwork['weights'][outputLayer]['bias'][i] = outNetwork['weights'][outputLayer]['bias'][i] / len(reverseMap[i])

    return outNetwork

def integer_nth_root(val, n):
    """
    Custom integer nth root script since rounding errors make the expected val**(1.0/n)
    fail some of the time.

    :param val: int Value to be nth rooted
    :param n: int nth root

    :raises TypeError: if val or n are not integers

    :returns: int nth root of val
    """
    if not (int(val) == val and int(n) == n):
        raise TypeError("Integer nth root supports integer values only")
    ret = int(val**(1./n))
    return ret + 1 if (ret + 1) ** n == val else ret

def getInKmerLen(inNetwork, alphabet=__canonical__):
    """
    Check kmer length used in original model by inferring from number of output units

    :param inNetwork: Dictionary Pre-existing network
    :param alphabet: list of characters DNA alphabet used in network

    :returns: int Length of kmers in pre-existing network.
    """
    numOutputs = inNetwork['layers'][-2]['size']

    # one of them is badKmer
    numOutputs = numOutputs - 1

    # kmerLen^len(alphabet) == numOutputs
    kmerLen = integer_nth_root(numOutputs, len(alphabet))

    return kmerLen

def expandModelAlphabet(inFilename, outFilename, outKmerLen, sequenceMotif=["",""], force=True):
    """
    Expand nanonet model file to include all possible outputs from an expanded model

    :param inFilename: String path to nanonet trained network JSON file
    :param outFilename: String path to output JSON file
    :param outKmerLen: int kmer length for output network (must be <= network kmer length)
    :param sequenceMotif: two element array of strings, canonical motif, modified motif
    :param force: Boolean value, force creation of extant file or not
    """
    if preventOverwrite(outFilename, force):
        return 1
    inNetwork = loadJson(inFilename)
    inKmerLen = getInKmerLen(inNetwork)
    alphabet = expandAlphabet(sequenceMotif, __canonical__)
    kmers, reverseMap = generateKmers(__canonical__, alphabet, inKmerLen, outKmerLen,
            sequenceMotif)
    outNetwork = createExpandedNetwork(inNetwork, kmers, reverseMap)
    saveJson(outFilename, outNetwork)

if __name__ == "__main__":
    """
    Command line interface.

    Usage: python expand_model_alphabet.py in_network.json out_network.json <kmer_length>
    <canonical_motif> <modified_motif>
    """
    args = sys.argv
    try:
        expandModelAlphabet(args[1], args[2], int(args[3]), [args[4], args[5]])
    except IndexError:
        # no sequence motif specified
        expandModelAlphabet(args[1], args[2], int(args[3]))
