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
from nanonet.util import all_kmers
import copy
import sys
import numpy as np

from seq_tools import unmodifySeq, __canonical__, expandAlphabet
from utils import loadJson, saveJson, preventOverwrite

def generateKmers(inAlpha, outAlpha, kmerLen, sequenceMotif):
    """
    Generate a list of kmers for the expanded alphabet, and a map from the modified kmers
    to their corresponding unmodified kmers, when there are any.

    :param inAlpha: List of characters, original alphabet
    :param outAlpha: List of characters, new alphabet
    :param kmerLen: int length of kmers
    :param sequenceMotif: two element array, canonical motif, modified motif

    :returns outKmers: array of strings, all kmers in new model
    :returns reverseMap: dictionary linking new kmer indices to original kmer strings where they exist, and None otherwise

    TODO: should we randomly initialise (rather than zero) for nonexistent kmers?
    """
    badKmer = 'X'*kmerLen

    # generate original kmers
    inKmers = all_kmers("".join(inAlpha), kmerLen)
    inKmers.append(badKmer)
    # create map of kmer to index in inKmers
    inKmersMap = {k:i for i,k in enumerate(inKmers)}
    # generate modified kmers
    outKmers = all_kmers("".join(outAlpha), kmerLen)
    outKmers.append(badKmer)

    # map each methylated kmer to its non-methylated equivalent
    reverseMap = dict()
    for i in range(len(outKmers)):
        k = unmodifySeq(outKmers[i], sequenceMotif)
        if inKmersMap.has_key(k):
            reverseMap[i] = inKmersMap[k]
        else:
            # kmer doesn't exist
            #reverseMap[i] = np.random.randint(len(inKmers))
            reverseMap[i] = None
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
            if reverseMap[j] is not None:
                inIdx = i * inNetwork['layers'][-3]['size'] + reverseMap[j]
                outIdx = i * outNetwork['layers'][-3]['size'] + j
                outNetwork['weights'][outputLayer]['input'][outIdx] = \
                        inNetwork['weights'][outputLayer]['input'][inIdx]

    biasSize = len(kmers)
    outNetwork['weights'][outputLayer]['bias'] = [0] * biasSize
    for i in range(biasSize):
        if reverseMap[i] is not None:
            outNetwork['weights'][outputLayer]['bias'][i] = \
                    inNetwork['weights'][outputLayer]['bias'][reverseMap[i]]

    return outNetwork

def run(inFilename, outFilename, kmer, sequenceMotif, force=True):
    """
    Expand nanonet model file to include all possible outputs from an expanded model

    :param inFilename: String path to nanonet trained network JSON file
    :param outFilename: String path to output JSON file
    :param kmer: int kmer length (must match network kmer length
    :param sequenceMotif: two element array of strings, canonical motif, modified motif
    :param force: Boolean value, force creation of extant file or not

    # TODO: can we infer kmer length?
    """
    if preventOverwrite(outFilename, force):
        return 1
    inNetwork = loadJson(inFilename)
    alphabet = expandAlphabet(sequenceMotif, __canonical__)
    kmers, reverseMap = generateKmers(__canonical__, alphabet, kmer,
            sequenceMotif)
    outNetwork = createExpandedNetwork(inNetwork, kmers, reverseMap)
    saveJson(outFilename, outNetwork)

def expandModelAlphabet(options):
    """
    Expand nanonet model file to include all possible outputs from an expanded model

    :param options: Namespace object from argparse
    """
    run(options.currenntTemplate, options.expandedTemplate, options.kmer,
            options.sequenceMotif, options.force)

if __name__ == "__main__":
    """
    Command line interface.

    Usage: python expand_model_alphabet.py in_network.json out_network.json <kmer_length>
    <canonical_motif> <modified_motif>
    """
    args = sys.argv
    run(args[1], args[2], int(args[3]), [args[4], args[5]])
