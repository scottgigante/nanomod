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

from seq_tools import unmodifySeq
from utils import loadJson, saveJson, callSubProcess

__alpha__ = ['A','G','C','T']

def expandAlphabet(alpha, sequenceMotif):
	expanded = list(alpha)
	for base in sequenceMotif[1]:
		if base not in expanded:
			expanded.append(base)
	return expanded

def generateKmers(inAlpha, outAlpha, kmerLen, sequenceMotif):
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
			reverseMap[i] = None
	return outKmers, reverseMap

def createExpandedNetwork(inNetwork, kmers, reverseMap):
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

def expandModelAlphabet(options):
	if callSubProcess("touch {}".format(options.expandedTemplate),
			options, newFile=options.expandedTemplate) == 1:
		return 1
	inNetwork = loadJson(options.currenntTemplate)
	alphabet = expandAlphabet(__alpha__, options.sequenceMotif)
	kmers, reverseMap = generateKmers(__alpha__, alphabet, options.kmer, 
			options.sequenceMotif)
	outNetwork = createExpandedNetwork(inNetwork, kmers, reverseMap)
	saveJson(options.expandedTemplate, outNetwork)