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
#
# TODO: store in memory rather than saving to disk
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 19 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

import os
import json
import logging
from multiprocessing import Pool
import re

from seq_tools import unmodifyFasta, __canonical__
from utils import makeDir, multiprocessWrapper, preventOverwrite

def indexRead(pattern, record, outputDir, force):
    outputFile = os.path.join(outputDir, record.id)
    if preventOverwrite(outputFile, force):
        return
    modIndex = dict()
    for m in pattern.finditer(str(record.seq)):
        i = m.start()
        modIndex[i] = record.seq[i]
    if modIndex != {}:
        with open(outputFile, 'wb') as handle:
            json.dump(modIndex, fp=handle)

def indexReadWrapper(args):
    return multiprocessWrapper(indexRead, args)

def indexAndCleanModifications(fastaFile, options):
    unmodifiedFastaFile = "{}.unmodified.fasta".format(options.outPrefix)
    outputDir = "{}.modIndex".format(options.outPrefix)
    if preventOverwrite(unmodifiedFastaFile, options.force): # TODO: what if indexing crashes?
        return unmodifiedFastaFile, outputDir

    makeDir(outputDir)

    logging.info("Cleaning fasta of modifications...")
    fasta = unmodifyFasta(fastaFile, unmodifiedFastaFile, options.sequenceMotif, options.threads)

    logging.info("Indexing modifications on reads...")
    p = Pool(options.threads)
    pattern = re.compile("[^{}]".format("|".join(__canonical__)))
    p.imap_unordered(indexReadWrapper, ([pattern, i, outputDir, options.force] for i in fasta))
    p.close()
    p.join()
    return unmodifiedFastaFile, outputDir



