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
# seq_tools.py: perform elementary DNA sequence manipulation.                  #
#                                                                              #
# TODO: include wildcard expansion                                             #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

from Bio import SeqIO, Seq
import os
from copy import copy
import numpy as np

__canonical__ = ['A','G','C','T']
wildcards = { 'W' : ['A', 'T'],
              'S' : ['C','G'],
              'M' : ['A','C'],
              'K' : ['G','T'],
              'R' : ['A','G'],
              'Y' : ['C','T'],
              'B' : ['G','C','T'],
              'H' : ['A','C','T'],
              'H' : ['A','G','T'],
              'V' : ['A','G','C'],
              'N' : ['A','G','C','T'] }

# explode wildcards in sequence motifs
# @param sequenceMotifs Sequences to be exploded
def explodeMotifs(sequenceMotifs):
    # NOT YET IMPLEMENTED
    exit(1)

def expandAlphabet(sequenceMotif, alphabet=__canonical__):
    expanded = list(alphabet)
    for base in sequenceMotif[1]:
        if base not in expanded:
            expanded.append(base)
    return expanded

def unmodifyFasta(inFile, outFile, sequenceMotif):
    """Load a fasta file and clean it of base modifications, returning the 
    original fasta as an array of records."""
    fasta = []
    with open(inFile, "rU") as inHandle, open(outFile, "w") as outHandle:
        for record in SeqIO.parse(inHandle, "fasta"):
            fasta.append(record)
            unmodifiedRecord = copy(record)
            unmodifiedRecord.seq = unmodifySeq(unmodifiedRecord.seq, 
                    sequenceMotif)
            SeqIO.write(unmodifiedRecord, outHandle, "fasta")
    return fasta

# create a dictionary linking read names to fast5 file paths
#
# @args fasta Filename of the fasta file created by poretools
# @return Dictionary linking read names to fast5 file paths
def loadRef(fasta):
    refs = dict()
    readsPath = '/'.join(fasta.split('/')[:-1])
    handle = open(fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        #print "%s with %d bases\n" % (record.id, len(record.seq))
        fast5Path = record.description.split(' ')[-1]
        refs[record.id] = os.path.abspath(fast5Path)
    handle.close()
    return refs

# load the reference genome
#
# @args filename Fasta filename of reference genome
# @return dictionary of SeqIO records corresponding to each contig of reference
def loadGenome(options, modified=False):
    genome = {}
    with open(options.genome, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contig = { 'id' : record.id }
            reverseSeq = Seq.reverse_complement(record.seq) 
            Seq.reverse_complement(reverseSeq)
            if modified:
                record.seq = modifySeq(options, record.seq)
                reverseSeq = modifySeq(options, reverseSeq)
            contig['seq'] = record.seq
            contig['reverseSeq'] = reverseSeq
            contig['record'] = record
            genome[record.id] = contig
    return genome

# get the kmer at a particular place in the genome
# @param genome the genome dictionary from loadGenome
# @param chromosome Contig id
# @param pos Position on forward reference
# @param kmer Kmer length
# @param forward Boolean marker of forward or reverse read
# @return Kmer as a string
def getKmer(genome, chromosome, pos, kmer, forward):
    stop_pos = pos + kmer
    if forward:
        seq = genome[chromosome]['seq']
        start_pos = pos
    else:
        seq = genome[chromosome]['reverseSeq']
        start_pos = -stop_pos
        stop_pos = -pos
    return str(seq[start_pos:stop_pos])

def randomPermuteSeq(seq, alphabet, rate):
	for i in range(len(seq)):
		if np.random.rand() < rate:
			seq = list(seq)
			seq[i] = alphabet[np.random.randint(len(alphabet))]
			seq = "".join(seq)
	return seq

# apply base modification to a sequence
#
# @args seq The sequence to be modified
# @return The modified sequence
def modifySeq(options, seq):
    # TODO: we do this twice - make it a method?    
    if type(seq) == Seq.Seq:
        returnBioSeq = True
        seq = str(seq)
    else:
        returnBioSeq = False
        
    seq = seq.replace(options.sequenceMotif[0],options.sequenceMotif[1])
    return Seq.Seq(seq) if returnBioSeq else seq

# replace a pattern only if at the start of a string
# http://code.activestate.com/recipes/577252-lreplace-and-rreplace-replace-the-beginning-and-en/
def lreplace(pattern, sub, string):
    return sub + string[len(pattern):] if string.startswith(pattern) else string

# replace a pattern only if at the end of a string
def rreplace(pattern, sub, string):
    return string[:-len(pattern)] + sub if string.endswith(pattern) else string

# undo base modification to a sequence
# assume positive site recognition at either end if possible
#
# @args seq The sequence to be unmodified
# @return The unmodified sequence
def unmodifySeq(seq, sequenceMotif):
    if type(seq) == Seq.Seq:
        returnBioSeq = True
        seq = str(seq)
    else:
        returnBioSeq = False
    
    # we have to check partial matches at either end
    pattern = sequenceMotif[1]
    sub = sequenceMotif[0]
    for i in range(1, len(pattern)):
        seq = lreplace(pattern[i:], sub[i:], seq)
        seq = rreplace(pattern[:i], sub[:i], seq)
    # now replace in the middle
    seq = seq.replace(pattern,sub)
    return Seq.Seq(seq) if returnBioSeq else seq
    
def getSeqDiff(s1, s2):
    """Get all positions where two sequences differ"""
    return [i for i in xrange(len(s1)) if s1[i] != s2[i]]
