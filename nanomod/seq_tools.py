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
import re
import multiprocessing
import functools
import itertools

__canonical__ = ['A','G','C','T']
__wildcards__ = { 'W' : ['A', 'T'],
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

def explodeMotifs(sequenceMotifs):
    """
    Explode wildcards in sequence motifs

    :param sequenceMotifs: Sequences to be exploded

    :returns: List of sequence motifs, exploded by DNA wildcards

    >>> explodeMotifs(['CCWGG', 'CMWGG'])
    [['CCAGG', 'CMAGG'], ['CCTGG','CMTGG']]

    NOT YET IMPLEMENTED: how do we deal with multiple sequence motifs???

    TODO: implement wildcard expansion
    """
    exit(1)

def expandAlphabet(sequenceMotif, alphabet=__canonical__):
    """
    Expand nucleotide alphabet based on sequence motif

    :param sequenceMotif: Two element array of strings: canonical motif, modified motif
    :param alphabet: Canonical DNA alphabet to use as a starting point

    >>> expandAlphabet(["CG", "MG"], ["A","C","G","T"])
    ["A","C","G","M","T"]
    """
    expanded = list(alphabet)
    for base in sequenceMotif[1]:
        if base not in expanded:
            expanded.append(base)
    return sorted(expanded)

def unmodifyRecord(record, sequenceMotif):
    unmodifiedRecord = copy(record)
    unmodifiedRecord.seq = unmodifySeq(unmodifiedRecord.seq, sequenceMotif)
    return record, unmodifiedRecord

def unmodifyFasta(inFile, outFile, sequenceMotif, threads):
    """
    Load a fasta file and clean it of base modifications, returning the
    original fasta as an array of records.

    :param inFile: String path to original fasta file
    :param outFile: String path to output fasta file
    :param sequenceMotif: Two-element array of strings, canonical motif, modified motif to clean

    :returns: list of records in original fasta file
    """
    pool = multiprocessing.Pool(threads)
    with open(inFile, "rU") as inHandle, open(outFile, "w") as outHandle:
        unmodifiedFasta = pool.imap_unordered(
                functools.partial(unmodifyRecord, sequenceMotif=sequenceMotif),
                SeqIO.parse(inHandle, "fasta"))
        for record, unmodifiedRecord in unmodifiedFasta:
            SeqIO.write(unmodifiedRecord, outHandle, "fasta")
            yield record
    pool.close()
    pool.join()

def loadRef(fasta):
    """
    Create a dictionary linking read names to fast5 file paths

    :param fasta: Filename of the fasta file created by poretools

    :returns: Dictionary linking read names to fast5 file paths
    """
    refs = dict()
    readsPath = '/'.join(fasta.split('/')[:-1])
    handle = open(fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        #print "%s with %d bases\n" % (record.id, len(record.seq))
        fast5Path = record.description.split(' ')[-1]
        refs[record.id] = os.path.abspath(fast5Path)
    handle.close()
    return refs

def loadGenome(genomeFile, sequenceMotif, modified=False):
    """
    Load a reference genome and store modified and reversed copies

    :param genomeFile: string Path to genome fasta file
    :param sequenceMotif: Two-element array of strings, canonical motif, modified motif
    :param modified: Boolean, whether or not the sequence is modified

    :returns: dictionary of dictionaries, each of which contains a SeqIO record, forward (perhaps modified) sequence and reversed sequence corresponding to each contig of reference
    """
    genome = {}
    with open(genomeFile, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contig = { 'id' : record.id }
            reverseSeq = Seq.reverse_complement(record.seq)
            Seq.reverse_complement(reverseSeq)
            if modified:
                record.seq = modifySeq(record.seq, sequenceMotif)
                reverseSeq = modifySeq(reverseSeq, sequenceMotif)
            contig['seq'] = record.seq
            contig['reverseSeq'] = reverseSeq
            contig['record'] = record
            genome[record.id] = contig
    return genome

def getKmer(genome, chromosome, pos, kmer, forward):
    """
    Get the kmer at a particular place in the genome.

    :param genome: the genome dictionary from loadGenome
    :param chromosome: string Contig id
    :param pos: int Position on forward reference
    :param kmer: int Kmer length
    :param forward: Boolean marker of forward or reverse read

    :returns: Desired kmer as a string
    """
    stop_pos = pos + kmer
    if forward:
        seq = genome[chromosome]['seq']
        start_pos = pos
    else:
        seq = genome[chromosome]['reverseSeq']
        start_pos = -stop_pos
        stop_pos = -pos
    return str(seq[start_pos:stop_pos])

def randomPermuteSeq(seq, alphabet=__canonical__, rate=0.001):
    """
    Randomly permute bases in a sequence to add noise and decrease overfitting.

    :param seq: string Sequence to be permuted
    :param alphabet: list of characters Alphabet to be drawn from
    :param rate: float Rate at which to randomly permute

    :returns: string The modified sequence
    """
    for i in range(len(seq)):
        if np.random.rand() < rate:
            seq = list(seq)
            seq[i] = alphabet[np.random.randint(len(alphabet))]
            seq = "".join(seq)
    return seq

def convertSeqToString(seq):
    """
    Check if a sequence is a string or Bio.Seq, and return the sequence as a string.

    :param seq: string or Bio.Seq.Seq The sequence to be converted

    :returns seq: string The converted string
    :returns converted: boolean True if the seq was converted, false otherwise
    """
    if type(seq) == Seq.Seq:
        converted = True
        seq = str(seq)
    else:
        converted = False
    return seq, converted

def modifySeq(seq, sequenceMotif):
    """
    Apply base modification to a sequence

    :param seq: string The sequence to be modified
    :param sequenceMotif: Two-element array of strings: canonical motif, modified motif

    :returns: string The modified sequence

    >>> modifySeq("CAGCGT", ["CG", "MG"])
    "CAGMGT"
    """
    seq, converted = convertSeqToString(seq)

    seq = seq.replace(sequenceMotif[0],sequenceMotif[1])
    return Seq.Seq(seq) if converted else seq

def unmodifySeq(seq, sequenceMotif):
    """
    Undo base modification to a sequence
    Assume positive site recognition at either end if possible


    :param seq: string The sequence to be unmodified
    :param sequenceMotif: Two-element array of strings: canonical motif, modified motif

    :returns: string The unmodified sequence

    >>> unmodifySeq("CAGMGM", ["CG", "MG"])
    "CAGCGC"
    """
    seq, converted = convertSeqToString(seq)

    # we have to check partial matches at either end
    pattern = sequenceMotif[1]
    sub = sequenceMotif[0]
    for i in range(1, len(pattern)):
        seq = re.sub("^" + pattern[i:], sub[i:], seq)
        seq = re.sub(pattern[:i] + "$", sub[:i], seq)
    # now replace in the middle
    seq = seq.replace(pattern,sub)
    return Seq.Seq(seq) if converted else seq

def getSeqDiff(s1, s2):
    """
    Get all positions where two sequences differ for sequences of equal length

    :param s1: string Sequence to be compared
    :param s2: string Another sequence to be compared

    :raises IndexError: if sequences are not the same length

    :returns: list of integer positions where sequences differ
    """
    if len(s1) != len(s2):
        raise IndexError("Sequences must be of same length")
    for i in xrange(len(s1)):
        if s1[i] != s2[i]:
            yield i
