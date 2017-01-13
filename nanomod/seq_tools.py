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
				record.seq = Seq.Seq(modifySeq(options, str(record.seq)))
				reverseSeq = Seq.Seq(modifySeq(options, str(reverseSeq)))
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

# apply CpG methylation to a sequence
#
# TODO: can we turn this into a more generalised method? how?
# TODO: maybe we should methylated the reference genome instead?
#
# @args seq The sequence to be methylated
# @return The methylated sequence
def modifySeq(options, seq):
	# simply replacing CG with MG leaves out sequences ending in C that could be methylated
	# as a result, we can't use the last base of the sequence effectively	
	return seq.replace(options.sequenceMotif[0],options.sequenceMotif[1])

# replace a pattern only if at the start of a string
# http://code.activestate.com/recipes/577252-lreplace-and-rreplace-replace-the-beginning-and-en/
def lreplace(pattern, sub, string):
    return sub + string[len(pattern):] if string.startswith(pattern) else string

# replace a pattern only if at the end of a string
def rreplace(pattern, sub, string):
    return string[:-len(pattern)] + sub if string.endswith(pattern) else string

# undo CpG methylation to a sequence
# assume positive site recognition at either end if possible
# allows wildcards UY
#
# @args seq The sequence to be demethylated
# @return The demethylated sequence
def unmodifySeq(seq, sequenceMotif):
	# we have to check partial matches at either end
	pattern = sequenceMotif[1]
	sub = sequenceMotif[0]
	for i in range(1, len(pattern)):
		seq = lreplace(pattern[i:], sub[i:], seq)
		seq = rreplace(pattern[:i], sub[:i], seq)
	return seq.replace(pattern,sub)
