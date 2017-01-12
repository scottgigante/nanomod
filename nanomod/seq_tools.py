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
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

from Bio import SeqIO
import os

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
def loadGenome(filename):
	genome = {}
	with open(filename, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			genome[record.id] = record
	return genome

# reverse complement a sequence
#
# @args seq The sequence to be reverse-complemented
# @return Reverse-complemented sequence
def reverseComplement(seq):
	# reverse
	seq = seq[::-1]	
	# complement
	# from http://stackoverflow.com/questions/6116978/python-replace-multiple-strings
	# in order to avoid double changes, use symbols first
	repls = ('A', '1'), ('T', '2'), ('C', '3'), ('G', '4')
	seq = reduce(lambda a, kv: a.replace(*kv), repls, seq)
	repls = ('1', 'T'), ('2', 'A'), ('3', 'G'), ('4', 'C')
	return reduce(lambda a, kv: a.replace(*kv), repls, seq)

# apply CpG methylation to a sequence
#
# TODO: can we turn this into a more generalised method? how?
# TODO: maybe we should methylated the reference genome instead?
#
# @args seq The sequence to be methylated
# @return The methylated sequence
def methylateSeq(seq):
	# simply replacing CG with MG leaves out sequences ending in C that could be methylated
	# as a result, we can't use the last base of the sequence effectively	
	return seq.replace("CG","MG")

# undo CpG methylation to a sequence
#
# @args seq The sequence to be demethylated
# @return The demethylated sequence
def demethylateSeq(seq):
	if seq[-1] == 'M':
		seq = seq[:-2] + 'C'
	return seq.replace("MG","CG")