from Bio import SeqIO
import os

def loadRef(fasta):
	refs = dict()
	readsPath = '/'.join(fasta.split('/')[:-1])
	handle = open(fasta, "rU")
	for record in SeqIO.parse(handle, "fasta"):
		#print "%s with %d bases\n" % (record.id, len(record.seq))
		fast5Path = record.description.split(' ')[-1]
		refs[record.id] = (fast5Path if os.path.isabs(fast5Path) else 
				os.path.join([readsPath, fast5Path]))
	handle.close()
	return refs

def loadGenome(filename):
	genome = {}
	with open(filename, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			genome[record.id] = record
	return genome
	
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

def methylateSeq(seq):
	# simply replacing CG with MG leaves out sequences ending in C that could be methylated
	# as a result, we can't use the last base of the sequence effectively	
	return seq.replace("CG","MG")