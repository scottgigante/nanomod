#!/usr/bin/python

# Usage: python list_bam_reads.py in1.bam [in2.bam ...] in.fasta

"""
A script for listing filenames of all reads in a bamfile
"""

import pysam
import sys
import os
import shutil
import multiprocessing
import functools
import itertools
from Bio import SeqIO

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
        refs[record.id] = os.path.basename(fast5Path)
    handle.close()
    return refs

bamfiles = sys.argv[1:-1]
fastafile = sys.argv[-1]
# the bamfile should be only those reads which align to the phage. Align to a chimeric chromosome and use bamtools split -in in.bam -reference
# the fastafile should be generated using poretools

refs = loadRef(fastafile)

for bamfile in bamfiles:
    with pysam.AlignmentFile(bamfile) as bam:
        print "\n".join([refs[read.query_name] for read in bam])


