#!/usr/bin/python

# Usage: python seq_length.py in.fasta | cut -f 2 | bash summary_stats.sh

from Bio import SeqIO
import sys
cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
 try:
  output_line = '%s\t%i' % \
(seq_record.id, len(seq_record))
  print output_line
 except Exception:
  pass
