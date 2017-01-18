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
# count_mismatch_indel.py: There must be a better way to do this. Checks indel #
# and mismatch error rates in a bam file to assess model quality.              #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 16 Jan 2017                                                            #
#                                                                              #
################################################################################

import sys
import csv
csv.field_size_limit(sys.maxsize)

def run(filename):
	with open(file, 'r') as bamfile:
		total_edit = 0
		total_score = 0
		total_insert = 0
		total_delete = 0
		total_read = total_genome = 0
		longest_read = ""
		longest_read_length = 0
		reads = 0
		fail_flags=[4]
		cigar_codes = ['S','H','M','D','I']
		bamfile.next() # header
		bamfile.next() # header
		for line in csv.reader(bamfile, delimiter='\t'):
			try:
				flag = int(line[1])
				if flag not in fail_flags:
					cigar = line[5]
					n = ""
					read_length = 0
					for i in cigar:
						if i not in cigar_codes:
							n = n + i
						elif i == 'D':
							total_delete += int(n)
							total_genome += int(n)
							n = ""
						elif i == 'I':
							total_insert += int(n)
							read_length += int(n)
							n = ""
						elif i == 'M':
							read_length += int(n)
							total_genome += int(n)
							n = ""
						else:
							n = ""
					edit_dist = int(line[11].split(':')[-1])
					score = int(line[13].split(':')[-1])
					total_edit += edit_dist
					total_score += score
					total_read += read_length
					print line[0] + "\t" + str(read_length)
					if read_length > longest_read_length:
						longest_read = line[0]
						longest_read_length = read_length
					reads += 1
			except IndexError:
				print line
	return { 'averageEditDistance' : float(total_edit)/reads, 
			 'averageAlignmentScore' : float(total_score)/reads,
			 'insertionPercent' : float(total_insert)/total_read,
			 'deletionPercent' : float(total_delete)/total_genome,
			 'averageReadLength' : float(total_read)/reads,
			 'longestRead' : longest_read,
			 'longestReadLength' : longest_read_length }
			

if __name__ == "__main__":
	args = sys.argv[1:]
	for file in args:
		print "{}: {}".format(file,run(file))

