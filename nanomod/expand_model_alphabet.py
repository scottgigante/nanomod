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
# expand_model_alphabet.py: builds a currennt network json for an enlarged     #
# base alphabet from an existing network: assumes that modified bases start    #
# with initialisation equal to the unmodified bases.                           #
#                                                                              #
# TODO: Make methylation / demethylation more general.                         #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 10 Jan 2017                                                            #
#                                                                              #
################################################################################

import json
from nanonet.util import all_kmers

from seq_tools import demethylateSeq

kmer_len = 5

model_orig = "models/r9_template.json"
model_new = "models/r9_template_methyl.json"

alphabet_orig = "AGCT"
alphabet_new = "AGCTM"

# generate kmers
bad_kmer = 'X'*kmer_len
kmers_orig = all_kmers(alphabet_orig, kmer_len)
kmers_orig.append(bad_kmer)
kmers_orig_map = {k:i for i,k in enumerate(kmers_orig)}
kmers_new = all_kmers(alphabet_new, kmer_len)
kmers_new.append(bad_kmer)

# map each methylated kmer to its non-methylated equivalent
reverse_map = dict()
for i in range(len(kmers_new)):
	k = demethylateSeq(kmers_new[i])
	if kmers_orig_map.has_key(k):
		reverse_map[i] = kmers_orig_map[k]
	else:
		# kmer doesn't exist
		reverse_map[i] = None

# load json
with open(model_orig, 'r') as in_fh:
	in_network = json.load(in_fh)
	
assert(in_network['layers'][-2]['size'] == len(kmers_orig))
assert(in_network['layers'][-1]['size'] == len(kmers_orig))
out_network = in_network.copy()

# expand softmax output and multiclass postoutput layer
out_network['layers'][-2]['size'] = len(kmers_new)
out_network['layers'][-1]['size'] = len(kmers_new)

# reallocate output layer weights according to reverse map
output_layer = out_network['layers'][-2]['name']

input_size = len(kmers_new) * out_network['layers'][-3]['size']
out_network['weights'][output_layer]['input'] = [0] * input_size
for i in range(out_network['layers'][-3]['size']):
	for j in range(len(kmers_new)):
		if reverse_map[j] is not None:
			in_idx = i * in_network['layers'][-3]['size'] + reverse_map[j]
			out_idx = i * out_network['layers'][-3]['size'] + j
			out_network['weights'][output_layer]['input'][out_idx] = \
					in_network['weights'][output_layer]['input'][in_idx]

bias_size = len(kmers_new)
out_network['weights'][output_layer]['bias'] = [0] * bias_size
for i in range(bias_size):
	if reverse_map[i] is not None:
		out_network['weights'][output_layer]['bias'][i] = \
				in_network['weights'][output_layer]['bias'][reverse_map[i]]

with open(model_new, 'w') as out_fh:
	json.dump(out_network, fp=out_fh, indent=4)