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
# create_serial_net.py: modifies the weights of a currennt json network to be  #
# integers in ascending order. Used to facilitate comparison pre- and post-    #
# conversion using pickle_to_currennt.py                                       #
#                                                                              #
# TODO: use argparse                                                           #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 10 Jan 2017                                                            #
#                                                                              #
################################################################################

import json
import sys

def reweightArray(net, layer, wgts):
	net['weights'][layer][wgts] = range(len(net['weights'][layer][wgts]))
	return net

def reweightLayer(net, layer):
	weights = ['input','bias','internal']
	for w in weights:
		net = reweightArray(net, layer, w)
	return net

def reweightNet(net):
	for k in net['weights'].keys():
		net = reweightLayer(net, k)
	return net

def readNet(filename):
	with open(filename, "r") as fh:
		net = json.load(fh)
	return net

def writeNet(net, filename):
	with open(filename, "w") as fh:
		json.dump(net, fp=fh, indent=4)

# main
in_filename = sys.argv[1]
out_filename = "{}.serial.json".format(filename) if len(sys.argv) < 2 else sys.argv[2]
net = readNet(filename)
net = reweightNet(net)
writeNet(net, out_filename)
