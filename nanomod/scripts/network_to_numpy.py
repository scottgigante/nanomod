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
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 22 Apr 2017                                                            #
#                                                                              #
################################################################################

# Usage: python network_to_numpy.py network_meta.json network1.json [network2.json] [...]

import numpy as np
import json
import re
import sys
from nanonet.currennt_to_pickle import network_to_numpy

def run(network, mod_meta):
    mod = json.load(open(network, 'r'))
    mod['meta'] = mod_meta
    mod = network_to_numpy(mod)
    np.save(re.sub("\.(json|jsn|autosave)$", ".npy", network), mod)

if __name__ == "__main__":
    meta = sys.argv[1]
    networks = sys.argv[2:]
    mod_meta = json.load(open(meta, 'r'))
    for n in networks:
        run(n, mod_meta)
