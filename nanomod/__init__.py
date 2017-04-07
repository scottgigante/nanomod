################################################################################
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
# __init__.py: initialise Python package.                                      #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

import subprocess
import os
import json
from version import __version__

__version_info__ = tuple([int(num) for num in __version__.split('.')])

__exe__ = {}
# list of executables required by this package
__prognames__ = { 'train' : [
                    'poretools',
                    'bwa',
                    'samtools',
                    'nanopolish',
                    'nanonettrain'
                    ],#, 'currennt']
                  'call' : [
                      'bwa',
                      'samtools',
                      'nanonetcall'
                      ] }
__config__ = "config.json"
__modes__ = ["skip", "stay", "step", "qscore", "random"]

# define path to an executable in global dictionary
#
# @args progname name of the executable to be run
# @return None
def initExecutable(progname):
    try:
        __exe__[progname] = os.path.abspath(os.environ[progname.upper()])
    except KeyError:
        __exe__[progname] = progname

# Check we can an executable by requesting its help text
#
# @args progname name of the executable to be run
# @return None
def checkExecutable(progname):
    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.call([__exe__[progname], '-h'], stdout=devnull,
                    stderr=devnull)
    except OSError:
        raise OSError(("Cannot execute {0}, it must be in your path as '{0}' or"
                " set via the environment variable '{1}'.").format(progname,
                progname.upper()))

def loadConfig(configName):
    with open(configName, 'r') as fh:
        config = json.load(fh)
        for key in config['exe']:
            __exe__[key] = config['exe'][key]

def saveConfig(configName, **config):
    with open(configName, 'w') as fh:
        json.dump(config, fp=fh, indent=4)

def init(command):
    configName = ".".join(["",command,__config__])
    if os.path.isfile(configName):
        loadConfig(configName)
    else:
        # initialise and check all executables
        for p in __prognames__[command]:
            initExecutable(p)
            checkExecutable(p)
        saveConfig(configName, exe=__exe__)
