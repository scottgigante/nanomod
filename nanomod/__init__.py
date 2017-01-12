################################################################################
#                                                                              #
# __init__.py: initialise Python package.                                      #
#                                                                              #
# This file is part of Nanomod. Say something about a GNU licence?             #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

import subprocess
import os
import json

__version__ = '0.0.1'
__version_info__ = tuple([int(num) for num in __version__.split('.')])

__exe__ = {}
# list of executables required by this package
__prognames__ = ['poretools', 'bwa', 'samtools', 'nanopolish', 'nanonettrain']#, 'currennt']
__config__ = ".config.json"

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

def loadConfig():
	with open(__config__, 'r') as fh:
		config = json.load(fh)
		for key in config['exe']:
			__exe__[key] = config['exe'][key]

def saveConfig(**config):
	with open(__config__, 'w') as fh:
		json.dump(config, fp=fh, indent=4)

def init():
	if os.path.isfile(__config__):
		loadConfig()
	else:
		# initialise and check all executables
		for p in __prognames__:
			initExecutable(p)
			checkExecutable(p)
		saveConfig(exe=__exe__)