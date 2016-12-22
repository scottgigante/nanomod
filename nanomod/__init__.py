import subprocess
import os

__version__ = '0.0.1'
__version_info__ = tuple([int(num) for num in __version__.split('.')])

__exe__ = {}
__prognames__ = ['poretools', 'bwa', 'samtools', 'nanopolish', 'nanonettrain']#, 'currennt']

def initExecutable(progname):
	# define path to an executable in global dictionary
	try:
		__exe__[progname] = os.path.abspath(os.environ[progname.upper()])
	except KeyError:
		__exe__[progname] = progname

def checkExecutable(progname):
    # Check we can an executable
    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.call([__exe__[progname], '-h'], stdout=devnull, 
            		stderr=devnull)
    except OSError:
        raise OSError(("Cannot execute {0}, it must be in your path as '{0}' or" 
        		" set via the environment variable '{1}'.").format(progname, 
        		progname.upper()))

for p in __prognames__:
	initExecutable(p)
	checkExecutable(p)