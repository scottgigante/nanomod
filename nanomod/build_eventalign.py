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
# build_eventalign.py: a simple script to run the eventalign pipeline, from    #
# Metrichor labelled reads to the final eventalign tsv file. Pipeline as       #
# described in nanopolish eventalign's README.                                 #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

import os
import glob
from multiprocessing import Pool
from functools import partial
import subprocess
import logging

from utils import callSubProcess, multiprocessWrapper, preventOverwrite
from . import __exe__

# custom use of subprocess for poretools to allow special return value for multiprocessing
#
# @args call The system call to be run
# @args options Namespace object from argparse
# @args idx The index of fasta file to be written
# @return Name of the fasta file that was written
def callPoretools(call, tempDir, idx):
    outfile = os.path.join(tempDir, "{}.fasta".format(idx))
    logging.debug("Poretools: {} files, output to {}.".format(len(call)-2, outfile))
    f = open(outfile, 'w')
    subprocess.call(call, stdout=f, shell=False)
    f.close()
    return outfile

def callPoretoolsWrapper(args):
    return multiprocessWrapper(callPoretools, args)

# runs poretools multithreaded on small chunks of fasta files
# 
# Could be deprecated by find . -name "*.fast5" | parallel "poretools fasta" 
#
# TODO: why can't this handle more than eight cores? Arbitrary.
#
# @args poretools Path to poretools
# @args options Namespace object from argparse
# @args output Path to final output fasta file
# @return None
def multithreadPoretools(poretools, tempDir, force, reads, output, filesPerCall=1000):
    # check first that we actually want to edit
    if preventOverwrite(output, force):
        return 1
        
    cwd = os.getcwd()
    os.chdir(reads)
    files = [os.path.join(cwd, reads, f) for f in glob.glob("*.fast5")]
    os.chdir(cwd)
    prefix = [poretools, "fasta", "--type", "fwd"]
    
    calls = [prefix + files[i:i + filesPerCall] for i in xrange(0, len(files),
            filesPerCall)]
    args = [[calls[i], tempDir, i] for i in range(len(calls))]
    
    pool = Pool(min(options.threads, 8)) # can't handle more than eight cores
    tempFastaList = pool.map(callPoretoolsWrapper, args)
    pool.close()
    pool.join()
    
    # write results to a single fasta
    with open(output, 'w') as outfile:
        for fa in tempFastaList:
            with open(fa) as infile:
                for line in infile:
                    outfile.write(line)
    
# main script. Build eventalign file according to nanopolish's pipeline.
#
# @args options Namespace object from argparse
# @return fastaFile Filename of fasta corresponding to reads
# @return eventalignFile Filename of eventalign tsv
def buildEventalign(options, reads, outPrefix):
    
    fastaFile = '{}.fasta'.format(outPrefix)
    
    if False:
        callSubProcess(('find {} -name "*.fast5" | parallel -j16 -X {} fasta ' 
                '--type fwd {}').format(reads, 
                __exe__['poretools'],"{}"), options.force, newFile=fastaFile, 
                outputFile=fastaFile)
    else:
        multithreadPoretools(__exe__['poretools'], options.tempDir, options.force, reads, fastaFile)
        poretoolsMaxFiles = 1000
        callSubProcess(('{} fasta --type fwd {}').format(__exe__['poretools'], 
                reads), options.force, newFile=fastaFile, outputFile=fastaFile)
    
    callSubProcess('{} index {}'.format(__exe__['bwa'], options.genome), 
            options.force, newFile="{}.bwt".format(options.genome))
            
    samFile = "{}.sam".format(outPrefix)
    callSubProcess('{} mem -x ont2d -t {} {} {}'.format(__exe__['bwa'], options.threads,
            options.genome, fastaFile), options.force, outputFile=samFile, newFile=samFile)
    
    sortedBamFile = "{}.sorted.bam".format(outPrefix)
    callSubProcess('samtools sort -o {} -O bam -@ {} -T nanomod {}'.format(sortedBamFile, options.threads, samFile), options.force, newFile = sortedBamFile)
    
    callSubProcess('{} index {}'.format(__exe__['samtools'], sortedBamFile), 
            options.force, newFile="{}.bai".format(sortedBamFile))
    
    eventalignFile = "{}.eventalign".format(outPrefix)
    callSubProcess(('{} eventalign -t {} --print-read-names -r {} -b {} -g {}' 
            ' --models {}').format(__exe__['nanopolish'], options.threads, 
            fastaFile, sortedBamFile, options.genome, options.nanopolishModels), 
            options.force, newFile=eventalignFile, outputFile=eventalignFile)
    
    return fastaFile, eventalignFile