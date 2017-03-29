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
import fnmatch

from utils import callSubProcess, multiprocessWrapper, preventOverwrite
from . import __exe__

def callPoretools(call, tempDir, idx):
    """
    Custom use of subprocess for poretools to allow special return value for
    multiprocessing

    :param call: The system call to be run
    :param tempDir: String path to directory for temporary files
    :param idx: int index of fasta file to be written to make unique naming system

    :returns: Name of the fasta file that was written
    """
    outfile = os.path.join(tempDir, "{}.fasta".format(idx))
    logging.debug("Poretools: {} files, output to {}.".format(len(call)-2, outfile))
    f = open(outfile, 'w')
    subprocess.call(call, stdout=f, shell=False)
    f.close()
    return outfile

def callPoretoolsWrapper(args):
    """
    Multiprocessing wrapper for callPoretools

    :param args: List of arguments for callPoretools

    :returns: return value from callPoretools
    """
    return multiprocessWrapper(callPoretools, args)

def multithreadPoretools(poretools, tempDir, force, reads, output, readLength=2000, filesPerCall=1000):
    """
    Runs poretools multithreaded on small chunks of fasta files

    Could be deprecated by find . -name "*.fast5" | parallel "poretools fasta"

    :param poretools: String Path to poretools
    :param tempDir: String path to directory for temporary files
    :param force: Boolean value, force creation or not
    :param output: String Path to final output fasta file
    :param filesPerCall: int number of fast5 files to combine in each call

    TODO: why can't this handle more than eight cores? Arbitrary.
    """
    # check first that we actually want to edit
    if preventOverwrite(output, force):
        return 1

    cwd = os.getcwd()
    os.chdir(reads)
    files = [os.path.join(cwd, reads, f) for f in glob.glob("*.fast5")]
    os.chdir(cwd)
    prefix = [poretools, "fasta", "--type", "fwd", "--min-length", str(readLength)]

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

def buildSortedBam(threads, genome, fastaFile, outPrefix, force, mapq=30):
    """
    Build an indexed sorted bam from a fasta file.

    :param threads: int number of threads
    :param genome: string path to genome fasta file
    :param fastaFile: string path to reads fasta file
    :param force: Boolean value, force creation or not
    :param mapq: int Minimum mapq to output (default: 30, same as bwa mem)

    :returns: string path to sorted bam file
    """
    samFile = "{}.sam".format(outPrefix)
    callSubProcess('{} mem -x ont2d -t {} -T {} {} {}'.format(__exe__['bwa'], threads,
            mapq, genome, fastaFile), force, outputFile=samFile, newFile=samFile)

    sortedBamFile = "{}.sorted.bam".format(outPrefix)
    callSubProcess('samtools sort -o {} -O bam -@ {} -T nanomod {}'.format(sortedBamFile, threads, samFile), force, newFile = sortedBamFile)
    os.remove(samFile)

    callSubProcess('{} index {}'.format(__exe__['samtools'], sortedBamFile),
            force, newFile="{}.bai".format(sortedBamFile))

    return sortedBamFile

def bamReadCount(bamfile):
    """
    Return a tuple of the number of mapped reads in a bam file

    From https://www.biostars.org/p/1890/

    :param bamfile: string Path to indexed bam file

    :returns: int Number of mapped reads
    """
    p = subprocess.Popen(['samtools', 'idxstats', bamfile], stdout=subprocess.PIPE)
    mapped = 0
    for line in p.stdout:
        rname, rlen, nm, nu = line.rstrip().split()
        mapped += int(nm)
    return mapped

def buildEventalign(options, reads, outPrefix):
    """
    Build eventalign file according to nanopolish's pipeline.

    :param options: Namespace object from argparse

    :returns fastaFile: String Filename of fasta corresponding to reads
    :returns eventalignFile: String Filename of eventalign tsv
    :returns readProp: float Proportion of original fast5 files contained in eventalign
    """

    # check how many files we have to begin with
    fast5Count = 0
    for root, dirnames, filenames in os.walk(reads):
        for filename in fnmatch.filter(filenames, '*.fast5'):
            fast5Count += 1
    logging.debug("Found {} fast5 files.".format(fast5Count))

    fastaFile = '{}.fasta'.format(outPrefix)
    # build fasta using poretools so we have index to fast5 files
    if True:
        try:
            # GNU parallel
            callSubProcess(('find {} -name "*.fast5" | parallel -j16 -X {} fasta '
                    '--type fwd --min-length {} $(pwd)/{}').format(reads,
                    __exe__['poretools'], options.readLength, "{}"), options.force, newFile=fastaFile,
                    outputFile=fastaFile)
        except Exception:
            # something went wrong - are we missing paralle?
            # use our custom multiprocessing script
            multithreadPoretools(__exe__['poretools'], options.tempDir, options.force, reads, fastaFile, options.readLength)
    else:
        # single threaded poretools - slow
        callSubProcess(('{} fasta --type fwd {}').format(__exe__['poretools'],
                reads), options.force, newFile=fastaFile, outputFile=fastaFile)

    # index genome using bwa
    callSubProcess('{} index {}'.format(__exe__['bwa'], options.genome),
            options.force, newFile="{}.bwt".format(options.genome))

    # build sorted bam file using bwa mem
    sortedBamFile = buildSortedBam(options.threads, options.genome, fastaFile, outPrefix, options.force, mapq=options.mappingQuality)
    mappedCount = bamReadCount(sortedBamFile)
    logging.debug("Mapped {} reads.".format(mappedCount))
    readProp = float(mappedCount) / fast5Count

    # run nanopolish eventalign
    eventalignFile = "{}.eventalign".format(outPrefix)
    callSubProcess(('{} eventalign -t {} --print-read-names -r {} -b {} -g {}'
            ' --models {}').format(__exe__['nanopolish'], options.threads,
            fastaFile, sortedBamFile, options.genome, options.nanopolishModels),
            options.force, newFile=eventalignFile, outputFile=eventalignFile)

    return fastaFile, eventalignFile, readProp
