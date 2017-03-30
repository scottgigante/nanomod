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
# embed_eventalign.py: use eventalign data to label individual fast5 reads     #
# in order to feed labelled fast5 files into neural network training.          #
#                                                                              #
# TODO: deal with secondary / auxiliary reads better                           #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

from Bio import SeqIO
import os
import pysam
import h5py
import numpy as np
import csv
import logging
import itertools
from numpy.lib.recfunctions import append_fields, drop_fields
from multiprocessing import Pool, cpu_count
from shutil import copyfile

from seq_tools import *
from utils import makeDir, multiprocessWrapper, preventOverwrite
from check_skip_stay_prob import getSkipStayConstraints

__max_read_len__ = 1000000

def getOutfile(fast5, options):
    """
    Get the path of the output file corresponding to a fast5 file.

    :param fast5: Basename of the fast5 file
    :param options: Namespace object from argparse

    :returns: Relative path to output fast5 file

    TODO: remove dependence on options
    """
    return os.path.join(options.outPrefix, fast5)

def writeFast5(options, events, fast5Path, initial_ref_index, last_ref_index,
        chromosome, forward, numSkips, numStays, readLength, kmer, genome):
    """
    Write a labelled fast5 file

    :param options: Namespace object from argparse
    :param events: Structured numpy array including event labels
    :param fast5Path: String relative path to original fast5 file
    :param initial_ref_index: int starting position of read on reference genome
    :param last_ref_index: int ending position of read on reference genome
    :param chromosome: String Name of the chromosome to which fast5 is aligned
    :param forward: Boolean value, True: forward read, False: reverse read
    :param numSkips: int number of skipped bases on genome
    :param numStays: int number of repeated events without moving along genome
    :param readLength: int length of labelled read
    :param kmer: int length of kmer labels
    :param genome: Dictionary containing genome records

    :returns: Boolean, True: high quality read, False: low quality read
    """

    # all events are considered good emissions - bad idea?
    events = append_fields(events, 'good_emission', events["kmer"] != 'X'*kmer)
    # remove move field - nanonet adds it
    events = drop_fields(events, "move")

    # median normalise
    # TODO: we do this twice - make this into a routine.
    # TODO: maybe swap to nanoraw instead of eventalign
    if not options.noNormalize:
        raw = np.concatenate([np.repeat(events["mean"][i], int(events["length"][i] * 4000)) for i in range(events.shape[0])])
        med = np.median(raw) # TODO: pull sample rate from uniqueglobalkey
        median_abs_deviation = np.median(abs(raw - med))
        events["mean"] = (events["mean"] - med) / median_abs_deviation
        events["stdv"] = events["stdv"] / median_abs_deviation

    # remove unlabelled events
    # long stretches of unlabelled events severely decrease accuracy
    events = events[events["kmer"] != 'X'*kmer]

    filename = fast5Path.split('/')[-1]
    newPath = getOutfile(filename, options)
    copyfile(fast5Path, newPath)
    logging.debug("Saving {}".format(filename))
    with h5py.File(newPath,'r+') as fast5:
        # TODO: if the read fails QC, should we just move on?
        attrs = fast5.get(("Analyses/Basecall_1D_000/Summary/basecall_1d"
                "_template")).attrs
        numEvents = float(attrs["called_events"])
        skipProb = attrs["num_skips"]/numEvents
        stayProb = attrs["num_stays"]/numEvents
        stepProb = 1-skipProb-stayProb
        readLength = attrs["sequence_length"]
        passQuality = (skipProb < options.constraints['maxSkips'] and
                stayProb < options.constraints['maxStays'] and
                stepProb > options.constraints['minSteps'] and
                readLength >= options.readLength)

        analysesGroup = fast5.get("Analyses")
        alignToRefGroup = analysesGroup.create_group("AlignToRef")

        # create events
        eventsGroup = alignToRefGroup.create_group(("CurrentSpaceMapped"
                "_template"))
        eventsGroup.create_dataset("Events", data=events)

        # create attrs
        summaryGroup = alignToRefGroup.create_group("Summary")
        attrsGroup = summaryGroup.create_group("current_space_map_template")
        attrs = attrsGroup.attrs
        attrs.create("genome_start", initial_ref_index)
        attrs.create("genome_end", last_ref_index)
        attrs.create("sequence_length", last_ref_index - initial_ref_index)
        attrs.create("ref_name", chromosome)
        attrs.create("direction", "+" if forward else "-")
        attrs.create("num_skips", numSkips)
        attrs.create("num_stays", numStays)
        attrs.create("num_events", numEvents)
        attrs.create("skip_prob", skipProb)
        attrs.create("stay_prob", stayProb)
        attrs.create("step_prob", stepProb)

        # create alignment group
        alignGroup = analysesGroup.create_group("Alignment")
        fastaGroup = alignGroup.create_group("Aligned_template")
        fasta = genome[chromosome]['record'].format("fasta")
        fastaGroup.create_dataset("Fasta",
                data=np.array(fasta, dtype='|S{}'.format(len(fasta))))

        # create attrs
        summaryGroup = alignGroup.create_group("Summary")
        attrsGroup = summaryGroup.create_group("genome_mapping_template")
        attrs = attrsGroup.attrs
        attrs.create("genome", chromosome)

    return passQuality

def processRead(options, idx, fast5Path, genome, modified, kmer, alphabet):
    """
    Step through eventalign data corresponding to a fast5 file and label events

    :param options: Namespace object from argparse
    :param idx: headers index of eventalign file
    :param fast5Path: path to original fast5 file
    :param genome: Dictionary containing genome records
    :param modified: Boolean value, whether or not this fast5 has modified motif
    :param kmer: int length of kmer labels
    :param alphabet: array of base labels in alphabet

    :returns: two-element array containing basename of fast5 file and boolean indicating read quality
    """

    fast5File = fast5Path.split('/')[-1]
    try:
        eventalign = np.load(os.path.join(options.tempDir, fast5File + ".npy"))
    except IOError:
        logging.warning("Failed to load {}.npy".format(os.path.join(options.tempDir, fast5File)))
        # why on earth?! we just created this file!
        return [fast5File, False]
    skip = False

    # open fast5 file
    try:
        fast5 = h5py.File(fast5Path, 'r')
        # r9 for now
        events = fast5.get("Analyses/Basecall_1D_000/BaseCalled_template/Events").value
        events = append_fields(events, 'kmer', ['X'*kmer] * events.shape[0])
        events = append_fields(events, 'seq_pos', [-1] * events.shape[0]) # is there a better option here? keep counting?
    except IOError:
        logging.warning("Failed to read {}".format(fast5Path))
        return [fast5File, False]
    except TypeError:
        logging.warning("Data corrupted in {}".format(fast5Path))
        return [fast5File, False]

    try:
        readLength = fast5.get("Analyses/Basecall_1D_000/Summary/basecall_1d_template").attrs["sequence_length"]
    except AttributeError:
        # attributes not found
        readLength = options.readLength
    initial_event_index = int(eventalign[0][idx['event_index']])
    initial_ref_index = int(eventalign[0][idx['ref_pos']])
    last_ref_index = initial_ref_index-1
    chromosome = eventalign[0][idx['contig']]
    start=False
    numSkips = 0
    numStays = 0

    # check strand direction
    i = 0
    while int(eventalign[0][idx['event_index']]) == int(eventalign[1][idx['event_index']]):
        i += 1
    forward = int(eventalign[i][idx['event_index']]) < int(eventalign[i + 1][idx['event_index']])

    for line in eventalign:

        current_ref_index = int(line[idx['ref_pos']])

        seq_pos_diff = abs(current_ref_index - last_ref_index)
        try:
            seq_pos = last_seq_pos + seq_pos_diff
        except NameError:
            # don't have a last_seq_pos, first event
            seq_pos = 0

        seq = getKmer(genome, chromosome, current_ref_index, kmer, forward)
        if options.rate > 0:
        	seq = randomPermuteSeq(seq, alphabet, options.rate)

        if last_ref_index == current_ref_index:
            numStays += 1
        elif abs(last_ref_index - current_ref_index) > 1:
            # count the size of the skip
            # skipping two bases at once counts as two skips
            numSkips += abs(last_ref_index - current_ref_index) - 1
        last_ref_index = current_ref_index
        last_seq_pos = seq_pos

        events["seq_pos"][int(line[idx['event_index']])] = seq_pos
        events["kmer"][int(line[idx['event_index']])] = seq

    fast5.close()
    try:
        pass_quality = writeFast5(options, events, fast5Path, initial_ref_index,
                last_ref_index, chromosome, forward, numSkips, numStays,
                readLength, kmer, genome)
    except Exception as e:
        logging.warning("Failed to save {}.".format(fast5File))
        print str(e)
        return [fast5File, False]
    return [fast5File, pass_quality]

def processEventalignWorker(options, eventalign, refs, idx, numReads, start, stop):
    """
    Write temporary .npy files to split eventalign tsv file into chunks
    corresponding to one fast5 file each - this allows multiprocessing later on.

    :param options: Namespace object from argparse
    :param eventalign: String filename of eventalign tsv
    :param refs: dictionary linking read names and fast5 file paths
    :param idx: dictionary linking headers with positions in the tsv
    :param start: int Starting line at which to begin processing (0-based)
    :param end: int Ending line at which to end processing (1-based)

    :returns filenames: array_like list of paths to fast5 files which have been processed
    :returns premadeFilenames: list of basenames of fast5 files which had already been processed prior to this call to Nanomod
    """

    with open(eventalign, 'r') as tsv:
        reader = csv.reader(tsv, delimiter="\t")
        reader.next() # headers
        current_read_name = ""
        skip=True
        started=False if start > 0 else True
        tmp=None
        filenames = set()
        premadeFilenames = set()
        n=0

        for i, line in reader.enumerate():

            if i < start:
                continue
            if skip and line[idx['read_name']] == current_read_name:
                continue
            elif line[idx['read_name']] != current_read_name:
                if not started:
                    # just started - could be half way through a read
                    started=True
                    skip=True
                    continue
                if tmp is not None and len(tmp) > 0:
                    try:
                        np.save(filename,tmp)
                        logging.debug("Saving {}.npy".format(filename))
                        n += 1
                        if numReads > 0 and n >= numReads:
                            break
                    except IOError as e:
                        logging.warning("Failed to write {}.npy: {}".format(filename, e))
                if end is not None and i > end:
                    # we're done!
                    break

                # new read, initialise
                tmp = []
                current_read_name = line[idx['read_name']]
                try:
                    fast5Path = refs[current_read_name]
                except KeyError:
                    skip = True
                    logging.warning(current_read_name + " not found.")
                    continue

                fast5Name = fast5Path.split('/')[-1]
                if line[idx['strand']] != 't':
                    # not template, skip
                    skip=True
                    continue

                outfile = getOutfile(fast5Name, options)
                filename = os.path.join(options.tempDir, fast5Name)
                if preventOverwrite(outfile, options.force):
                    skip=True
                    premadeFilenames.add(fast5Name)
                    continue
                elif preventOverwrite(filename + ".npy", options.force):
                    skip=True
                    filenames.add(fast5Path)
                    continue

                filenames.add(fast5Path)
                skip=False

            tmp.append(line)

        # last one gets missed
        # TODO: should we exclude fail reads earlier to save time?
        if tmp is not None and len(tmp) > 0:
            np.save(filename,tmp)
            logging.debug("Saving {}.npy".format(filename))

    return filenames, premadeFilenames

def processEventalign(options, eventalign, refs):
    """
    Write temporary .npy files to split eventalign tsv file into chunks
    corresponding to one fast5 file each - this allows multiprocessing later on.

    :param options: Namespace object from argparse
    :param eventalign: String filename of eventalign tsv
    :param refs: dictionary linking read names and fast5 file paths

    :returns filenames: array_like list of paths to fast5 files which have been processed
    :returns idx: dictionary linking headers with positions in the tsv
    :returns premadeFilenames: list of basenames of fast5 files which had already been processed prior to this call to Nanomod
    :returns kmer: int length of kmer labels
    """
    if not os.path.exists(options.tempDir):
        os.makedirs(options.tempDir)

    with open(eventalign, 'r') as tsv:

        # open file
        reader = csv.reader(tsv, delimiter="\t")
        headers = reader.next()
        idx={}
        idx['contig'] = headers.index("contig")
        idx['ref_pos'] = headers.index("position")
        idx['ref_kmer'] = headers.index("reference_kmer")
        idx['model_kmer'] = headers.index("model_kmer")
        idx['read_name'] = headers.index("read_name")
        idx['strand'] = headers.index("strand")
        idx['event_index'] = headers.index("event_index")
        idx['mean'] = headers.index("event_level_mean")
        idx['stdv'] = headers.index("event_stdv")
        idx['length'] = headers.index("event_length")

        line = reader.next()
        kmer=len(line[idx['ref_kmer']])
        # TODO: can we do this without reopening?

    if options.threads > 1 and options.numReads <= 0:
        nlines = sum(1 for _ in open(eventalign, 'r'))
        pool = Pool(options.threads)
        data = pool.map(processEventalignWorkerWrapper,
                [[options, eventalign, refs, idx, -1, i, i + __max_read_len__] for i in reversed(xrange(0, nlines, __max_read_len__))])
        pool.close()
        pool.join()
        filenames, premadeFilenames = tuple(itertools.chain(*i) for i in itertools.izip(*data))
    else:
        filenames, premadeFilenames = processEventalignWorker(options, eventalign, refs, idx, options.numReads, 0, None)
    return filename, idx, premadeFilenames, kmer

def checkPremade(options, filename):
    """
    Check quality of fast5 files generated in previous runs of nanomod

    :param options: Namespace object from argparse
    :param filename: Name of premade fast5 file

    TODO: remove dependence on options
    """
    try:
        with h5py.File(getOutfile(filename, options), 'r') as fh:
            attrs = fh.get(("Analyses/Basecall_1D_000/Summary/basecall_1d"
                    "_template")).attrs
            numEvents = float(attrs["called_events"])
            skipProb = attrs["num_skips"]/numEvents
            stayProb = attrs["num_stays"]/numEvents
            stepProb = 1-skipProb-stayProb
            readLength = attrs["sequence_length"]
            passQuality = (skipProb < options.constraints['maxSkips'] and
                    stayProb < options.constraints['maxStays'] and
                    stepProb > options.constraints['minSteps'] and
                    readLength >= options.readLength)
    except IOError:
        logging.warning("Failed to read {}".format(filename))
        passQuality = False
    except AttributeError:
        # can't find attributes
        logging.info("Attributes missing on {}".format(filename))
        passQuality = False
    return [filename, passQuality]

def writeTrainfiles(options, trainData, outPrefix):
    """
    Write text files for training and validation consisting of fast5 basenames

    Small text files are used for selective training of a smaller dataset.

    :param options: Namespace object from argparse
    :param trainData: Dataset consisting or two-element arrays of filename and boolean indicating read quality, where reads marked with True are included in the small dataset
    """
    trainFilename = outPrefix + ".train.txt"
    valFilename = outPrefix + ".val.txt"

    with open(trainFilename, 'w') as trainFile, open(valFilename, 'w') as valFile, open(trainFilename + ".small", 'w') as smallTrainFile, open(valFilename + ".small", 'w') as smallValFile:

        # write trainFile / valFile headers
        trainFile.write("#filename\n")
        valFile.write("#filename\n")
        smallTrainFile.write("#filename\n")
        smallValFile.write("#filename\n")

        for filename, pass_quality in trainData:
            if os.path.isfile(getOutfile(filename, options)):
                if np.random.rand() > options.valFraction:
                    logging.debug("Training set{}: {}".format(" - Selected" if pass_quality else "", filename))
                    trainFile.write(filename + "\n")
                    if pass_quality:
                        smallTrainFile.write(filename + "\n")
                else:
                    logging.debug("Validation set{}: {}".format(" - Selected" if pass_quality else "", filename))
                    valFile.write(filename + "\n")
                    if pass_quality:
                        smallValFile.write(filename + "\n")

def processReadWrapper(args):
    """
    Multiprocessing wrapper for processRead

    :param args: List of arguments for processRead

    :returns: return value from processRead
    """
    return multiprocessWrapper(processRead, args)


def processEventalignWorkerWrapper(args):
    """
    Multiprocessing wrapper for processEventalignWorker

    :param args: List of arguments for processEventalignWorker

    :returns: return value from processEventalignWorker
    """
    return multiprocessWrapper(processEventalignWorker, args)

def checkPremadeWrapper(args):
    """
    Multiprocessing wrapper for checkPremade

    :param args: List of arguments for checkPremade

    :returns: return value from checkPremade
    """
    return multiprocessWrapper(checkPremade, args)

def embedEventalign(options, fasta, eventalign, reads, outPrefix, modified):
    """
    Embed labels from an eventalign run into fast5 files for training by nanonettrain

    :param options: Namespace object from argparse
    :param fasta: String filename corresponding to fasta file containing reads (MUST be produced by poretools in order to retrieve fasta names)
    :param eventalign: String filename of eventalign tsv output
    :param reads: String path to fast5 files
    :param outPrefix: String prefix for output files
    :param modified: Boolean value, whether or not reads have modified motif
    """
    output = "{}.train.txt.small".format(outPrefix)
    if preventOverwrite(output, options.force):
        return 1
    makeDir(options.outPrefix)

    logging.info("Loading fast5 names and references...")
    refs = loadRef(fasta)
    genome = loadGenome(options, modified)

    if "random" in options.selectMode:
        # no point populating constrained files, we will overwrite later
        options.constraints = {'maxSkips' : 0,
            'maxStays' : 0,
            'minSteps' : 1 }
    else:
        logging.info("Calculating skip/stay count constraints...")
        options.constraints = getSkipStayConstraints(reads,
                options.dataFraction, options.selectMode, options.readLength)
        logging.debug(str(options.constraints))

    logging.info("Splitting eventalign into separate files...")
    filenames, idx, premadeFilenames, kmer = writeTempFiles(options, eventalign, refs)
    pool = Pool(options.threads)
    alphabet = expandAlphabet(options.sequenceMotif)

    logging.info("Embedding labels into fast5 files...")
    trainData = pool.map(processReadWrapper,
            [[options, idx, i, genome, modified, kmer, alphabet] for i in filenames])

    logging.info("Adding data for premade fast5 files...")
    trainData.extend(pool.map(checkPremadeWrapper,
            [[options, i] for i in premadeFilenames]))
    pool.close()
    pool.join()

    logging.info("Saving datasets for {} total fast5 files...".format(len(trainData)))
    writeTrainfiles(options, trainData, outPrefix)
