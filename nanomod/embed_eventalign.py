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

def processEventalignWorker(options, eventalign, refs, idx, numReads, start, stop, genome, modified, kmer, alphabet):
    """fast5Path, genome, modified, kmer, alphabet
    Write temporary .npy files to split eventalign tsv file into chunks
    corresponding to one fast5 file each - this allows multiprocessing later on.

    :param options: Namespace object from argparse
    :param eventalign: String filename of eventalign tsv
    :param refs: dictionary linking read names and fast5 file paths
    :param idx: dictionary linking headers with positions in the tsv
    :param start: int Starting line at which to begin processing (0-based)
    :param end: int Ending line at which to end processing (1-based)
    :param genome: Dictionary containing genome records
    :param modified: Boolean value, whether or not this fast5 has modified motif
    :param kmer: int length of kmer labels
    :param alphabet: array of base labels in alphabet

    :returns trainData: array of two-element arrays containing basename of fast5 file and boolean indicating read quality
    :returns premadeFilenames: list of basenames of fast5 files which had already been processed prior to this call to Nanomod
    """

    with open(eventalign, 'r') as tsv:
        reader = csv.reader(tsv, delimiter="\t")
        reader.next() # headers
        current_read_name = ""
        skip=True
        started=False if start > 0 else True
        events=None
        filenames = set()
        trainData = list()
        premadeFilenames = set()
        n=0

        for i, line in enumerate(reader):

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

                #########################
                # save previous read
                #########################
                if events is not None:
                    try:
                        passQuality = writeFast5(options, events, fast5Path,
                            initial_ref_index, last_ref_index, chromosome, forward,
                            numSkips, numStays, readLength, kmer, genome)
                        trainData.append([fast5Name, passQuality])
                        filenames.add(fast5Name)
                        events = None
                        n += 1
                        if options.numReads > 0 and n >= options.numReads:
                            break
                    except Exception as e:
                        logging.warning("Failed to save {}.".format(fast5Name))
                        logging.warning(str(e))
                if stop is not None and i > stop:
                    # we're done!
                    break

                ############################
                # new read, initialise
                ############################
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
                elif fast5Name in filenames or fast5Name in premadeFilenames:
                    skip = True
                    continue

                outfile = getOutfile(fast5Name, options)
                filename = os.path.join(options.tempDir, fast5Name)
                if preventOverwrite(outfile, options.force):
                    skip=True
                    premadeFilenames.add(fast5Name)
                    continue

                # open fast5 file
                try:
                    fast5 = h5py.File(fast5Path, 'r')
                    # r9 for now
                    events = fast5.get("Analyses/Basecall_1D_000/BaseCalled_template/Events").value
                    events = append_fields(events, 'kmer', ['X'*kmer] * events.shape[0])
                    events = append_fields(events, 'seq_pos', [-1] * events.shape[0]) # is there a better option here? keep counting?
                except IOError:
                    logging.warning("Failed to read {}".format(fast5Path))
                    skip=True
                    continue
                except TypeError:
                    logging.warning("Data corrupted in {}".format(fast5Path))
                    skip=True
                    continue

                try:
                    readLength = fast5.get("Analyses/Basecall_1D_000/Summary/basecall_1d_template").attrs["sequence_length"]
                except AttributeError:
                    # attributes not found
                    readLength = options.readLength
                initial_event_index = int(line[idx['event_index']])
                initial_ref_index = int(line[idx['ref_pos']])
                last_ref_index = initial_ref_index-1
                chromosome = line[idx['contig']]
                numSkips = 0
                numStays = 0
                forward = None
                seq_pos = None
                last_seq_pos = None

                filenames.add(fast5Path)
                skip=False

            ########################
            # process line
            ########################

            # check strand direction
            if forward is None:
                if int(line[idx['event_index']]) == initial_ref_index:
                    continue
                else:
                    forward = initial_ref_index < int(line[idx['event_index']])
            current_ref_index = int(line[idx['ref_pos']])

            # generate seq pos
            seq_pos_diff = abs(current_ref_index - last_ref_index)
            try:
                seq_pos = last_seq_pos + seq_pos_diff
            except TypeError:
                # don't have a last_seq_pos, first event
                seq_pos = 0

            # pull seq from genome
            seq = getKmer(genome, chromosome, current_ref_index, kmer, forward)
            if options.rate > 0:
                seq = randomPermuteSeq(seq, alphabet, options.rate)

            # count skips and stays
            if last_ref_index == current_ref_index:
                numStays += 1
            elif abs(last_ref_index - current_ref_index) > 1:
                # count the size of the skip
                # skipping two bases at once counts as two skips
                numSkips += abs(last_ref_index - current_ref_index) - 1
            last_ref_index = current_ref_index
            last_seq_pos = seq_pos

            # add data to events table
            events["seq_pos"][int(line[idx['event_index']])] = seq_pos
            events["kmer"][int(line[idx['event_index']])] = seq

        # last one gets missed
        # TODO: can we exclude fail reads earlier to save time?
        if events is not None:
            try:
                passQuality = writeFast5(options, events, fast5Path,
                    initial_ref_index, last_ref_index, chromosome, forward,
                    numSkips, numStays, readLength, kmer, genome)
                trainData.append([fast5Name, passQuality])
            except Exception as e:
                logging.warning("Failed to save {}.".format(fast5File))
                logging.warning(str(e))

    return trainData, premadeFilenames

def processEventalign(options, eventalign, refs, genome, modified, alphabet):
    """
    Write temporary .npy files to split eventalign tsv file into chunks
    corresponding to one fast5 file each - this allows multiprocessing later on.

    :param options: Namespace object from argparse
    :param eventalign: String filename of eventalign tsv
    :param refs: dictionary linking read names and fast5 file paths
    :param genome: Dictionary containing genome records
    :param modified: Boolean value, whether or not this fast5 has modified motif
    :param alphabet: array of base labels in alphabet

    :returns trainData: array of two-element arrays containing basename of fast5 file and boolean indicating read quality
    :returns premadeFilenames: list of basenames of fast5 files which had already been processed prior to this call to Nanomod
    """

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
        threads = max(min(nlines / __max_read_len__, options.threads), 1)
        pool = Pool(threads)
        splitRange = list(reversed(xrange(0, nlines, nlines / threads)))
        if nlines % threads != 0:
            # drop the overhang - better to have one long thread than an extra short one
            splitRange = splitRange[1:]
        splitRange = [nlines] + splitRange

        data = pool.map(processEventalignWorkerWrapper,
                [[options, eventalign, refs, idx, -1, splitRange[i+1], splitRange[i], genome, modified, kmer, alphabet] for i in range(len(splitRange)-1)])
        pool.close()
        pool.join()
        trainData, premadeFilenames = tuple(itertools.chain(*i) for i in itertools.izip(*data))
    else:
        trainData, premadeFilenames = processEventalignWorker(options, eventalign, refs, idx, options.numReads, 0, None, genome, modified, kmer, alphabet)
    return trainData, premadeFilenames

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
    seenFiles = set()

    with open(trainFilename, 'w') as trainFile, open(valFilename, 'w') as valFile, open(trainFilename + ".small", 'w') as smallTrainFile, open(valFilename + ".small", 'w') as smallValFile:

        # write trainFile / valFile headers
        trainFile.write("#filename\n")
        valFile.write("#filename\n")
        smallTrainFile.write("#filename\n")
        smallValFile.write("#filename\n")

        for filename, pass_quality in trainData:
            if filename in seenFiles:
                continue
            seenFiles.add(filename)
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
    :param eventalign: subprocess.Process Process containing eventalign output
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
    genome = loadGenome(options.genome, options.sequenceMotif, modified)

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
    alphabet = expandAlphabet(options.sequenceMotif)
    trainData, premadeFilenames = processEventalign(options, eventalign, refs, genome, modified, alphabet)

    pool = Pool(options.threads)
    logging.info("Adding data for premade fast5 files...")
    itertools.chain(trainData, pool.map(checkPremadeWrapper,
            [[options, i] for i in premadeFilenames]))
    pool.close()
    pool.join()

    logging.info("Saving datasets...")
    writeTrainfiles(options, trainData, outPrefix)
