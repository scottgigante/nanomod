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
# select_data_fraction.sh: use shell scripts to efficiently select a random    #
# number of lines (as well as the header) from a text file. This is used to    #
# randomly select a fraction of the labelled fast5 files to feed into nanonet  #
# for training.                                                                #
#                                                                              #
# @args $1 The fraction, as a float, of reads to be conserved                  #
# @args $2 Name of the training input file                                     #
# @args $3 Name of the validation input file                                   #
# @return $2.small The reduced training output file                            #
# @return $3.small The reduced validation output file                          #
#                                                                              #
# TODO: Allow one or both input files to be empty.                             #
# TODO: Allow a random seed / just take head instead of shuf                   #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

# @
FRAC=$1
TRAIN_IN=$2
VAL_IN=$3

TRAIN_OUT="$TRAIN_IN.small"
VAL_OUT="$VAL_IN.small"

# count lines
TRAIN_N=$(cat $TRAIN_IN | wc -l)
VAL_N=$(cat $VAL_IN | wc -l)

# calculate fraction
TRAIN_FRAC=$(echo "($TRAIN_N * $FRAC)/1" | bc)
VAL_FRAC=$(echo "($VAL_N * $FRAC)/1" | bc)

# headers
head -n 1 $TRAIN_IN > $TRAIN_OUT
head -n 1 $VAL_IN > $VAL_OUT

# random selection
tail -n +2 $TRAIN_IN | shuf -n $TRAIN_FRAC - >> $TRAIN_OUT
tail -n +2 $VAL_IN | shuf -n $VAL_FRAC - >> $VAL_OUT
