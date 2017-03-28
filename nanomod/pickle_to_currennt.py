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
# pickle_to_current.py: reverse engineer a CURRENNT network file from the      #
# .npy network as saved by Nanonet's currennt_to_pickle.py.                    #
#                                                                              #
# Author: Scott Gigante                                                        #
# Contact: gigante.s@wehi.edu.au                                               #
# Date: 06 Jan 2017                                                            #
#                                                                              #
################################################################################

#!/usr/bin/env python

import argparse
import json
import sys
import numpy as np
import logging

from utils import preventOverwrite

def parse_layer_input(in_network):
    """
    Create input layer from nanonet nn object

    :param in_network: nn network object from nanonet

    :returns layer: dictionary Layer description
    """
    size = in_network.meta['n_features']
    layer = { "name" : "input",
              "type" : "input",
              "size" : size }
    return layer

def parse_layer_feedforward(layer, layer_type, idx):
    """
    Parse an nn FeedForward object.

    :param layer: FeedForward layer object from nanonet

    :returns l: dictionary Layer description
    :returns weights: dictionary Layer weights
    """
    size = layer.out_size
    l = get_layer_dict(layer_type, idx, size)
    weights = dict()
    weights['input'] = layer.W.transpose().flatten()
    weights['bias'] = layer.b
    weights['internal'] = np.ndarray(shape=(0,0))
    return l, weights


def parse_layer_feedforward_tanh(layer, layer_type, idx):
    return parse_layer_feedforward(layer, layer_type, idx)


def parse_layer_feedforward_sigmoid(layer, layer_type, idx):
    return parse_layer_feedforward(layer, layer_type, idx)


def parse_layer_feedforward_linear(size, weights):
    return parse_layer_feedforward(layer, layer_type, idx)


def parse_layer_softmax(layer, layer_type, idx):
    """
    Parse an nn softmax object.

    :param layer: softmax layer object from nanonet

    :returns l: dictionary Layer description
    :returns weights: dictionary Layer weights
    """
    size = layer.out_size
    l = get_layer_dict(layer_type, idx, size)
    weights = dict()
    weights['input'] = layer.W.transpose().flatten()
    weights['bias'] = layer.b
    weights['internal'] = np.ndarray(shape=(0,0))
    return l, weights


def parse_layer_multiclass(in_network):
    """
    Create multiclass output layer from nanonet nn object

    :param in_network: nn network object from nanonet

    :returns layer: dictionary Layer description
    """
    size = len(in_network.meta['kmers'])
    layer = { "name": "postoutput",
              "type": "multiclass_classification",
              "size": size }
    return layer


def parse_layer_blstm(layer, layer_type, idx):
    """
    Parse an nn BLSTM object.

    :param layer: BLSTM layer object from nanonet

    :returns l: dictionary Layer description
    :returns weights: dictionary Layer weights
    """
    # get layer description
    size = layer.layers[0].size * 2
    l = get_layer_dict(layer_type, idx, size)

    # get raw weights from forward and reverse
    iM1, lM1, bM1, pM1 = get_lstm_weights(layer.layers[0])
    iM2, lM2, bM2, pM2 = get_lstm_weights(layer.layers[1].layer)

    # make empty arrays
    weights = dict()
    wgts_input = np.empty([iM1.shape[0], 2, iM1.shape[1], iM1.shape[2]], dtype=np.float32)
    wgts_bias = np.empty([bM1.shape[0], 2, bM1.shape[1]], dtype=np.float32)
    wgts_internalMat = np.empty([lM1.shape[0], 2, lM1.shape[1], lM1.shape[2]], dtype=np.float32)
    wgts_internalPeep = np.empty([pM1.shape[0], 2, pM1.shape[1]], dtype=np.float32)

    # fill arrays
    wgts_input[:, 0, :, :] = iM1
    wgts_bias[:, 0, :] = bM1
    wgts_internalMat[:, 0, :, :] = lM1
    wgts_internalPeep[:, 0, :] = pM1

    wgts_input[:, 1, :, :] = iM2
    wgts_bias[:, 1, :] = bM2
    wgts_internalMat[:, 1, :, :] = lM2
    wgts_internalPeep[:, 1, :] = pM2

    # transform
    weights['input'] = wgts_input.reshape(-1,4,wgts_input.shape[1],
            wgts_input.shape[2],
            wgts_input.shape[3]).transpose(1,2,3,4,0).reshape(-1)
    weights['bias'] = wgts_bias.reshape(-1)
    weights['internal'] = np.append(
            wgts_internalMat.transpose(0,1,3,2).reshape(-1),
            wgts_internalPeep.reshape(-1))
    return l, weights

def parse_layer_lstm(size, weights):
    """
    Parse an nn LSTM object.

    :param layer: LSTM layer object from nanonet

    :returns l: dictionary Layer description
    :returns weights: dictionary Layer weights
    """
    l = get_layer_dict(layer_type, idx, layer.size)

    weights=dict()
    iM, bM, lM, pM = get_lstm_weights(layer)

    # transform
    weights['input'] = iM.transpose(0,2,1).reshape(-1)
    weights['bias'] = bM.reshape(-1)
    weights['internal'] = np.append(lM.transpose(0,2,1).reshape(-1),
            pM.reshape(-1))
    return l, weights

def get_lstm_weights(layer):
    """
    Retrieve LSTM weights from the nn LSTM object.

    :param layer: LSTM layer object from nanonet

    :returns iW: list of floats Input weights
    :returns lW: list of floats Internal weights
    :returns b: list of floats Bias weights
    :returns p: list of floats Peek weights
    """
    iW = layer.iW.reshape((-1,1,layer.size))
    lW = layer.lW.reshape(layer.size,4,layer.size).transpose(1,0,2)
    b = layer.b.reshape((4,-1))
    p = layer.p
    return iW, lW, b, p

__layer_name_dict__ = {'BiRNN' : 'blstm',
                    'FeedForward' : 'feedforward_tanh',
                    'SoftMax' : 'softmax' }

__layer_parse_func_dict__ = {'input' : parse_layer_input,
              'blstm' : parse_layer_blstm,
              'feedforward_tanh' : parse_layer_feedforward_tanh,
              'feedforward_logistic' : parse_layer_feedforward_sigmoid,
              'feedforward_identity' : parse_layer_feedforward_linear,
              'lstm' : parse_layer_lstm,
              'blstm' : parse_layer_blstm,
              'softmax' : parse_layer_softmax,
              'multiclass_classification' : parse_layer_multiclass}

def parse_layer(layer, idx):
    """
    Parse a layer object from nanonet.

    :param layer: Layer object from nanonet
    :param idx: int Layer index

    :returns l: dictionary Layer description
    :returns w: dictionary Layer weights
    """
    layer_type = get_layer_type(layer)
    if not layer_type in __layer_parse_func_dict__:
        logging.error('Unsupported layer type {}.\n'.format(layer_type))
        exit(1)
    l, w = __layer_parse_func_dict__[layer_type](layer, layer_type, idx)

    #de-numpythonise
    for key in w:
        w[key] = w[key].tolist()
    return l, w

def get_layer_type(layer):
    """
    Find the CURRENNT layer type corresponding to a Layer object from nanonet

    :param layer: Layer object from nanonet

    :returns: string Layer type

    >>> get_layer_type(nanonet.nn.BiRNN())
    "blstm"
    """
    return __layer_name_dict__[type(layer).__name__]

def get_layer_name(layer_type, idx):
    """
    Create a name for a given layer.

    :param layer_type: string Type of layer
    :param idx: int Layer index

    :returns: string Layer name

    >>> get_layer_name("blstm", 5)
    "blstm_5"
    """
    return "{}_{}".format(layer_type, idx)

def get_layer_dict(layer_type, idx, size):
    """
    Create a basic layer dictionary.

    :param layer_type: string Type of layer to be created
    :param idx: int Layer index
    :param size: int Number of units in layer

    :returns: dictionary representing layer

    >>> get_layer_dict("blstm", 5, 128)
    { "name" : "blstm_5", "type" : "blstm", "size" : 128, "bias" : 1.0 }

    TODO: why bias 1?
    """
    return { "name" : get_layer_name(layer_type, idx),
             "type" : layer_type,
             "size" : size,
             "bias" : 1.0 }

def numpy_to_network(in_network):
    """
    Transform a numpy representation of a network into a json representation.

    :param in_network: nn object from nanonet to be converted

    :returns: dictionary Coverted network
    """
    layers = list()
    weights = dict()

    # numpy strips input layer
    layers.append(parse_layer_input(in_network))

    i=0
    for layer in in_network.layers:
        l, w = parse_layer(layer, i)
        layers.append(l)
        weights[l['name']] = w
        i += 1

    # numpy strips postoutput layer
    layers.append(parse_layer_multiclass(in_network))

    return { "layers" : layers,
             "weights" : weights,
             "meta" : in_network.meta }

def runConvertPickle(in_filename, out_filename):
    """
    Convert 'pickled' nn object from nanonet to a JSON file.

    This script reverse engineered from currennt_to_pickle.py in nanonet.

    :param in_filename: string Path to pickled network (.npy file)
    :param out_filename: string Path to output network (.json file)

    :returns: int Exit code: 1 if aborted, 0 if run successfully
    """
    try:
        in_network = np.load(in_filename).item()
    except:
        logging.error('Failed to read from {}.\n'.format(in_filename))
        return 1

    network = numpy_to_network(in_network)
    with open(out_filename, 'w') as out_network:
        json.dump(network, fp=out_network, indent=4)
    return 0

def convertPickle(options):
    """
    Run tool from nanomod.

    :param options: Namespace object from argparse

    :returns: int Exit code: 1 if aborted, 0 if run successfully
    """
    # TODO: don't need to pass options to this function
    if preventOverwrite(options.currenntTemplate, options.force) == 1:
        return 1
    logging.info("Converting nanonet network {} to currennt network {}.".format(options.nanonetTemplate, options.currenntTemplate))
    return runConvertPickle(options.nanonetTemplate, options.currenntTemplate)

def parse_args():
    """
    Parse arguments for running from command line.

    :returns: Namespace object from argparse
    """
    parser = argparse.ArgumentParser(
        description='Convert currennt json network file into pickle. Makes assumptions about meta data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('input', help='File containing current network')
    parser.add_argument('output', help='Output pickle file')
    return parser.parse_args()

if __name__ == '__main__':
    """
    This tool can be run directly from the command line.

    Example: python pickle_to_currennt.py in_network.npy out_network.json
    """
    args = parse_args()
    runConvertPickle(args.input, args.output)

