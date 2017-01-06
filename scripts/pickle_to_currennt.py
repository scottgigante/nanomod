################################################################################
#                                                                              #
# pickle_to_current.py: reverse engineer a CURRENNT network file from the      #
# .npy network as saved by Nanonet's currennt_to_pickle.py.                    #
#                                                                              #
# This file is part of Nanomod. Say something about a GNU licence?             #
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

def get_parser():
    parser = argparse.ArgumentParser(
        description='Convert currennt json network file into pickle. Makes assumptions about meta data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('input', help='File containing current network')
    parser.add_argument('output', help='Output pickle file')
    return parser

def parse_layer_input(in_network):
	size = in_network.meta['n_features']
	layer = { "name" : "input",
			  "type" : "input",
			  "size" : size }
	return layer

def parse_layer_feedforward(layer, layer_type, idx):
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
	size = layer.out_size
	l = get_layer_dict(layer_type, idx, size)
	weights = dict()
	weights['input'] = layer.W.transpose().flatten()
	weights['bias'] = layer.b
	weights['internal'] = np.ndarray(shape=(0,0))
	return l, weights


def parse_layer_multiclass(in_network):
	size = len(in_network.meta['kmers'])
	layer = { "name": "postoutput",
			  "type": "multiclass_classification",
			  "size": size }
	return layer


def parse_layer_blstm(layer, layer_type, idx):
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
	weights['input'] = wgts_input.reshape(4, wgts_input.shape[0]/4, 
			wgts_input.shape[1], wgts_input.shape[2], 
			wgts_input.shape[3]).transpose(1,2,4,3,0).reshape(-1)
	weights['bias'] = wgts_bias.reshape(-1)
	weights['internal'] = np.append(
			wgts_internalMat.transpose(0,1,3,2).reshape(-1), 
			wgts_internalPeep.reshape(-1))
	return l, weights

def parse_layer_lstm(size, weights):
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
	# retrieve weights from nn
	iW = layer.iW.reshape((-1,1,layer.size))
	lW = layer.lW.reshape(layer.size,4,layer.size).transpose(1,0,2)
	b = layer.b.reshape((4,-1))
	p = layer.p
	return iW, lW, b, p

LAYER_DICT = {'input' : parse_layer_input,
              'blstm' : parse_layer_blstm,
              'feedforward_tanh' : parse_layer_feedforward_tanh,
              'feedforward_logistic' : parse_layer_feedforward_sigmoid,
              'feedforward_identity' : parse_layer_feedforward_linear,
              'lstm' : parse_layer_lstm,
              'blstm' : parse_layer_blstm,
              'softmax' : parse_layer_softmax,
              'multiclass_classification' : parse_layer_multiclass}


def parse_layer(layer, idx):
	layer_type = get_layer_type(layer)
	if not layer_type in LAYER_DICT:
		sys.stderr.write('Unsupported layer type {}.\n'.format(layer_type))
		exit(1)
	l, w = LAYER_DICT[layer_type](layer, layer_type, idx)
	
	#de-numpythonise
	for key in w:
		w[key] = w[key].tolist()
	return l, w

LAYER_NAME_DICT = {'BiRNN' : 'blstm',
				   'FeedForward' : 'feedforward_tanh',
				   'SoftMax' : 'softmax' }

def get_layer_type(layer):
	return LAYER_NAME_DICT[type(layer).__name__]

def get_layer_name(layer_type, idx):
	return "{}_{}".format(layer_type, idx)

def get_layer_dict(layer_type, idx, size):
	return { "name" : "{}_{}".format(layer_type, idx),
			 "type" : layer_type,
			 "size" : size,
			 "bias" : 1.0 } # why bias 1?

def numpy_to_network(in_network):
    """Transform a numpy representation of a network into a json
    representation.
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


if __name__ == '__main__':
    args = get_parser().parse_args() 

    try:
        in_network = np.load(args.input).item()
    except:
        sys.stderr.write('Failed to read from {}.\n'.format(args.input))
        exit(1)

    network = numpy_to_network(in_network)
    with open(args.output, 'w') as out_network:
    	json.dump(network, fp=out_network, indent=4)
