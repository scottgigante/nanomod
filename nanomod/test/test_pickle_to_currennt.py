import numpy as np
import json
from nanonet.currennt_to_pickle import network_to_numpy

from ..scripts.pickle_to_currennt import numpy_to_network

NANONET_NN = "nanomod/test/data/r9_template.npy"
SERIAL_JSON = "nanomod/test/data/serial.json"

def assert_network_equivalence(n1, n2):
	assert(len(n1.layers) == len(n2.layers))
	for i in range(len(n1.layers)):
		assert_layer_equivalence(n1.layers[i], n2.layers[i])
	return 0

def assert_layer_equivalence(l1, l2):
	assert(type(l1) == type(l2))
	try:
		# BiRNN
		assert(len(l1.layers) == len(l2.layers))
		for i in range(len(l1.layers)):
			assert_layer_equivalence(l1.layers[i], l2.layers[i])
	except AttributeError:
		try:
			# Reverse
			assert_layer_equivalence(l1.layer, l2.layer)
		except AttributeError:
			# Base layer
			assert((l1.b == l2.b).all())
			try:
				# Feedforward
				assert((l1.W == l2.W).all())
			except AttributeError:
				# LSTM
				assert(l1.size == l2.size)
				assert((l1.iW == l2.iW).all())
				assert((l1.lW == l2.lW).all())
				assert((l1.p == l2.p).all())
	return 0

def test_numpy_conversion():
	in_network = np.load(NANONET_NN).item()
	converted_network = numpy_to_network(in_network)
	out_network = network_to_numpy(converted_network)

	return assert_network_equivalence(in_network, out_network)

def assert_serial_network(net):
	for layer in net['weights']:
		for weights in net['weights'][layer]:
			previous = -1
			for w in net['weights'][layer][weights]:
				assert(w > previous)
				previous = w
	return 0

def test_json_conversion():
	with open(SERIAL_JSON, 'r') as fh:
		in_network = json.load(fh)
	converted_network = network_to_numpy(in_network)
	out_network = numpy_to_network(converted_network)

	return assert_serial_network(out_network)

