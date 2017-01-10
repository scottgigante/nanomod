import json

def reweightArray(net, layer, wgts):
	net['weights'][layer][wgts] = range(len(net['weights'][layer][wgts]))
	return net

def reweightLayer(net, layer):
	net = reweightArray(net, layer, 'input')
	net = reweightArray(net, layer, 'bias')
	net = reweightArray(net, layer, 'internal')
	return net

def reweightNet(net):
	for k in net['weights'].keys():
		net = reweightLayer(net, k)
	return net

def readNet(filename):
	with open(filename, "r") as fh:
		net = json.load(fh)
	return net

def writeNet(net, filename):
	with open(filename, "w") as fh:
		json.dump(net, fp=fh, indent=4)

# main
net = readNet("models/default_template.json")
net = reweightNet(net)
writeNet(net, "models/serial.json")