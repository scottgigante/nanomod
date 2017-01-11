import sys

from ..test.test_pickle_to_currennt import *

def run(argv):
	test_numpy_conversion()
	test_json_conversion()

if __name__ == "__main__":
	run(sys.argv)