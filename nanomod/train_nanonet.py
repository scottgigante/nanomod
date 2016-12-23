from utils import callSubProcess
from . import __exe__

def trainNanonet(options):
	
	callSubProcess(("{0} --train {1} --train_list {1}.train.txt --val {1} "
				   "--val_list {1}.val.txt --output {1} --kmer_length {2}"
				   "--parallel_sequences {3} --workspace {4} --cache_path {4}"
				   "--model {5} --cuda").format(__exe__['nanonettrain'], 
				   options.outPrefix, options.kmer, options.parallelSequences, 
				   options.tempDir, options.currenntTemplate)
	# return ???