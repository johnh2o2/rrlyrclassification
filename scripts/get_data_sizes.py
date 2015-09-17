import os, sys
import numpy as np


def find_size_of_tfalcs(dirname):
	dirsize=0.0
	if os.path.exists(dirname):
		for f in os.listdir(dirname):
			if not "tfalc" in f: continue
			full_fname = os.path.join(dirname, f)
			dirsize+=os.path.getsize(full_fname) * 1E-9
	print "%-75s: %.3f GB"%(dirname, dirsize)
	return dirsize

if __name__ == '__main__':
	import cPickle as pickle
	field_list = pickle.load(open("field_info2.pkl", 'rb'))
	fields =  [ field for field in field_list ]
	#fields = [ fields[i] for i in range(3) ]

	dir_sizes = {}

	for i, field in enumerate(fields):
		if i < 137: continue
		print "finding size of field %s (%d of %d)"%(field, i+1, len(fields))
		dir_sizes[field] = find_size_of_tfalcs(field_list[field])
	all_dir_sizes = [ dir_sizes[field] for field in dir_sizes ]
	pickle.dump(dir_sizes, open("/home/jhoffman/fieldsizes.log", 'wb'))
	print "===================="
	print "Mean:   %.3e GB"%(np.mean(all_dir_sizes))
	print "Std:    %.3e GB"%(np.std(all_dir_sizes))
	print "Total:  %.3e GB"%(sum(all_dir_sizes))
	print "        ------------"

