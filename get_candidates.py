import os, sys, glob

from settings import *
if RUNNING_ON_DELLA:
	import matplotlib as mpl
	mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import utils.feature_selection as fs
import utils.readhatlc as rhlc
from utils.miscutils import *
from utils.featureutils import *
import cPickle as pickle
from contextlib import closing
from paramiko import SSHConfig, SSHClient
from mpi4py import MPI

#fields_to_analyze = [ '145' ]
overwrite=True

# Get rank and size of mpi process
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
ROOT = (rank == 0)


logprint(" get_candidates: MPI: rank = %d size = %d"%(rank, size))

# What iteration are we on?
iteration = get_iteration_number()
logprint(" get_candidates: Looks like we're on iteration %d"%(iteration))

dest_path = LCCACHE
if not os.path.isdir(dest_path): 
	raise Exception("The dest_path %s does not exist on the system."%(dest_path))


# Get ID's to skip.
bad_ids = load_bad_ids(iteration)

# Get list of labeled hatids (obviously avoid these...)
if not os.path.exists(get_labeled_hatids_fname()): raise Exception(" Cannot find labeled hatids file %s"%(get_labeled_hatids_fname()))
labeled_hatids = np.loadtxt(get_labeled_hatids_fname(), dtype=np.dtype([('hatid', 'S15'), ('iter', np.int_), ('label', 'S15')]))['hatid'].tolist()
#print len(labeled_hatids)
# Get HAT ID's that are located in the fields that we're analyzing (and that aren't already labeled)
hatids = [ hatid for hatid in hatid_field_list if get_field_of(hatid) in fields_to_analyze and not hatid in labeled_hatids ]


#hatids = labeled_hatids 
#print len(hatids)

if not NFILES_MAX is None and len(hatids) > NFILES_MAX:
	logprint(" get_candidates: only using first %d hatids."%(NFILES_MAX))
	hatids = [ hatids[i] for i in range(NFILES_MAX) ]

assert(len(np.unique(hatids)) == len(hatids)) # Make sure there are no repeats!

# Load classifier
logprint(" get_candidates: Loading classifier: %s"%(get_classifier_fname(iteration)))
classifier = pickle.load(open(get_classifier_fname(iteration), 'rb'))

logprint(" get_candidates: getting keylist...")

keylists=None
# Obtain keylists for relevant field(s); this gives you info about filters & stuff for each obs.
if ROOT: 
	logprint(" get_candidates: Loading keylists")
	keylists = {field : load_keylist(field) for field in fields_to_analyze }
	
logprint("               : broadcasting keylist dictionary")
keylists = comm.bcast(keylists, root=0)
logprint("               : done.")

# Master/slave workload distribution to make & test features
if ROOT:
	BAD_IDS, CANDIDATES = [],[]

	# Generate features
	logprint(" get_candidates: generating features!")
	results = msl.master(hatids)
	for ID, status in results:
		if status is False and not ID in bad_ids:
			BAD_IDS.append(ID)
	
	# Classify ids
	logprint(" get_candidates: testing hatids!")
	results = msl.master(hatids)
	for ID, scores, status in results:
		if status is None and not ID in bad_ids:
			BAD_IDS.append(ID)
		elif status:
			CANDIDATES.append(ID)
		pickle.dump(scores, open(get_scores_fname(ID, iteration), 'wb'))
	
	# Print results
	print "%d total hatids"%(len(hatids))
	print "%d bad ids"%(len(BAD_IDS))
	print "%d candidates"%(len(CANDIDATES))

	# Save results
	pickle.dump(BAD_IDS, open(get_bad_ids_fname(iteration), 'wb'))
	pickle.dump(CANDIDATES, open(get_candidate_fname(iteration), 'wb'))
	
else:
	msl.slave(lambda hatid : (hatid, generate_features(hatid, field=hatid_field_list[hatid])))#, keylist=keylists[hatid_field_list[hatid]])))
	msl.slave(lambda hatid : (hatid, ) + test_hatid(hatid, model_prefix, min_score, min_frac_above_min_score, iteration, N=nmc))
	

logprint("               : done.")
