import os, sys, glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import utils.feature_selection as fs
import utils.readhatlc as rhlc
from utils.miscutils import *
from utils.featureutils import *
from settings import *
import cPickle as pickle
from contextlib import closing
from paramiko import SSHConfig, SSHClient
from mpi4py import MPI

overwrite=True
field_to_analyze = '219'

# Get rank and size of mpi process
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

logprint(" get_candidates: MPI: rank = %d size = %d"%(rank, size))

# What iteration are we on?
iteration = get_iteration_number()

bad_ids = load_bad_ids(iteration)

remote_path = field_info[field_to_analyze]

dest_path = LCCACHE
if not os.path.isdir(dest_path): 
	logprint("%s is not a dir. Making it one."%(dest_path),all_nodes=True)
	os.makedirs(dest_path)

# load parameters to setup ssh connection
logprint(" get_candidates: setting up ssh connection...")
ssh, sftp = open_ssh_connection()
logprint("               : done.")

# Get HAT ID's that are located in this field
field_hatids = get_field_ids(field_to_analyze, sftp)

if not NFILES_MAX is None and len(all_files) > NFILES_MAX:
	hatids = [ field_hatids[i] for i in range(NFILES_MAX) ]

assert(len(np.unique(hatids)) == len(hatids)) # Make sure there are no repeats!

# Download things!
if ROOT:
	results = msl.master(hatids)
else:
	msl.slave(lambda hatid : sftp.get( get_remote_lc_filename(hatid), get_raw_lc_fname(hatid)))

# Load classifier
logprint(" get_candidates: Loading classifier: %s"%(get_classifier_fname(iteration)))
classifier = BaggedModel()
classifier.load(get_classifier_fname(iteration))

logprint(" get_candidates: getting keylist...")


# Obtain keylist for field; this gives you info about filters & stuff for each obs.
if rank == 0: 
	keylist = load_keylist(field_to_analyze)
	
	
close_ssh_connection(ssh, sftp)

logprint("               : broadcasting keylist dictionary")
keylist = comm.bcast(keylist, root=0)
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
	for ID, status in results:
		if status is None and not ID in bad_ids:
			BAD_IDS.append(ID)
		elif status:
			CANDIDATES.append(ID)

	# Print results
	print "%d total hatids"%(len(hatids))
	print "%d bad ids"%(len(BAD_IDS))
	print "%d candidates"%(len(CANDIDATES))

	# Save results
	pickle.dump(BAD_IDS, open(get_bad_ids_fname(iteration), 'wb'))
	pickle.dump(CANDIDATES, open(get_candidate_fname(iteration), 'wb'))

else:
	msl.slave(lambda hatid : (hatid, generate_features(hatid, field=field_to_analyze, keylist=keylist)))
	msl.slave(lambda hatid : (hatid, test_hatid(hatid, model_prefix, min_score, min_frac_above_min_score)))


logprint("               : done.")