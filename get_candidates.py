import matplotlib.pyplot as plt
import readhatlc as rhlc
import numpy as np
import os, sys, glob
from scipy.interpolate import interp1d
import feature_selection as fs
from utils import *
from featureutils import *
from settings import *
import cPickle as pickle
from sklearn.metrics import confusion_matrix
from contextlib import closing
from paramiko import SSHConfig, SSHClient
from mpi4py import MPI

min_score = 0.05
min_frac_above_min_score = 0.3

# TODO:
# 1. Implement function to *build* a field's keylist / 2mass file on the remote system.
# 2. Add checks to make sure relevant files exist on remote system.
# 3. Add a logger
# 4. Improve load distribution for downloading lightcurves

# Get rank and size of mpi process
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# What iteration are we on?
iteration = 1
while not os.path.exists(get_classifier_fname(iteration)): iteration += 1

# load a dictionary that contains path information for each field
field_info = pickle.load(open('field_info.dict', 'rb'))
remote_path = field_info[field_to_analyze]['path']

# specify hostname to connect to and the remote/local paths
hostname = 'phn1'
dest_path = LCCACHE


# load parameters to setup ssh connection
config = SSHConfig()
with open(os.path.expanduser('~/.ssh/config')) as config_file:
    config.parse(config_file)
d = config.lookup(hostname)

# connect, find all tfalc lightcurves, download them to the local drive
available_hatids = []

with closing(SSHClient()) as ssh:
    ssh.load_system_host_keys() #NOTE: no AutoAddPolicy() 
    ssh.connect(d['hostname'], username=d.get('user'))
    with closing(ssh.open_sftp()) as sftp:
        # cd into remote directory
        sftp.chdir(remote_path)
        # cd to local destination directory
        os.chdir(dest_path)
        # download all files in it to destdir directory

        
        all_files = stfp.listdir()

        nfiles_per_core = len(all_files)/size
        while nfiles_per_core*size < len(all_files): nfiles_per_core+=1

        # ideally this should be a master-slave type workload distribution, but we'll do
        # the naive equal-sized chunk thing for now (to save some coding time)
        #files = comm.scatter()
        files = [ all_files[i] for i in range(rank*nfiles_per_core, min(len(all_files), (rank+1)*nfiles_per_core)) ]

        for filename in files:
        	# Skip non-tfalc lightcurves
        	if not 'tfalc' in filename: continue

    		# if more than 1 period in filename, skip.
    		if len(filename.split('.')) > 2: continue
    		hatid, junk = filename.split('.')

    		available_hatids.append(hatid)
        	sftp.get(filename, filename)


assert(len(np.unique(available_hatids)) == available_hatids) # Make sure there are no repeats!

# Load classifier
classifier = pickle.load(open(get_classifer_fname(iteration), 'rb'))

bad_ids = []

if rank == 0:
	# Get keylist + 2mass data.
	keylist = {}
	with closing(SSHClient()) as ssh:
	    ssh.load_system_host_keys() #NOTE: no AutoAddPolicy() 
	    ssh.connect(d['hostname'], username=d.get('user'))
	    with closing(ssh.open_sftp()) as sftp:
	    	sftp.chdir(get_keylist_dir(field_to_analyze))
	        # cd to local destination directory
	        os.chdir(get_local_keylist_dir())

	        # I should double check that the keylist.txt file actually exists on the remote system...
	        stfp.get("keylist.txt", get_local_keylist_fname(field_to_analyze))

	        stfp.chdir(get_remote_2mass_dir())
	        os.chdir(get_local_2mass_dir())

	        stfp.get(get_remote_2mass_fname(field_to_analyze), get_local_2mass_fname(field_to_analyze))

	keylist_data = np.loadtxt(get_local_keylist_fname(field_to_analyze), dtype = keylist_dt)
	for kl in klist_data:
		keylist[kl['TSTF']] = { }
		for k in keylist_dt.names:
			if k == 'TSTF': continue
			keylist[kl['TSTF']][k] = kl[k]

	twomass_dict = {}
	twomass_data = np.loadtxt(get_local_2mass_fname(field_to_analyze), dtype=twomass_dt)
	for tm in twomass_data:
		twomass_dict[tm['hatid']] = {}
		for c in twomass_dt.names:
			if c == 'hatid': continue
			twomass_dict[tm['hatid']][c] = tm[c]

	# BROADCAST KEYLIST DATA INSTEAD OF JUST PICKLING IT.

	pickle.dump(keylist, "%s-dict.pkl"%(get_local_keylist_fname(field_to_analyze)))
	pickle.dump(twomass_dict, "%s-dict.pkl"%(get_local_2mass_fname(field_to_analyze)))

keylist = pickle.load(open("%s-dict.pkl"%(get_local_keylist_fname(field_to_analyze)), 'rb'))
twomass_dict = pickle.load(open("%s-dict.pkl"%(get_local_2mass_fname(field_to_analyze)), 'rb'))

def is_candidate(scores):
	nabove = sum([ s for s in scores if s > min_score ])
	if float(nabove)/float(len(scores)) > min_frac_above_min_score: return True
	return False

candidates = []
for ID in available_hatids:
	feat_fname = hat_features_fname(ID, model_name)
	if not os.path.exists(feat_fname) or overwrite:
		lc = load_tfalc("%s/%s.tfalc"%(dest_path, ID))
		lc = add_keylist_data(lc, keylist)
		lc = add_2mass(lc, twomass_dict[ID])
		features = fs.get_features(lc, save_pcov=True, pcov_file=get_pcov_file(ID))
	else:
		features = LoadFeatures(ID)
	# Load/make features
	#features = LoadFeatures(ID)

	# If features is None, this is a bad ID
	if features is None:
		bad_ids.append(ID)
		continue

	# Obtain MC scores
	scores = score_features(features, pcov_file=get_pcov_file(ID), iteration=iteration)

	# Mark ID if it's a candidate
	if is_candidate(scores): candidates.append(ID)

# gather everything to the root node.
BAD_IDS = np.ravel(comm.gather(bad_ids), root=0)
CANDIDATES = np.ravel(comm.gather(candidates), root=0)
AVAILABLE_HATIDS = np.ravel(comm.gather(available_hatids, root=0))

if rank == 0:
	print "%d available hatids"%(len(AVAILABLE_HATIDS))
	print "%d bad ids"%(len(BAD_IDS))
	print "%d candidates"%(len(CANDIDATES))
	pickle.dump(BAD_IDS, open('%s/bad_ids.pkl', 'wb'))
	pickle.dump(CANDIDATES, open('%s/candidates.pkl', 'wb'))
