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
# TODO:
# 1. Implement function to *build* a field's keylist / 2mass file on the remote system.
# 2. Add checks to make sure relevant files exist on remote system.
# 3. Add a logger
# 4. Improve load distribution for downloading lightcurves

# Get rank and size of mpi process
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

logprint(" get_candidates: MPI: rank = %d size = %d"%(rank, size))

# What iteration are we on?
iteration = get_iteration_number()

# load a dictionary that contains path information for each field
logprint(" get_candidates: loading field_info")
field_info = pickle.load(open('field_info.dict', 'rb'))
remote_path = field_info[field_to_analyze]['path']
logprint("               : done.")

# specify hostname to connect to and the remote/local paths
#ssh_hostname = 'phn1' # Specified in settings.py
dest_path = LCCACHE
if not os.path.isdir(dest_path): 
	logprint("%s is not a dir. Making it one."%(dest_path),all_nodes=True)
	os.makedirs(dest_path)


# load parameters to setup ssh connection
logprint(" get_candidates: setting up ssh connection...")
config = SSHConfig()
with open(os.path.expanduser('~/.ssh/config')) as config_file:
    config.parse(config_file)
d = config.lookup(ssh_host_name)
logprint("               : done.")

# connect, find all tfalc lightcurves, download them to the local drive
available_hatids = []
logprint(" get_candidates: transferring files...")
with closing(SSHClient()) as ssh:
    ssh.load_system_host_keys() #NOTE: no AutoAddPolicy() 
    ssh.connect(d['hostname'], username=d.get('user'))
    with closing(ssh.open_sftp()) as sftp:
        # cd into remote directory
        sftp.chdir(remote_path)
        # cd to local destination directory
        os.chdir(dest_path)
        # download all files in it to destdir directory

        logprint("               : getting a list of all files in directory %s"%(remote_path))
        raw_all_files = sftp.listdir()
        logprint("               : %d files listed"%(len(raw_all_files)))
        all_files = []
        for filename in raw_all_files:

        	# Skip non-tfalc lightcurves
        	if not 'tfalc' in filename: 
        		#logprint("               : file %s has no tfalc"%(filename))	
        		continue

    		# if more than 1 period in filename, skip.
    		if len(filename.split('.')) > 2: 
    			#logprint("               : file %s has more than one period!"%(filename))
    			continue
        	#logprint("               : [file %s added]"%(filename))

        	all_files.append(filename)

        if not NFILES_MAX is None and len(all_files) > NFILES_MAX:
        	all_files = [ all_files[i] for i in range(NFILES_MAX) ]

        logprint(" get_candidates: len(all_files) = %d"%(len(all_files)))
        nfiles_per_core = len(all_files)/size
        while nfiles_per_core*size < len(all_files): nfiles_per_core+=1

        logprint(" get_candidates: nfiles_per_core = %d"%(nfiles_per_core))
        # ideally this should be a master-slave type workload distribution, but we'll do
        # the naive equal-sized chunk thing for now (to save some coding time)
        # files = comm.scatter()
        files = [ all_files[i] for i in range(rank*nfiles_per_core, min([len(all_files), (rank+1)*nfiles_per_core])) ]

        logprint(" get_candidates: len(files) = %d"%(len(files)))
        for filename in files:
        	
    		hatid, junk = filename.split('.')

    		available_hatids.append(hatid)
        	sftp.get(filename, "%s/%s"%(dest_path,filename))

logprint(" get_candidates: %d == %d?"%(len(np.unique(available_hatids)), len(available_hatids)))

assert(len(np.unique(available_hatids)) == len(available_hatids)) # Make sure there are no repeats!

# Load classifier
logprint(" get_candidates: Loading classifier: %s"%(get_classifier_fname(iteration)))
classifier = BaggedModel()
classifier.load(get_classifier_fname(iteration))
#classifier = pickle.load(open(get_classifier_fname(iteration), 'rb'))

bad_ids = []

logprint(" get_candidates: getting keylist and twomass data...")
if rank == 0:
	# Get keylist + 2mass data.
	keylist = {}
	with closing(SSHClient()) as ssh:
	    ssh.load_system_host_keys() #NOTE: no AutoAddPolicy() 
	    ssh.connect(d['hostname'], username=d.get('user'))
	    with closing(ssh.open_sftp()) as sftp:
	    	logprint("               : %s = (remote) keylist_dir"%(get_keylist_dir(field_to_analyze)))
	    	sftp.chdir(get_keylist_dir(field_to_analyze))

	    	logprint("               : %s = (local) keylist_dir"%(get_local_keylist_dir()))
	        # cd to local destination directory
	        os.chdir(get_local_keylist_dir())

	        klfilesize = sftp.stat("keylist.txt").st_size
	        logprint("               : remote keylist.txt filesize "+`klfilesize`)

	        # I should double check that the keylist.txt file actually exists on the remote system...
	        sftp.get("keylist.txt", get_local_keylist_fname(field_to_analyze))

	        logprint("               : chdir to remote_2massdir;")
	        sftp.chdir(get_remote_2mass_dir())
	        logprint("               : chdir to local_2massdir;")
	        os.chdir(get_local_2mass_dir())


	        sftp.get(get_remote_2mass_fname(field_to_analyze), get_local_2mass_fname(field_to_analyze))

	keylist_data = np.loadtxt(get_local_keylist_fname(field_to_analyze), dtype = keylist_dt)
	for kl in keylist_data:
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

	# TODO: YOU SHOULD BROADCAST KEYLIST DATA INSTEAD OF JUST PICKLING IT.

	pickle.dump(keylist, open("%s-dict.pkl"%(get_local_keylist_fname(field_to_analyze)), 'wb'))
	pickle.dump(twomass_dict, open("%s-dict.pkl"%(get_local_2mass_fname(field_to_analyze)), 'wb'))

#available_hatids = [ "HAT-173-0000034", "HAT-173-0000016", "HAT-173-0000037" ]

####
keylist = pickle.load(open("%s-dict.pkl"%(get_local_keylist_fname(field_to_analyze)), 'rb'))
twomass_dict = pickle.load(open("%s-dict.pkl"%(get_local_2mass_fname(field_to_analyze)), 'rb'))
logprint("               : done.")

def is_candidate(scores):

	nabove = sum([ 1 for s in scores if s > min_score ])
	#print min_score, nabove
	if float(nabove)/float(len(scores)) > min_frac_above_min_score: return True
	return False

candidates = []
logprint(" get_candidates: generating features!")
for ID in available_hatids:
	if ID in bad_ids: continue
	logprint("               : %s"%(ID))
	feat_fname = hat_features_fname(ID, model_prefix)
	if not os.path.exists(feat_fname) or overwrite:
		logprint("                  loading tfalc")
		lc = load_tfalc("%s/%s.tfalc"%(dest_path, ID))
		logprint("                  adding keylist data")
		lc = add_keylist_data(lc, keylist)
		logprint("                  adding twomass data")
		lc = add_2mass(lc, twomass_dict[ID])
		logprint("                  generating features!")
		features = fs.get_features(lc, save_pcov=True, pcov_file=get_pcov_file(ID))
	else:
		logprint("                  loading features")
		features = LoadFeatures(ID)
	# Load/make features
	#features = LoadFeatures(ID)

	# If features is None, this is a bad ID
	if features is None:
		logprint("                  no features...")
		bad_ids.append(ID)
		continue

	# Obtain MC scores
	scores = score_features(features, pcov_file=get_pcov_file(ID), iteration=iteration)

	logprint("                  max score = %.3e!"%(max(scores)))

	# Mark ID if it's a candidate
	if is_candidate(scores): 
		logprint("              ***IS A CANDIDATE!!!!!***")
		candidates.append(ID)
logprint("               : done.")
# gather everything to the root node.
BAD_IDS = np.ravel(comm.gather(bad_ids, root=0))
CANDIDATES = np.ravel(comm.gather(candidates, root=0))
AVAILABLE_HATIDS = np.ravel(comm.gather(available_hatids, root=0))

if rank == 0:
	print "%d available hatids"%(len(AVAILABLE_HATIDS))
	print "%d bad ids"%(len(BAD_IDS))
	print "%d candidates"%(len(CANDIDATES))
	pickle.dump(BAD_IDS, open('%s/bad_ids_iter%04d.pkl'%(model_output_dir, iteration), 'wb'))
	pickle.dump(CANDIDATES, open(get_candidate_fname(iteration), 'wb'))
