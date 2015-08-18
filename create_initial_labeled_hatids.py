import os, sys
from utils.miscutils import *
from utils.featureutils import *
from settings import *
import cPickle as pickle
from mpi4py import MPI
from contextlib import closing
import masterslave as msl

# MPI parameters
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
ROOT = (rank == 0)

# Make the directory for this model if it doesn't already exist
if not os.path.exists(model_output_dir): os.makedirs(model_output_dir)

# Check if the labeled HAT ids file already exists. If so, grumble and exit
if os.path.exists(get_labeled_hatids_fname()):
	print "%s already exists. You should delete %s if you're sure you want to overwrite it."%(get_labeled_hatids_fname(), get_labeled_hatids_fname())
	sys.exit()

logprint(" create_initial_labeled_hatids: Getting full list of gcvs cross matches")
# Get the full list of GCVS cross matches
gcvs_crossmatches 					= np.loadtxt(full_gcvs_match_cat_fname , dtype=match_cat_dt)
gcvs_ids                            = gcvs_crossmatches['id'].tolist()
categories 							= GetCategoriesForEachHATID(gcvs_crossmatches)
good_gcvs_hatids_fname 				= "%s/good_gcvs_hatids.list"%(parent_dir)


if overwrite or not os.path.exists(good_gcvs_hatids_fname):
	# Make list of "good" ID's if one doesn't already exist (or if we're redoing it)
	def is_good(hatid):
		logprint("                              : %s"%(hatid), all_nodes=True)
		if hatid in bad_ids: return False
		if categories[hatid] not in vartypes_to_classify: return False
		return True

	if size == 1:
		# Serial version -- prune out ID's that aren't in bad_ids or are an irrelevant variable star type.
		logprint("                              : Now pruning out bad/irrelevant GCVS sources")
		hatids = []
		for hatid in gcvs_ids:
			if is_good(hatid): hatids.append(hatid)

		logprint("      (saving results to %s)"%(good_gcvs_hatids_fname))
		pickle.dump(hatids, open(good_gcvs_hatids_fname, 'wb'))	

	elif ROOT:
		# Master/slave
		logprint("                              : Now pruning out bad/irrelevant GCVS sources")
		results = msl.master(gcvs_ids)
		hatids = []
		logprint("      (done -- now just deleting things)")
		for hatid, r in results:
			if r: hatids.append(hatid)
		logprint("      (saving results to %s)"%(good_gcvs_hatids_fname))
		pickle.dump(hatids, open(good_gcvs_hatids_fname, 'wb'))	

	else:
		msl.slave(lambda hatid : (hatid, is_good(hatid)))

else:
	# Otherwise load the saved list.
	logprint("                              : Found hatids fname! Reading it.")
	hatids = pickle.load(open(good_gcvs_hatids_fname, 'rb'))
	if not isinstance(hatids, list): hatids = hatids.tolist()
	logprint("                              : done.")

def prune_and_save(results):
	# Add to "BAD_ID's" if there was a problem with the ID
	BAD_IDS = []
	logprint("                              : Pruning out bad ID's")
	for ID, status in results:
		if not status is True: BAD_IDS.append(ID)
	logprint("                              : Writing labeled ID's")
	# Write initial labeled ids!
	if os.path.exists(get_labeled_hatids_fname()): raise Exception("File %s already exists; delete first and then try again."%(fname))

	f = open(get_labeled_hatids_fname(), 'w')
	for ID in hatids:
		if ID in BAD_IDS: continue
		f.write("%-20s%-10i%-20s\n"%(ID, 0, categories[ID]))
	f.close()
	logprint("                               : Done.")

if size == 1:
	
	#logprint("                               : opening ssh connection", all_nodes = True)
	
	# Download full tfalc lightcurves
	#ssh, sftp = open_ssh_connection()
	#logprint("                               : Downloading lcs.")
	dl = lambda hatid : load_full_tfalc_from_scratch(hatid)#, ssh=ssh, sftp=sftp)
	for hatid in hatids: dl(hatid)

	logprint("                              : Now generating features!")
	# Generate features for hatids
	results = []
	for hatid in hatids: results.append((hatid, generate_features(hatid)))

	prune_and_save(results)
	#close_ssh_connection(ssh, sftp)

elif ROOT:
	logprint("                              : Downloading lcs.")
	# Download full tfalc lightcurves
	msl.master(hatids)

	logprint("                              : Now generating features!")
	# Generate features for hatids
	results = msl.master(hatids)
	prune_and_save(results)

else:
	# Download full tfalc lightcurves
	#logprint("                                     : opening ssh connection", all_nodes = True)
	#ssh, sftp = open_ssh_connection()	
	logprint("                                     : Done -- starting slave operation", all_nodes=True)
	msl.slave(lambda hatid : load_full_tfalc_from_scratch(hatid, save_full_lc=True))
	#close_ssh_connection(ssh, sftp)

	# Generate features!
	msl.slave(lambda hatid : (hatid, generate_features(hatid)))

print "Done!"


