import os, sys
import utils.miscutils as mutils
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
download = False

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
	
	logprint("                              : Pruning out bad ID's")
	
	BAD_IDS = [ h for h,r in results if not r ]
	pickle.dump(BAD_IDS, open(get_bad_ids_fname(0), 'wb'))

	logprint("                              : Writing labeled ID's")
	# Write initial labeled ids!
	if os.path.exists(get_labeled_hatids_fname()): raise Exception("File %s already exists; delete first and then try again."%(fname))

	f = open(get_labeled_hatids_fname(), 'w')
	for ID in hatids:
		if ID in BAD_IDS: continue
		if not ID in categories: raise Exception("Cannot find category for %s."%(ID))
		if categories[ID] is None: continue
		f.write("%-20s%-10i%-20s\n"%(ID, 0, categories[ID]))
	f.close()


	logprint("                               : Done.")

add_twomass_info_field('gcvs')
tmi = mutils.twomass_info_for_field['gcvs']
def tmdat(hatid):
	if not hatid in tmi: return None
	return tmi[hatid]
dl = lambda hatid : load_full_tfalc_from_scratch(hatid, twomass_dat=tmdat(hatid), save_full_lc=True, force_redo=False)

if size == 1:
	if download:
		for hatid in hatids: dl(hatid)

	logprint("                              : Now generating features!")

	# Generate features for hatids
	results = []
	for hatid in hatids: results.append((hatid, generate_features(hatid)))

	logprint("                              : %d OK sets of features!"%(len([ h for h,r in results if r ])))

	prune_and_save(results)
	#close_ssh_connection(ssh, sftp)

elif ROOT:
	logprint("                              : Fixing lcs.")
	# Download full tfalc lightcurves
	if download : msl.master(hatids)

	logprint("                              : Now generating features!")
	# Generate features for hatids
	results = msl.master(hatids)

	logprint("                              : %d OK sets of features!"%(len([ h for h,r in results if r ])))

	# Save!
	prune_and_save(results)

else:
	# Load full tfalcs from scratch
	if download : msl.slave(dl)

	# Generate features!
	msl.slave(lambda hatid : (hatid, generate_features(hatid)))

print "Done!"


