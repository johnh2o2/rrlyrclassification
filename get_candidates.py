print "Loading libraries..."
import os, sys, glob, argparse
print "   settings"
from settings import *
if RUNNING_ON_DELLA:
	import matplotlib as mpl
	mpl.use('Agg')
print "   matplotlib"
import matplotlib.pyplot as plt
print "   numpy"
import numpy as np
print "   interp1s"
from scipy.interpolate import interp1d
print "   fs from utils.feature_selection"
import utils.feature_selection as fs
print "   rhlc from utils.readhatlc"
import utils.readhatlc as rhlc
print "   miscutils"
from utils.miscutils import *
print "   featureutils"
from utils.featureutils import *
print "   the rest!"
import cPickle as pickle
from mpi4py import MPI


# Get rank and size of mpi process
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
ROOT = (rank == 0)
logprint(" get_candidates: MPI: rank = %d size = %d"%(rank, size))

# INPUT ==============================================================
# --------------------------------------------------------------------

# Define command line options
parser = argparse.ArgumentParser(description="(1) processes listed HATNet lightcurves; (2) score them with classifier; (3) pick the best candidates")

parser.add_argument("hatids", nargs='*', help="One or more hatids to process")
parser.add_argument("--dont-save", action='store_false', dest='save', help="Don't save anything from this run")
parser.add_argument("--fields", nargs='*',help="Field(s) to process")
parser.add_argument("--list",nargs='*', help="Path to pickled list of hatids")
parser.add_argument("--rm-orig", action='store_true', help='Remove original .tfalc lightcurve')
parser.add_argument("--force-redo", action='store_true', help='Not exactly sure what this does...')
parser.add_argument("--iteration", type=int, default=get_iteration_number(), help='Not exactly sure what this does...')

# Parse command line arguments
args = parser.parse_args()
iteration = args.iteration
def safe_save(obj, fname):
	if args.save: pickle.dump(obj, open(fname, 'wb'))

# build list of hatids to process
hatids_to_process = []

# Add fields
if not args.fields is None:
	for f in args.fields:
		fname = get_hatids_in_field_fname(f)
		list_of_ids = pickle.load(open(fname, 'rb'))
		logprint( " process: Found %d hatids in field %s"%(len(list_of_ids), f))
		hatids_to_process.extend(list_of_ids)

# Add lists
if not args.list is None:
	for l in args.list:
		if not os.path.exists(l):
			raise Exception( " cannot find %s"%(l))
		
		list_of_ids = pickle.load(open(l, 'rb'))
		logprint(" get_candidates: Found %d hatids in list %s"%(len(list_of_ids), l))
		hatids_to_process.extend(list_of_ids)

# Add specific hatids
if not args.hatids is None:
	logprint(" process: You gave %d hatids on the command line"%(len(args.hatids)))
	for hatid in args.hatids:
		hatids_to_process.append(hatid)

# ====================================================================
# ------------------------------------------------------------- /INPUT

# ORGANIZE ===========================================================
# --------------------------------------------------------------------

# Set the global hatid_field_list variable
hatid_field_list = set_hatid_field_list(all_fields)

# Set the global `fields_to_analyze` list.
set_fields_to_analyze(all_fields)

# Remove any repeats
hatids_to_process = np.unique(hatids_to_process)

# Load the hatid <-> field mapping
def split_hatids_into_fields(hatids):
	'''
	splits a list of hatids into a dictionary of lists
	
	returns { field1 = [ hatids in field 1 ], ... }
	'''
	nofield_hatids = [ hatid for hatid in hatids if not hatid in hatid_field_list ]
	nofield_hatids.extend([ hatid for hatid in hatids if hatid in hatid_field_list and hatid_field_list[hatid] is None ])

	fields = np.unique( [ hatid_field_list[hatid] for hatid in hatids if not hatid in nofield_hatids ] )

	return { f : [ hatid for hatid in hatids if hatid_field_list[hatid] == f ] for f in fields if not f is None }, nofield_hatids

logprint(" get_candidates: sorting hatids by field")

# Split the hatids into their respective fields to save time
field_split, nofield_hatids = split_hatids_into_fields(hatids_to_process)

# Save the hatids with no known field
if len(nofield_hatids) > 0:
	logprint(" get_candidates: %d / %d hatids have no field"%(len(nofield_hatids), len(hatids_to_process)))
	i = 0
	nofield_fname_gen = lambda i : "nofield_hatids_%d.pkl"%(i)
	while os.path.exists(nofield_fname_gen(i)): i += 1
	nofield_fname = nofield_fname_gen(i)
	logprint(" get_candidates: Saving list of hatids with no known field to %s"%(nofield_fname))
	safe_save(nofield_hatids, nofield_fname)

# Load keylist and twomass information for fields
keylists, twomass = {}, {}
if ROOT:
	BAD_IDS = []
	BAD_IDS.extend(nofield_hatids)

	logprint(" get_candidates: Loading keylists.")
	keylists_ = msl.master(field_split.keys())
	logprint(" get_candidates: Loading 2mass information")
	twomass_ = msl.master(field_split.keys())

	# Merge the lists of dicts
	for kl in keylists_:
		field = kl.keys()[0]
		keylists[field] = kl[field]
	for tm in twomass_:
		field = tm.keys()[0]
		twomass[field] = tm[field]
else:
	msl.slave(lambda field : { field : load_keylist(field) })
	msl.slave(lambda field : { field : twomass_info_file(field)} )


# Broadcast keylist and twomass information
keylists = comm.bcast(keylists, root=0)
twomass = comm.bcast(twomass, root=0)

# Prune out the hatids that don't have keylist or 2mass information.
for field in field_split.keys():

	# hatids in field.
	hatids = field_split[field]

	# Does the field have a keylist?
	if field not in keylists or keylists[field] is None:
		logprint(" get_candidates: keylist for field %s is None; %d hatids are going to be labeled BAD"%(field, len(hatids)))
		BAD_IDS.extend(hatids)
		continue

	# Does the field have 2mass information available?
	if not field in twomass or twomass[field] is None:
		logprint(" get_candidates: 2mass info not available for field %s; %d hatids are going to be labeled BAD"%(field, len(hatids)))
		BAD_IDS.extend(hatids)
		continue

	# Are any hatids missing from the 2mass information?
	nno2mass = 0
	for hatid in hatids:
		if hatid not in BAD_IDS and (hatid not in twomass[field] or twomass[field][hatid] is None):
			nno2mass += 1
			BAD_IDS.append(hatid)

	# Print the number of hatids discarded because of missing 2mass information
	if nno2mass > 0: logprint(" get_candidates: 2mass info not available for %d hatids in field %s"%(nno2mass, field))

# ====================================================================
# ---------------------------------------------------------- /ORGANIZE


# PROCESS ============================================================
# --------------------------------------------------------------------

CANDIDATES = []
# Now process the hatids
for field in field_split:
	hatids = field_split[field]
	logprint(" get_candidates: Handling the %d hatids in field %s"%(len(hatids), field))

	# Master/slave workload distribution to make & test features
	if ROOT:

		# Generate lightcurves
		logprint("               :  processing lightcurve!")
		results = msl.master(hatids)
		for ID, status in results:
			if status is False and not ID in bad_ids:
				BAD_IDS.append(ID)

		# Generate features
		logprint("               :  generating features!")
		results = msl.master(hatids)
		for ID, status in results:
			if status is False and not ID in bad_ids:
				BAD_IDS.append(ID)

		# Classify ids
		logprint("               :  scoring hatids!")
		results = msl.master(hatids)
		for ID, scores, status in results:
			if status is None and not ID in bad_ids:
				BAD_IDS.append(ID)
			elif status:
				CANDIDATES.append(ID)
			safe_save(scores, get_scores_fname(ID,iteration))
		
		
	else:
		msl.slave(lambda hatid : (hatid, load_full_tfalc_from_scratch(hatid, field=None, 
				keylist_dat=keylists[field], twomass_dat=twomass[field][hatid], save_full_lc=True, 
				min_observations=5, delete_raw_lc=args.rm_orig, force_redo = args.force_redo)))
		msl.slave(lambda hatid : (hatid, generate_features(hatid)))
		msl.slave(lambda hatid : (hatid, ) + test_hatid(hatid, model_prefix, min_score, min_frac_above_min_score, iteration, N=nmc))

	logprint("               : done.")

# ====================================================================
# ----------------------------------------------------------- /PROCESS

# Save results.
safe_save(BAD_IDS, get_bad_ids_fname(iteration))
safe_save(CANDIDATES, get_candidate_fname(iteration))

logprint("               : %d / %d ids were 'bad'"%(len(BAD_IDS), len(hatids_to_process)))
logprint("               : %d / %d ids were 'candidates'"%(len(CANDIDATES), len(hatids_to_process)))
