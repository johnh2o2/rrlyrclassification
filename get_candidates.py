import os, sys, glob, argparse
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
from mpi4py import MPI
logprint(" get_candidates: loaded libraries", all_nodes=True)

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
parser.add_argument('--parallelize-by-hatid', action='store_false', dest='parallelize_fields', help="Pass this flag if num_hatids >>> num_fields")

# Parse command line arguments
args = parser.parse_args()
iteration = args.iteration
parallelize_fields = args.parallelize_fields

def safe_save(obj, fname):
	if args.save: pickle.dump(obj, open(fname, 'wb'))

# build list of hatids to process
hatids_to_process = []

# Add fields
if not args.fields is None:
	for f in args.fields:
		fname = get_hatids_in_field_fname(f)
		list_of_ids = pickle.load(open(fname, 'rb'))
		logprint( " get_candidates: Found %d hatids in field %s"%(len(list_of_ids), f))
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
if not args.hatids is None and len(args.hatids) > 0:
	logprint(" get_candidates: You gave %d hatids on the command line"%(len(args.hatids)))
	for hatid in args.hatids:
		hatids_to_process.append(hatid)


# ====================================================================
# ------------------------------------------------------------- /INPUT

# ORGANIZE ===========================================================
# --------------------------------------------------------------------

# Set the global hatid_field_list variable
hatid_field_list = set_hatid_field_list(all_fields)

logprint(" get_candidates: loaded hatid_field_list", all_nodes = True)

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

	return { f : [ hatid for hatid in hatids if not hatid in nofield_hatids and hatid_field_list[hatid] == f ] for f in fields if not f is None }, nofield_hatids

logprint(" get_candidates: sorting hatids by field")

# Split the hatids into their respective fields to save time
field_split, nofield_hatids = split_hatids_into_fields(hatids_to_process)

# Start a list of 'bad' hatids
BAD_IDS = []
BAD_IDS.extend(nofield_hatids)

# Save the hatids with no known field
if len(nofield_hatids) > 0:
	logprint(" get_candidates: %d / %d hatids have no field"%(len(nofield_hatids), len(hatids_to_process)))
	i = 0
	nofield_fname_gen = lambda i : "nofield_hatids_%d.pkl"%(i)
	while os.path.exists(nofield_fname_gen(i)): i += 1
	nofield_fname = nofield_fname_gen(i)
	logprint(" get_candidates: Saving list of hatids with no known field to %s"%(nofield_fname))
	safe_save(nofield_hatids, nofield_fname)

# ====================================================================
# ---------------------------------------------------------- /ORGANIZE


# PROCESS ============================================================
# --------------------------------------------------------------------

def process_field_serial(field):
	hatids = field_split[field]
	bad_ids_in_field = [ hatid for hatid in BAD_IDS ]
	candidates_in_field = []

	logprint(" get_candidates: Handling the %d hatids in field %s"%(len(hatids), field), all_nodes=True)
	
	# Load keylist for field + check that its not None
	keylist = load_keylist(field)
	if keylist is None:
		logprint("               : keylist for field %s is None; %d hatids are going to be labeled BAD"%(field, len(hatids)), all_nodes=True)
		bad_ids_in_field.extend(hatids)
		return bad_ids_in_field, candidates_in_field

	# Load twomass info for field + check that its not None
	twomass = twomass_info_file(field)
	if twomass is None:
		logprint("               : twomass info for field %s is None; %d hatids are going to be labeled BAD"%(field, len(hatids)), all_nodes=True)
		bad_ids_in_field.extend(hatids)
		return bad_ids_in_field, candidates_in_field

	# How many hatids do not have available twomass information?
	nno2mass = 0
	for hatid in hatids:
		if hatid not in bad_ids_in_field and (hatid not in twomass or twomass[hatid] is None):
			nno2mass += 1
			bad_ids_in_field.append(hatid)

	# Remove bad hatids
	hatids = [ hatid for hatid in hatids if not hatid in bad_ids_in_field ]

	# Print the number of hatids discarded because of missing 2mass information
	if nno2mass > 0: logprint(" get_candidates: 2mass info not available for %d hatids in field %s"%(nno2mass, field), all_nodes=True)

	# Skip this field if there are no good hatids.
	if len(hatids) == 0: 
		logprint(" get_candidates: there are no more viable hatids in field %s"%(field), all_nodes=True)
		return bad_ids_in_field, candidates_in_field

	for hatid in hatids:
		if hatid in bad_ids_in_field: continue

		# Generate lightcurve ... 
		lc = load_full_tfalc_from_scratch(hatid, field=None, 
				keylist_dat=keylist, twomass_dat=twomass[hatid], save_full_lc=True, 
				min_observations=5, delete_raw_lc=args.rm_orig, force_redo = args.force_redo)

		#  ... and check that the lc is OK
		if lc is None:
			bad_ids_in_field.append(hatid)
			continue

		# Generate features ...
		feats = generate_features(hatid)

		# ... and check that the features are ok
		if feats is None:
			bad_ids_in_field.append(hatid)
			continue

		# Score the hatid ...
		scores, status = test_hatid(hatid, model_prefix, min_score, min_frac_above_min_score, iteration, N=nmc)

		# ... and check that the scoring went OK
		if status is None:
			bad_ids_in_field.append(hatid)
		elif status:
			candidates_in_field.append(hatid)
		safe_save(scores, get_scores_fname(hatid,iteration))
		
	return bad_ids_in_field, candidates_in_field
def process_field_parallel(field):
	
	# Get a list of hatids in this field
	hatids = field_split[field]
	bad_ids_in_field = [ hatid for hatid in BAD_IDS ]
	candidates_in_field = []
	
	logprint(" get_candidates: Handling the %d hatids in field %s"%(len(hatids), field))

	# Load keylist and twomass information for this field
	keylist, twomass = None, None
	if ROOT:
		logprint(" get_candidates: Loading keylist...")
		keylist = load_keylist(field)
		logprint(" get_candidates: Loading twomass...")
		twomass = twomass_info_file(field)

	# Broadcast keylist + 2mass
	logprint(" get_candidates: broadcasting keylist and twomass")
	keylist = comm.bcast(keylist, root=0)
	twomass = comm.bcast(twomass, root=0)

	# Test twomass/keylist info
	if keylist is None:
		logprint(" get_candidates: keylist for field %s is None; %d hatids are going to be labeled BAD"%(field, len(hatids)))
		bad_ids_in_field.extend(hatids)
		return bad_ids_in_field, candidates_in_field
	if twomass is None:
		logprint(" get_candidates: 2mass info not available for field %s; %d hatids are going to be labeled BAD"%(field, len(hatids)))
		bad_ids_in_field.extend(hatids)
		return bad_ids_in_field, candidates_in_field

	# Are any hatids missing from the 2mass information?
	nno2mass = 0
	for hatid in hatids:
		if hatid not in bad_ids_in_field and (hatid not in twomass or twomass[hatid] is None):
			nno2mass += 1
			bad_ids_in_field.append(hatid)

	# Remove bad hatids
	hatids = [ hatid for hatid in hatids if not hatid in bad_ids_in_field ]

	# Print the number of hatids discarded because of missing 2mass information
	if nno2mass > 0: logprint(" get_candidates: 2mass info not available for %d hatids in field %s"%(nno2mass, field))

	# Skip this field if there are no good hatids.
	if len(hatids) == 0: 
		logprint(" get_candidates: there are no more viable hatids in field %s"%(field))
		return bad_ids_in_field, candidates_in_field

	# Master/slave workload distribution to make & test features
	if ROOT:
		logprint(" get_candidates: Handling the %d hatids in field %s"%(len(hatids), field))

		# Generate lightcurves
		logprint("               :  processing lightcurve!")
		results = msl.master(hatids)
		for ID, status in results:
			if status is False and not ID in bad_ids_in_field:
				bad_ids_in_field.append(ID)

		# Generate features
		logprint("               :  generating features!")
		results = msl.master(hatids)
		for ID, status in results:
			if status is False and not ID in bad_ids_in_field:
				bad_ids_in_field.append(ID)
				

		# Classify ids
		logprint("               :  scoring hatids!")
		results = msl.master(hatids)
		for ID, scores, status in results:
			if status is None and not ID in bad_ids_in_field:
				bad_ids_in_field.append(ID)	
			elif status:
				candidates_in_field.append(ID)
			safe_save(scores, get_scores_fname(ID,iteration))
		

		
	else:
		msl.slave(lambda hatid : (hatid, load_full_tfalc_from_scratch(hatid, field=None, 
				keylist_dat=keylist, twomass_dat=twomass[hatid], save_full_lc=True, 
				min_observations=5, delete_raw_lc=args.rm_orig, force_redo = args.force_redo)))
		msl.slave(lambda hatid : (hatid, generate_features(hatid)))
		msl.slave(lambda hatid : (hatid, ) + test_hatid(hatid, model_prefix, min_score, min_frac_above_min_score, iteration, N=nmc))

	logprint("               : done.")
	return bad_ids_in_field, candidates_in_field

# Choose the method of parallelization 
CANDIDATES = []
if parallelize_fields:
	logprint(" get_candidates: parallelizing FIELDS (nfields ~ nhatids)")
	if ROOT:
		bad_ids, cands = msl.master(field_split.keys())
		for bids, cs in zip(bad_ids, cands):
			BAD_IDS.extend(bids)
			CANDIDATES.extend(cands)
	else:
		msl.slave(process_field_serial)
else:
	logprint(" get_candidates: parallelizing HATIDS (nhatids >>> nfields)")
	for field in field_split:
		bad_ids, cands = process_field_parallel(field)
		if ROOT:
			BAD_IDS.extend(bad_ids)
			CANDIDATES.extend(cands)

# ====================================================================
# ----------------------------------------------------------- /PROCESS

# Save results.
if ROOT:
	safe_save(BAD_IDS, get_bad_ids_fname(iteration))
	safe_save(CANDIDATES, get_candidate_fname(iteration))

	logprint("               : %d / %d ids were 'bad'"%(len(BAD_IDS), len(hatids_to_process)))
	logprint("               : %d / %d ids were 'candidates'"%(len(CANDIDATES), len(hatids_to_process)))
