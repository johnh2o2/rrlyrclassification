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
from contextlib import closing
from paramiko import SSHConfig, SSHClient
from mpi4py import MPI

# Get rank and size of mpi process
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
ROOT = (rank == 0)


logprint(" process: MPI: rank = %d size = %d"%(rank, size))

# ====================================================================
# INPUT ==============================================================

# Define command line options
parser = argparse.ArgumentParser(description="Process lightcurves for a given field, hatid, or list of hatids")

parser.add_argument("hatids", nargs='*', help="One or more hatids to process")
parser.add_argument("--fields", nargs='*',help="Field(s) to process")
parser.add_argument("--list",nargs='*', help="Path to pickled list of hatids")
parser.add_argument("--make-features", action='store_true', help='Make features for each hatid')
parser.add_argument("--rm-orig", action='store_true', help='Remove original .tfalc lightcurve')
parser.add_argument("--force-redo", action='store_true', help='Not exactly sure what this does...')

# Parse command line arguments
args = parser.parse_args()

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
		logprint(" process: Found %d hatids in list %s"%(len(list_of_ids), l))
		hatids_to_process.extend(list_of_ids)

# Add specific hatids
if not args.hatids is None:
	logprint(" process: You gave %d hatids on the command line"%(len(args.hatids)))
	for hatid in args.hatids:
		hatids_to_process.append(hatid)
# ====================================================================
# ============================================================= /INPUT


# Remove any repeats
hatids_to_process = np.unique(hatids_to_process)

logprint(" process: Loading list of ALL hatids")
# Load the hatid <-> field mapping
hatid_field_list = load_hatid_field_list(all_fields)


def split_hatids_into_fields(hatids):
	'''
	splits a list of hatids into a dictionary of lists
	
	returns { field1 = [ hatids in field 1 ], ... }
	'''
	nofield_hatids = [ hatid for hatid in hatids if not hatid in hatid_field_list ]
	nofield_hatids.extend([ hatid for hatid in hatids if hatid in hatid_field_list and hatid_field_list[hatid] is None ])

	fields = np.unique( [ hatid_field_list[hatid] for hatid in hatids if not hatid in nofield_hatids ] )

	return { f : [ hatid for hatid in hatids if hatid_field_list[hatid] == f ] for f in fields if not f is None }, nofield_hatids

logprint(" process: sorting hatids by field")
# Split the hatids into their respective fields to save time
field_split, nofield_hatids = split_hatids_into_fields(hatids_to_process)

# Save the "bad" hatids
if len(nofield_hatids) > 0:
	logprint(" process: %d / %d hatids have no field"%(len(nofield_hatids), len(hatids_to_process)))
	i = 0
	nofield_fname_gen = lambda i : "nofield_hatids_%d.pkl"%(i)
	while os.path.exists(nofield_fname_gen(i)): i += 1
	nofield_fname = nofield_fname_gen(i)
	logprint(" process: Saving list of hatids with no known field to %s"%(nofield_fname))
	pickle.dump(nofield_hatids, open(nofield_fname, 'wb'))

keylists, twomass = {}, {}
# Load keylist and twomass information for fields
if ROOT:
	BAD_IDS = []
	BAD_IDS.extend(nofield_hatids)

	logprint(" process: Loading keylists.")
	keylists_ = msl.master([ field for field in field_split])
	logprint(" process: Loading 2mass information")
	twomass_ = msl.master([ field for field in field_split ])

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

# Now process the hatids
for field in field_split:
	hatids = field_split[field]
	logprint(" process: Handling the %d hatids in field %s"%(len(hatids), field))

	# Master/slave workload distribution to make & test features
	if ROOT:

		# Generate lightcurves
		logprint("        :  processing lightcurve!")
		results = msl.master(hatids)
		for ID, status in results:
			if status is False and not ID in bad_ids:
				BAD_IDS.append(ID)

		
		if not args.make_features: continue

		# Generate features
		logprint("        :  generating features!")
		results = msl.master(hatids)
		for ID, status in results:
			if status is False and not ID in bad_ids:
				BAD_IDS.append(ID)

		
	else:

		msl.slave(lambda hatid : (hatid, load_full_tfalc_from_scratch(hatid, field=None, 
				keylist_dat=keylists[field], twomass_dat=twomass[field][hatid], save_full_lc=True, 
				min_observations=5, delete_raw_lc=args.rm_orig, force_redo = args.force_redo)))
		
		if not args.make_features: continue
		msl.slave(lambda hatid : (hatid, generate_features(hatid)))

	logprint("        : done.")

bad_ids_fname = "bad_ids_from_process_script_%d.pkl"%(i)
pickle.dump(BAD_IDS, open(bad_ids_fname, 'wb'))
logprint("        : COMPLETELY DONE! %d / %d ids were 'bad' and stored in %s"%(len(BAD_IDS), len(hatids_to_process), bad_ids_fname))

