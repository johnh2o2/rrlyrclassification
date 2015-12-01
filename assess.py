import settings
import numpy as np
import utils.readhatlc as rhlc
import os, sys, gzip
import masterslave as msl
import utils.miscutils as mutils
import matplotlib.pyplot as plt
import utils.featureutils as futils
import cPickle as pickle
from time import time
import argparse
from mpi4py import MPI
from math import *


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
ROOT = (rank == 0)

print_keylist_commands = False

def time_function(func, *args, **kwargs):
	t0 = time()
	result = func(*args, **kwargs)
	dt = time() - t0
	return t0, result

all_fields = \
[
	'145', '219', '216', '214', '215', '212', '213', '210', '211', '093', 
	'095', '096', '133', '132', '136', '135', '134', '138', '161', '288', 
	'341', '342', '163', '285', '284', '287', '286', '123', '124', '125', 
	'126', '127', '128', '269', '268', '378', '416', '294', '292', '293', 
	'377', '376', '318', '199', '198', '195', '194', '311', '196', '191', 
	'190', '193', '192', '115', '088', '277', '142', '143', '207', '206', 
	'366', '364', '365', '362', '300', '247', '384', '240', '242', '388', 
	'389', '100', '101', '249', '248', '162', '241', '432', '431', '144', 
	'258', '259', '177', '176', '175', '174', '173', '257', '170', '203', 
	'182', '183', '186', '187', '184', '185', '188', '063', '062', '317', 
	'316', '168', '315', '164', '165', '166', '167', '160', '222', '221', 
	'220', '314', '151', '150', '153', '152', '155', '154', '159', '238', 
	'239', '235', '236', '146', '147', '089', '205', '204', '140', '141', 
	'209', '087', '086', '148', '149', '357', '267' 
]

# Does the lightcurve have a known field?
def has_known_field(hatid):
	return ( not get_field_of(hatid) is None )

# Does the lightcurve exist on the remote server?
def remote_lightcurve_exists(hatid):
	return os.path.exists(get_remote_tfalc_filename(hatid))

# Does the raw lightcurve exist?
def raw_lightcurve_exists(hatid):
	return os.path.exists(get_raw_lc_fname(hatid))

# Does the gzipped lightcurve exist?
def gzipped_lightcurve_exists(hatid):
	return os.path.exists(get_lc_fname(hatid))

# Is the gzipped lightcurve none?
def gzipped_lightcurve_is_none(hatid):
	return (pickle.load(gzip.open(get_lc_fname(hatid), 'rb')) is None)

# How many observations does the raw lightcurve have?
def num_obs_raw(hatid):
	lc = mutils.load_tfalc(get_raw_lc_fname(hatid))
	if lc is None: return np.nan
	if not 'BJD' in lc: return np.nan
	return len(lc['BJD'])

# Does the keylist exist for the field?
def does_keylist_exist(hatid):
	field = get_field_of(hatid)
	if field is None: return False
	kl = load_keylist(field)
	if kl is None: return False
	else: return True

# How many observations are in the keylist for the field?
def num_obs_in_keylist(hatid):
	fname = get_lc_fname(hatid)
	if not os.path.exists(fname): return np.nan
	lc = pickle.load(gzip.open(fname, 'rb'))
	if lc is None: return np.nan
	if not 'BJD' in lc: return np.nan
	return len(lc['BJD'])

# Does twomass info exist for field?
def field_has_2mass_info(hatid):
	field = get_field_of(hatid)
	if field is None: return False
	tm = twomass_info_file(field)
	if tm is None: return False
	
# Does the lightcurve have twomass info?
def has_2mass_info(hatid):
	field = get_field_of(hatid)
	if field is None: return False
	tm = twomass_info_file(field)
	if tm is None: return False
	if hatid in tm: return True

def addtxt(a, b):
	return "%s\n%s"%(a, b)

def get_keylist_commands(fields):
	a = ""
	a = addtxt(a, "export HATPIPE=/home/hatuser/HATpipebin")
	a = addtxt(a, "export PATH=${PATH}:/home/hatuser/HATpipebin/bin")

	for field in fields:
		a = addtxt(a, "00_keylist_gen.sh %s 0 > keylist_%s_newscript.txt"%( field, field ))

	for field in fields:
		a = addtxt(a, "gen_keylist.sh %s > keylist_%s_oldscript.txt"%( field, field ))
	return a 

# Can we find additional keylist info for the lightcurve using Joel's script?
# Can we find twomass info for the lightcurve using the 2mass__ script?

def assess_gzipped_lc(gzipped_fname):
	# try to load gzipped lc
	if os.path.exists(gzipped_fname):
		tgzlc_open, fstream = time_function(gzip.open, gzipped_fname, 'rb')
		tgzlc_load, lc = time_function(pickle.load, fstream)
		t_gzipped_lc_open = tgzlc_open + tgzlc_load
	else:
		t_gzipped_lc_open, lc = 0, None


	a = {
		'gzipped_lc_is_none' : lc is None,	
		'time_to_load_gzipped_lc' : t_gzipped_lc_open
	}

	if not a['gzipped_lc_is_none']:
		if not 'BJD' in lc:
			a['nobs_processed'] = None
		else:
			a['nobs_processed'] = len(lc['BJD'])

	return a

def assess_raw_lc(raw_fname):

	# try to load raw lightcurve
	if os.path.exists(raw_fname):
		traw, raw_lc = time_function(mutils.load_tfalc, raw_fname)
	else:
		traw, raw_lc = 0, None

	a = {
		'raw_lc_is_none' : raw_lc is None,
		'time_to_load_raw_lc' : traw
	}

	if not a['raw_lc_is_none']:
		if not 'BJD' in raw_lc: 
			a['nobs_raw'] = None,
		else:
			a['nobs_raw'] = len(raw_lc['BJD'])
	return a

def assess_features(features_fname):
	# try to load features
	features = None
	if os.path.exists(features_fname):
		t1, fstream = time_function(open, features_fname, 'rb')
		t2, features = time_function(pickle.load, fstream)
		tfeatures = t1 + t2
	else:
		tfeatures = 0
	a = {
		'time_to_load_features' : tfeatures,
		'features_are_none' : features is None
	}
	if not a['features_are_none']:
		a.update({
				'vmag'			: features['V'],
				'std'			: features['std'],
				'max_lsp_power' : features['raw_lsp_peak1_power']
			})
	return a
		

def assess_hatid(hatid, keylist, twomass, load_lightcurves=False, load_features=True):
	raw_fname = settings.get_raw_lc_fname(hatid)
	gzipped_fname = settings.get_lc_fname(hatid)
	features_fname = settings.hat_features_fname(hatid)

	a = {
		'raw_fname' : raw_fname,
		'gzipped_fname' : gzipped_fname,
		'features_fname' : features_fname,

		'has_raw_lc' : os.path.exists(raw_fname),
		'has_gzipped_lc' : os.path.exists(gzipped_fname),
		'has_features' : os.path.exists(features_fname),
		'has_2mass_info' : hatid in twomass,
	}

	if load_features: 
		a.update(assess_features(features_fname))
		
	if load_lightcurves: 
		a.update(assess_raw_lc(raw_fname))
		a.update(assess_gzipped_lc(gzipped_fname))

	return a 
def strdict(d):
	a = ""
	for key, value in d.iteritems():
		if a == "": a = "{0} = {1}".format(key, value)
		else: a = "{0}; {1} = {2}".format(a, key, value)
	return a

def field_assessment(field, assess_each_hatid=False, load_lightcurves=False, load_features=True):
	mutils.logprint(" assessing field %s"%(field), all_nodes=True)
	
	# Load hatids in that field
	mutils.logprint("     %s: loading hatids "%(field), all_nodes=True)
	hatids_fname = os.path.join(settings.hatids_in_fields_dir, "hatids_in_field_%s.list"%(field))
	if os.path.exists(hatids_fname):
		thatid_list_open, fstream = time_function(open, hatids_fname, 'rb')
		thatid_list_pload, hatids = time_function(pickle.load, fstream)
	else:
		hatids = None
		thatid_list_open, thatid_list_pload = 0, 0

	mutils.logprint("     %s: loading keylist "%(field), all_nodes=True)
	# Load keylist for the field
	tkeylist, keylist = time_function(mutils.load_keylist, field)

	mutils.logprint("     %s: loading twomass "%(field), all_nodes=True)
	# Load the twomass_info for the field
	ttwomass, twomass_info = time_function(settings.twomass_info_file,field)

	mutils.logprint("     %s: setting up assessment. "%(field), all_nodes=True)
	# Record some basic assessments of the field
	field_assessment = {
		'directory' : mutils.field_info[field],
		'has_hatid_field_list' : os.path.exists(hatids_fname),
		'field_has_2mass_info' : not (twomass_info is None),
		'does_keylist_exist' : not keylist is None,
		'hatids' : None,
		'num_hatids' : None,
		'hatid_list_is_ok' : isinstance(hatids, list) and len(hatids) > 0,
		'time_to_load_hatid_list' : thatid_list_open + thatid_list_pload,
		'time_to_load_keylist' : tkeylist,
		'time_to_load_2mass' : ttwomass,
	}

	if field_assessment['hatid_list_is_ok']: 
		field_assessment['num_hatids'] = len(hatids)

	# if any of the (twomass info, keylist, hatid list) are missing, this is a 'bad' field and we're done.
	bads = { key : not value for key, value in field_assessment.iteritems()\
					 if key in [ 'hatid_list_is_ok', 'field_has_2mass_info', 'does_keylist_exist']}
	if any(bads): 

		mutils.logprint("     %s: %s"%(field, strdict(bads)), all_nodes=True)
		return field_assessment

	# Or, if we're not looking at each individual hatid, we're done.
	if not assess_each_hatid:
		mutils.logprint("     %s: done (not assessing indiv. hatids)"%(field), all_nodes=True)
		return field_assessment

	mutils.logprint("     %s: assessing indiv. hatids")
	field_assessment['hatids'] = { hatid : assess_hatid(hatid, load_lightcurves=load_lightcurves, load_features=load_features) for hatid in hatids }
	mutils.logprint("     %s: DONE!")

	return field_assessment

assessment_dir =  os.path.join(settings.SCRATCH, 'assessments')
def make_assessment_dir():
	if not os.path.isdir(assessment_dir) and ROOT: os.makedirs(assessment_dir)
fld_assmnt_name = lambda field : 'assessment_field_%s.pkl.gz'%(field)
field_assessment_fname = lambda field : os.path.join(assessment_dir, fld_assmnt_name(field))


if __name__ == '__main__':

	# Handle arguments
	parser = argparse.ArgumentParser(description="assess the lightcurves in a given field (or set of fields)")
	parser.add_argument("--all-fields", action='store_true', help="analyze all known fields (%d of them); overrides --fields"%(len(all_fields)))
	parser.add_argument("--fields", nargs='*',help="Field(s) to process")
	parser.add_argument("--read-lightcurves", action='store_true', help="Reads raw and gzipped lightcurves to give more info; takes much longer!")
	parser.add_argument("--dont-read-features", action='store_false', dest='read_features', help="Don't read features; this may save time but you dont get as many features")
	parser.add_argument("--output-root", help="string distinguishing this particular round of assessments" )
	parser.add_argument("--redo", action="store_true", help="Redo assessments even if there's an existing file (overwrites existing files!)")
	args = parser.parse_args()

	# Get fields
	if args.all_fields:
		fields = all_fields
	else:
		fields = args.fields

	if fields is None:
		raise Exception("No fields specified!")

	# Set up filenames
	if not args.output_root is None:
		if '/' in args.output_root: 
			field_assessment_fname = lambda field : "%s/%s"%(args.output_root, fld_assmnt_name(field))
		else:
			fld_assmnt_name_new = lambda field : "%s_%s"%(args.output_root, fld_assmnt_name(field))
			field_assessment_fname = lambda field : os.path.join(assessment_dir, fld_assmnt_name_new(field))

	save_field_assessment = lambda field, assessment : pickle.dump(assessment, gzip.open(field_assessment_fname(field), 'wb'))

	# Prints and doesn't return things!
	def process_and_print(field, **kwargs):
		a = field_assessment(field, **kwargs)
		mutils.logprint("     %s: Saving!!"%(field))
		save_field_assessment(field, a)
		#return a

	if ROOT:
		make_assessment_dir()
		results = msl.master(fields)

	else:
		msl.slave(process_and_print)




