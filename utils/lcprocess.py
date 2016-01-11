'''

utilities for processing HAT lightcurves into gzipped python dictionaries

John Hoffman 

'''

import os, io, sys, gzip
import numpy as np
import cPickle as pickle
from time import time
import logging
from config import *

logger = logging.getLogger(__name__)

#### ERROR CODES ####
errlist = [ 
	'OK' , 'RAW_LC_FILE_NOT_FOUND', 'TOO_FEW_OBSERVATIONS', 'PROBLEM_READING_RAW_LC_FILE', 'NO_HATID_LIST_FOUND', 
	'NO_RAW_TWOMASS_FILE_FOUND', 'NO_RAW_KEYLIST_FILE_FOUND', 'NO_RAW_LCSTAT_FILE_FOUND', 'NO_KEYLIST_DICTIONARY_FOUND',
	'NO_TWOMASS_DICTIONARY_FOUND', 'NO_LCSTAT_DICTIONARY_FILE_FOUND', 'NO_MISSING_EXPOSURES_FILE_FOUND', 
	'UNKNOWN_LIGHTCURVE_FORMAT', 'KEYLIST_CONTAINS_NULL_VALUES', 'CANNOT_NPLOADTXT_RAW_LIGHTCURVE', 'CANNOT_FIND_LCDIR' 
]


def add_keylist_data(lc, keylist, **kwargs):
	""" adds keylist data to the lightcurve (lc) dictionary.
		necessary for getting the dates of each exposure
	"""
	cols = [ 'FLT', 'EXP', 'BJD' ]

	# HATnorth sources have a slightly different format than HATsouth sources;
	#  north TSTFC = <TSTF>_<CCD>
	#  south TSTF  = <TSTF>

	if 'HAT' in lc: 
		is_hatnorth_source = True
	else:
		is_hatnorth_source = False
		
	# initialize
	for col in cols: lc[col] = []

	missing_exposures = set()
	# For each measurement
	for i in range(len(lc['TF1'])):

		# Obtain CCD/exposure codes
		if is_hatnorth_source:
			EXP, CCD = lc['TSTFC'][i].split('_')
		else:
			EXP = lc['TSTF'][i]

		not_available = (EXP not in kl)
		if not_available: missing_exposures.add(EXP)
			
		# Add info from keylist!
		for col in cols:
			if not_available:
				lc[col].append('?')
			else:
				lc[col].append(kl[EXP][col])

	# Convert to numpy arrays
	for col in cols:
		lc[col] = np.array(lc[col])

	return 'OK', missing_exposures

def add_twomass_data(lc, twomass, **kwargs):
	""" adds twomass information to 'lc' (lightcurve dictionary) 
	"""
	lc['info'] = { key : value for key,value in twomass.iteritems() if key in kwargs['twomass_columns_to_use']  }
	#lc['ra'] 		= twomass['ra']
	#lc['dec'] 		= twomass['dec']
	#lc['mags'] 		= [ twomass[x] for x in ['Vmag', 'Rmag', 'Imag', 'jmag', 'hmag', 'kmag'] ]
	#lc['ugriz'] 	= [ twomass[x] for x in ['umag', 'gmag', 'rmag', 'imag', 'zmag' ] ]
	#lc['ndet'] 		= len(lc['RJD'])
	
	#lc['hatstations'] = np.unique(lc['STF'])
	#lc['filters'] 	= [ flt for flt in np.unique(lc['FLT']) ]

	return 'OK'

def load_raw_keylist(field, **kwargs):
	""" returns a numpy ndarray containing the keylist data
		for each exposure of the field
	"""

	# load things as usual
	filename = get_keylist_filename(field, **kwargs)
	if not os.path.exists(filename): return None, 'NO_RAW_KEYLIST_FILE_FOUND'

	
	# Some of the keylists are corrupted. This try/except will 
	# avoid breaking things unless necessary
	try:
		keylist = np.loadtxt(filename, dtype = keylist_dt)
		
	except ValueError, e:
		if kwargs['clean_nulls_from_keylist']:
			# TRY to clean out NULLs -- this will throw away 
			# otherwise good exposures and should be avoided.
			logger.warning("[field %s]; keylist file likely contains NULLs :( " %(field))
			logger.info("cleaning keylist")
			cleaned_keylist = clean_nulls(filename)
			logger.info("done.")
			keylist = np.loadtxt(io.StringIO(cleaned_keylist), dtype = keylist_dt)
		else:
			return None, 'KEYLIST_CONTAINS_NULL_VALUES'
		
	return keylist, 'OK'

def make_hatid_list(field, **kwargs):
	""" generates and saves a list of hatids located in the
		directory for a given field

		returns list of hatids and OK or None and an error msg
	"""

	lcdir = get_lc_dir(field, **kwargs) 
	if not os.path.isdir(lcdir): return None, 'CANNOT_FIND_LCDIR'

	# Get a list of all files in lcdir
	files = os.listdir(lcdir)
	# Prune out ones that arent .tfalc lightcurves
	files = [ f for f in files if is_a_lightcurve(f) ]
	# Translate these filenames to hatids
	hatids = [ f.split('.')[0] for f in files ]
	# Prune out ones that dont have the right length
	hatids = set([ hatid for hatid in hatids if len(hatid) == 15 ])

	pickle.dump(hatids, open(get_hatid_list_filename(field, **kwargs), 'wb'))

	return hatids, 'OK'



def combine_all_eras(hatid, **kwargs):
	""" combine all available lightcurves for a given hatid
		into a single file
	"""
	pass

def read_raw_lc_file(lc_filename, **kwargs):
	""" reads a HatNet .tfalc lightcurve.
		returns a numpy ndarray
	"""

	# Read in raw lightcurve text
	lcfile = raw_lc_open_func(lc_filename, 'r')
	rawtxt = lcfile.read()
	lines = rawtxt.split('\n')
	
	# Determine number of columns
	i = 0
	while '#' in lines[i]: i+=1
	ncols = len(lines[i].split())

	# Use ncols to determine the correct dtype (HS or HN)
	if ncols == len(tfalc_hs_arr): dtype = tfalc_hs_dt
	elif ncols == len(tfalc_hn_arr): dtype = tfalc_hn_dt
	else: return None, 'UNKNOWN_LIGHTCURVE_FORMAT'

	# Attempt to load with np.loadtxt
	try: 
		lcraw = np.loadtxt(io.StringIO(rawtxt), dtype=dtype)
		return lcraw, OK
	except ValueError, e:
		return None, 'CANNOT_NPLOADTXT_RAW_LIGHTCURVE'


def load_raw_lc(hatid, lcdir, **kwargs):
	""" loads a raw (gzipped ascii) lightcurve (usually tfalc).
	"""
	raw_lc_filename = get_raw_lc_filename(lcdir, hatid)
	if not os.path.exists(raw_lc_filename):
		return None, 'RAW_LC_FILE_NOT_FOUND'

	return read_raw_lc_file(raw_lc_filename, **kwargs)


def process_lightcurve(hatid, keylist, twomass, lcdir, **kwargs):
	""" processes a HAT lightcurve, or returns an error code if 
		a known problem is encountered. Returns OK if processed 
		successfully
	"""
	# Read the raw tfalc file
	logger.info("[%s] in process_lightcurve"%(hatid))
	logger.info("[%s] loading raw lightcurve"%(hatid))
	lc, msg = load_raw_lc(raw_lc_fname, **kwargs)
	if not msg is 'OK': return msg

	logger.info("[%s] converting raw lightcurve to dictionary"%(hatid))
	# convert to a dictionary
	lc = { c : lc[c].tolist() for c in lc and c in kwargs['lightcurve_columns_to_use'] }

	if kwargs['add_keylist_data']:
		logger.info("[%s] adding keylist data"%(hatid))
		# Add keylist data
		msg, new_missing_exposures = add_keylist_data(lc, keylist)
		if len(new_missing_exposures) > 0:
			logger.info("[%s] %d exposures thrown out (missing)"%(hatid, len(new_missing_exposures)))
		if not msg is 'OK': return msg

	if kwargs['add_twomass_data']:
		# Add twomass data
		logger.info("[%s] adding twomass data"%(hatid))
		msg = add_twomass_data(lc, twomass)
		if not msg is 'OK': return msg

	# Save the processed lightcurve
	logger.info("[%s] saving processed lightcurve")
	if 'save' in kwargs and kwargs['save']:
		msg = save_processed_lc(lc, get_processed_lc_filename(lcdir, hatid))

	# Delete the raw lc file if requested
	if 'save' in kwargs and kwargs['save'] and kwargs['delete_raw_lc_after_processing']:
		logger.info("[%s] deleting raw lightcurve (%s)"%(hatid, raw_lc_filename))
		os.remove(raw_lc_filename)

	return OK

def process_field(field, **kwargs):
	""" Processes the hatids for an entire HAT field. Returns a
		dictionary of error codes (indexed by hatid) for lightcurves
		that encountered problems during processing.
	"""
	logger.info("in process_field")

	lcdir = get_lc_dir(field, **kwargs) 

	if kwargs['add_keylist_data']:
		logger.info("fetching keylist")
		keylist, msg = fetch_keylist(field, **kwargs)
		logger.info("fetching missing_exposures")
		missing_exposures, msg = fetch_missing_exposures(field, **kwargs)

	if kwargs['add_twomass_data']:
		logger.info("fetching twomass data")
		twomass, msg = fetch_twomass(field, **kwargs)

	if not 'hatids' in kwargs:
		logger.info("fecthing hatids")
		hatids, msg = fetch_hatids(field, **kwargs)
		logger.info("%d hatids found for field %s"%(len(hatids), field))
	else:
		hatids = kwargs['hatids']

	problems = {}
	if kwargs['use_mpi_on_single_field']:
		from mpi4py import MPI
		import masterslave as msl

		comm = MPI.COMM_WORLD
		size = comm.Get_size()
		rank = comm.Get_rank()
		ROOT = (rank == 0)

		problems = None
		if ROOT:
			mpi_formatter  = logging.Formatter('NODE %-4d: %s'%(rank, format_string))
			console.setFormatter(mpi_formatter)

			results = msl.master(hatids)
			problems = { hatid : result for hatid, result in results if not result is 'OK' }

		else:
			msl.slave(lambda hatid : (hatid, process_lightcurve(hatid, keylist, twomass, lcdir,save=True, **kwargs)))

		problems = comm.bcast(problems, root=0)
		return problems

	else:
		for hatid in hatids:
			msg = process_hatid(hatid, keylist, twomass, lcdir, **kwargs)
			if not msg is 'OK': 
				logger.warning('%s has problem: %s'%(hatid, msg))
				problems[hatid] = msg

	return problems

def check_field(field, **kwargs):
	""" Makes sure all relevant data exists for a given hat field.
	"""
	logger.info('in check_field')
	lcdir = get_lc_dir(field, **kwargs) 
	logger.info('got lcdir: %s'%(lcdir))
	if not os.path.isdir(lcdir):
		return 'CANNOT_FIND_LCDIR'

	if kwargs['add_keylist_data']:
		logger.info('fetching keylist')
		keylist, msg = fetch_keylist(field, **kwargs)
		logger.info('done; result is %s'%(msg))
		if not msg is 'OK': return msg

		logger.info('fetching missing exposures')
		missing_exposures, msg = fetch_missing_exposures(field, **kwargs)
		logger.info('done; result is %s'%(msg))
	
	if kwargs['add_twomass_data']:
		logger.info('fetching twomass information')
		twomass, msg = fetch_twomass(field, **kwargs)
		logger.info('done; result is %s'%(msg))
		if not msg is 'OK': return msg

	if kwargs['use_lcstat_data']:
		logger.info('fetching lcstat data')
		lcstat, msg = fetch_lcstat(field, **kwargs)
		logger.info('done; result is %s'%(msg))
		if not msh is 'OK': return msg

	logger.info('fetching hatids')
	hatids, msg = fetch_hatids(field, **kwargs)
	logger.info('done; result is %s.')
	if not msg is 'OK': return msg
	
def action_functions(field, **kwargs):
	""" python dictionary of actions to be taken 
		for a series of recoverable error codes
		(e.g. error codes that are a result of missing processed
		files.)
	"""
	return {
				'NO_HATID_LIST_FOUND'             : lambda : make_hatid_list(field, **kwargs)
				'NO_KEYLIST_DICTIONARY_FOUND'     : lambda : make_keylist_dict(field, **kwargs)
				'NO_TWOMASS_DICTIONARY_FOUND'     : lambda : make_twomass_dict(field, **kwargs)
				'NO_LCSTAT_DICTIONARY_FILE_FOUND' : lambda : make_lcstat_dict(field, **kwargs)
			}

def preprocess_field(field, **kwargs):
	""" recursively generates files necessary for processing
		a HAT field.
	"""
	logger.info("in preprocess_field; checking %s"%(field))
	msg = check_field(field, **kwargs)
	logger.info("checking is done: result is %s"%(msg))

	if not msg in actions: 
		if not msg is 'OK': logger.error("field %s : %s"%(field, msg))
		return msg
	else: 
		action_functions(field, **kwargs)[msg]()
		preprocess_field(field, **kwargs)

def log_problems(problems):
	# Output information about each field
	for field, hatid_problems in problems.iteritems():

		# Count of number of hatid with each problem
		counts = { problem : 
					len([ hatid for hatid, p in hatid_problems.iteritems() if p == problem ]) \
						for problem in np.unique(hatid_problems.values()) }

		logger.info("Field %s had %d hatids that failed: "%(field, len(hatid_problems)))
		for problem, number in counts.iteritems():
			logger.info("    (field %s): %-10d hatids had problem %-25s"%(field, number, problem))


if __name__ == '__main__':

	# Create argument parser
	import argparse
	parser = argparse.ArgumentParser(desc='Process HAT lightcurves')
	parser.add_argument('--fields', nargs="*", help="fields to process")
	parser.add_argument('--config-file', 
		help="Path to a .ini file readable by the configparse python library; must contain lcprocess header.")
	parser.add_argument('--hatids', nargs="*", help="hatids to process")
	parser.add_argument('--hatid-list-pkl', help="Path to pickled file containing list of hatids to process")
	parser.add_argument('--hatid-list', help="Path to ascii file containing list of hatids to process")
	parser.add_argument('--logfile', help="Path to save logging results (optional)")

	# Allow user to specify arguments for each recognized setting
	add_settings_to_parser(parser, **default_settings)

	args = parser.parse_args()

	# Add file handler for log if needed
	if args.logfiles not None:
		file_handler = logging.FileHandler(args.logfile)
		logger.addHandler(file_handler)

	# initialize settings (as defualt)
	settings = { key : value['default'] for key, value in default_settings.iteritems() }

	# Update settings if the user specified a configuration file
	if not args.config_file is None:
		import ConfigParser

		# Read config file
		configfile = ConfigParser.ConfigParser()
		configfile.read(args.config_file)
		config_settings =  translate_configparser(configfile, 'lcprocess')

		# Update settings
		settings.update(config_settings)

	# Now build dictionary of fields to process, with a set() of hatids
	# for each field.
	process_dict = {}

	# Add entire fields
	if not args.fields is None:
		for field in args.fields: 
			add_hatids_to_dict(process_dict, fetch_hatids(field, **kwargs))
	
	# Add specific hatids
	if not args.hatids is None:
		add_hatids_to_dict(process_dict, args.hatids)

	# Add a pickled hatid list
	if not args.hatid_list_pkl is None:
		add_hatids_to_dict(process_dict, pickle.load(open(args.hatid_list_pkl, 'rb')))

	# Add an ascii list of hatids
	if not args.hatid_list is None:
		add_hatids_to_dict(process_dict, load_ascii_list(args.hatid_list))

	# Do things in parallel if requested to
	if settings['use_mpi'] and not settings['use_mpi_on_single_field']:
		# Import mpi/masterslave libraries
		from mpi4py import MPI
		import masterslave as msl

		# Get information about this node
		comm = MPI.COMM_WORLD
		size = comm.Get_size()
		rank = comm.Get_rank()
		ROOT = (rank == 0)

		# Change the formatter so that we know what each node is doing
		mpi_formatter  = logging.Formatter('NODE %-4d: %s'%(rank, format_string))
		console.setFormatter(mpi_formatter)

		if ROOT:
			all_fields = process_dict.keys()

			# Preprocess all of the fields
			preprocess_results = msl.master(all_fields)

			# Only deal with the fields that passed
			good_fields = []
			for field,result in zip(all_fields, preprocess_results):
				if not result is 'OK':
					logger.error("Field %s did not preprocess successfully: error was %s"%(field, result))
					continue
				good_fields.append(field)

			# Test if there are any fields that did OK or not
			if len(good_fields) == 0:
				logger.error("No preprocessed fields left to process.")
				sys.exit()

			# Process each preprocessed field
			process_results = msl.master(good_fields)

			# Combine all of the problematic hatids
			problems = { field : result for field, result in zip(good_fields, process_results) if not result is 'OK' }

			log_problems(problems)

		else:
			msl.slave(lambda field : preprocess_field(field, **settings))
			msl.slave(lambda field : process_field(field, hatids=process_dict[field], **settings))

	else:
		problems = {}
		for field, hatids in process_dict.iteritems():
			logger.info("preprocessing field %s"%(field))
			msg = preprocess_field(field, **settings)
			if not msg is 'OK':
				logger.error('preprocessing field %s has failed with message %s'%(field, msg))
				continue
			logger.info("processing %d hatids in field %s"%(len(hatids), field))
			problems[field] = process_field(field, hatids=hatids, **settings)

		log_problems(problems)
















