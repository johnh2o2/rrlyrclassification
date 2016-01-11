import numpy as np
import os, gzip
import logging

format_string = '%(asctime)s %(name)-15s [%(funcName)-30s]: %(levelname)-8s %(message)s'

console = logging.StreamHandler()
formatter = logging.Formatter(format_string)
console.setFormatter(formatter)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(console)


directory_structure = {
	'SCRATCH' : {
		'LCCACHE' : [ 'ERAS', 'FIELDS' ],
		'INFO'    : [ 'HATID_LISTS', 'KEYLISTS', 'LCSTAT', 'TWOMASS' ]
	}
}
default_settings = {
	'use_mpi_on_single_field' : dict(
		default = True, 
		help = 'This should be true only if number of fields is small compared to the number of available processors',
		action = 'store_true'
	),

	'use_mpi' : dict(
		default= True, 
		help = 'Recommended. Defaults to assigning a field to each processor until all fields have been analyzed.',
		action = 'store_true'
	),

	'add_twomass_data' : dict(
		default = True, 
		help = 'Adds color data from the twomass catalog as well as estimates for other colors',
		action = 'store_true'
	),

	'add_keylist_data' : dict(
		default = False, 
		help = 'Adds keylist information (BJD, CCD, etc) for each exposure. Not recommended.',
		action='store_true'
	),

	'dustmap' : dict(
		default = None,
		help = 'Does nothing at the moment. Might be better to do this when generating features.'
	),

	'use_lcstat_data' : dict(
		default = False, 
		action='store_true',
		help = 'Uses statistics computed for each hatid in the field. Not available or relevant at this stage.'
	),

	'SCRATCH' : dict(
		default='/Users/jah5/Documents/Gaspar/rrlyr_classification/SCRATCH_NEW', 
		help = 'Location of scratch directory (place to store lightcurves and other data)'
	),

	'deal_with_multiple_lcs' : dict(
		default = False, 
		action='store_true',
		help = 'Doesn\'t do anything yet; in the future this will combine all of the lightcurves for a single object'
	),

	'era_locations_file' : dict(
		default = '/Users/jah5/Documents/Gaspar/work/new_era_locations.dict', 
		help = 'Location to find the remote directory locations of each era'
	),

	'redo_all' : dict( 
		default = False, 
		help = 'Redo everything we can'
	),

	'delete_raw_lc_after_processing' : dict(
		default = False, 
		help = 'Deletes original (unprocessed) lightcurve -- saves space'
	),

	'lightcurve_columns_to_use' : dict(
		default = [ 'RJD', 'TF1', 'TF2', 'TF3', 'TSTF', 'TSTFC', 'FLD', 'FLT' ],
		nargs = '*',
		help = 'which columns of the lightcurve to keep'
	),
	'twomass_columns_to_use' : dict(
		default = [ 'ra','dec','Vmag','Rmag','Imag','jmag',
					'hmag','kmag', 'umag' ,'gmag','rmag','imag','zmag' ],
		nargs = '*',
		help = 'which columns of the twomass information to keep'
	)
}

##############################################

############## NUMPY DATA TYPES ##############

keylist_dt_arr = [ ('TSTF', 'S8'), ('EXP', float), ('BJD', float), ('FLT', str),
		('unknown1', 'S20'), ('unknown2', 'S20'), ('unknown3', 'S20'), 
		('unknown4', 'S20'), ('unknown5', 'S20'), ('unknown6', 'S20'), 
		('unknown7', 'S20'),
		# The S20's are really float's, but some keylists 
		# have NULL strings just shoved in random rows.
		# That makes things break. So, to avoid this, 
		# we treat the irrelevant parts of the keylist as strings.
	]
keylist_dt = np.dtype(keylist_dt_arr)
twomass_dt_arr = [ 	
	('hatid', 'S15'), ('ra', float), ('dec', float), 
	('jmag', float), ('jerr', float), ('hmag', float), ('herr', float),
	('kmag', float), ('kerr', float),
	('flags', 'S3'), ('Bmag', float), ('Vmag', float), ('Rmag', float), ('Imag', float),
	('umag', float), ('gmag', float), ('rmag', float), ('imag', float), ('zmag', float),
	('hatfield', int), ('hatfield objid', int)
]
twomass_dt = np.dtype(twomass_dt_arr)
lcstat_dt_arr = [
	('hatid')
	#TODO
]

tfalc_hn_arr = [
	 ('TSTF', 'S20'), ('RJD', float), ('IM1', float), ('IE1', float), 
	 ('IQ1', 'S2'),   ('IM2', float), ('IE2', float), ('IQ2', 'S2' ), 
	 ('IM3', float),  ('IE3', float), ('IQ3', 'S2'),  ('RM1' , float), 
	 ('RM2', float),  ('RM3', float), ('EP1', float), ('EP2', float), 
	 ('EP3', float) , ('TF1' , float), ('TF2', float), ('TF3', float)
   ]
tfalc_hn_dt = np.dtype(tfalc_hn_arr)

tfalc_hs_arr = [
	('HAT','S15'), ('TSTFC','S10'), ('FLD','S15'), ('BJD',float), ('IM1',float), 
	('IE1',float), ('IQ1','S1' ), ('IM2'  ,float), ('IE2',float), ('IQ2','S1' ), 
	('IM3',float), ('IE3',float), ('IQ3','S1'), ('RM1',float), ('RM2'  ,float), 
	('RM3',float), ('EP1',float), ('EP2',float), ('EP3',float), ('TF1',float), 
	('TF2'  ,float), ('TF3',float), ('XCC',float), ('YCC',float),  ('BGV',float), 
	('BGE'  ,float), ('FSV',float), ('FDV',float), ('FKV',float),
	('IHA',float), ('IZD'  ,float), ('RJD',float)
]
tfalc_hs_dt = np.dtype(tfalc_hs_arr)

##################################################################
#				   miscellanious utilities					     #
##################################################################
convert_ndarray_to_dict = \
		lambda x, index_column : { 
			row[index_column] : { c : row[c] for c in x.dtype.names \
										if not c == index_column } 
		}

def clean_nulls(filename):
	""" removes lines with NULL values
		and returns the full text
	"""
	cleaned = ""
	with open(filename, 'r') as f:
		for line in f:
			if 'NULL' in line: continue
			cleaned = "%s%s\n"%(cleaned, line)
	return cleaned

def make_dict_from_file(field, load_raw, save_dict, index_column, **kwargs):
	""" loads a file as a numpy ndarray, and then converts this
		to a python dict indexed by 'indexed_column' and containing
		all other entries (except those specified by 'index_column')
	"""
	data, msg = load_raw(field, **kwargs)
	if not msg is OK: return None, msg

	data_dict = convert_ndarray_to_dict(data, index_column)

	save_dict(field, data_dict, **kwargs)

	return data_dict, OK

def load_raw_file(field, get_filename, dtype, errmsg, **kwargs):
	""" loads a file as a numpy ndarray with the specified dtype
		will return None, errmsg if the file does not exist.
	"""
	filename = get_filename(field)
	if not os.path.exists(filename): return None, errmsg

	data = np.loadtxt(get_filename(field, **kwargs), dtype=dtype)

	return data, OK

def fetch_processed_file(field, get_filename, errmsg, **kwargs):
	""" loads the python dictionary associated with a processed
		file (e.g. twomass, keylist, etc.)
	"""
	dict_filename = get_filename(field, 'dict')
	if not os.path.exists(dict_filename): return None, errmsg

	return pickle.load(open(dict_filename, 'rb'))

def add_settings_to_parser(parser, **kwargs):
	""" adds arguments for each setting in kwargs
		to an argparse.ArgumentParser object
	"""
	for variable, argparse_kwargs in kwargs.iteritems():
		parser.add_argument("--%s"%(variable), **argparse_kwargs)
		
def load_ascii_list(filename, comment='#'):
	with open(filename, 'r') as f:
		txt = f.read()
		lines = [ [ c.split() for c in line.split(comment) ][0] for line in lines ]
		lines = [ line for line in lines if len(line) > 0 ]
		if len(lines[0]) == 1:
			lines = [ line[0] for line in lines ]
		return lines
		
def add_hatids_to_dict(d, hatids):
	for hatid in hatids:
		field = field_of(hatid)
		if not field in d: d[field] = set([hatid])
		else: d[field].add(hatid)

config_translation = {
	'True' : True,
	'False' : False,
	'None' : None
}
def translate_configparser(parser, section):
	
	translated = {}
	for key, value in parser.items(section):
		if value in config_translation:
			translated[key] = config_translation[value]
		elif ',' in value:
			translated[key] = value.split(',')
		else:
			translated[key] = value
	for key, value in translated.iteritems():
		logger.info("Read {0} as {1} from config file.".format(key, value))
	return translated

##################################################################
# FILENAMES AND OTHER FUNCTIONS
##################################################################

raw_lc_open_func = gzip.open
processed_lc_open_func = gzip.open

save_processed_lc = lambda lc, outfile : pickle.dump(lc, processed_lc_open_func(outfile, 'wb'))

get_missing_exposures_filename = lambda field, **kwargs      	: os.path.join(kwargs['SCRATCH'], 'INFO', 'KEYLISTS', 'missing_exposures_%s.set'%(field))
get_hatid_list_filename        = lambda field, **kwargs      	: os.path.join(kwargs['SCRATCH'], 'INFO', 'HATID_LISTS', "hatids_in_field_%s.list"%(field))
get_keylist_filename           = lambda field, ext  , **kwargs 	: os.path.join(kwargs['SCRATCH'], 'INFO', 'KEYLISTS', "keylist_%s.%s"%(field, ext))
get_lcstat_filename            = lambda field, ext  , **kwargs 	: os.path.join(kwargs['SCRATCH'], 'INFO', 'LCSTAT', 'lcstat_%s.%s'%(field, ext))
get_twomass_filename           = lambda field, ext  , **kwargs 	: os.path.join(kwargs['SCRATCH'], 'INFO', 'TWOMASS',  "twomass_%s.%s"%(field, ext))
get_raw_lc_filename            = lambda lcdir, hatid, **kwargs 	: os.path.join(lcdir, "%s.tfalc.gz"%(hatid))
get_processed_lc_filename      = lambda lcdir, hatid, **kwargs 	: os.path.join(lcdir, "%s-processed.pkl.gz"%(hatid))
get_lc_dir                     = lambda field, **kwargs 		: os.path.join(kwargs['SCRATCH'], 'LCCACHE', 'FIELDS', field)

load_raw_twomass 		= lambda field, **kwargs : load_raw_file(field, get_twomass_filename, twomass_dt, NO_RAW_TWOMASS_FILE_FOUND, **kwargs)
load_raw_lcstat 		= lambda field, **kwargs : load_raw_file(field, get_lcstat_filename , keylist_dt, NO_RAW_LCSTAT_FILE_FOUND , **kwargs)

save_twomass_dict       = lambda field, twomass_dict, **kwargs : pickle.dump(twomass_dict, open( get_twomass_filename( field, 'dict', **kwargs ), 'wb' ))
save_keylist_dict       = lambda field, keylist_dict, **kwargs : pickle.dump(keylist_dict, open( get_keylist_filename( field, 'dict', **kwargs ), 'wb' ))
save_lcstat_dict        = lambda field, lcstat_dict , **kwargs : pickle.dump(lcstat_dict , open( get_lcstat_filename(  field, 'dict', **kwargs ), 'wb' ))

make_twomass_dict		= lambda field, **kwargs : make_dict_from_file(field, load_raw_twomass, save_twomass_dict, 'hatid', **kwargs)
make_keylist_dict 		= lambda field, **kwargs : make_dict_from_file(field, load_raw_keylist, save_keylist_dict, 'hatid', **kwargs)
make_lcstat_dict  		= lambda field, **kwargs : make_dict_from_file(field, load_raw_lcstat , save_lcstat_dict , 'hatid', **kwargs)

keylist_dict_filename 	= lambda field, **kwargs : get_keylist_filename(field, 'dict', **kwargs )
twomass_dict_filename 	= lambda field, **kwargs : get_twomass_filename(field, 'dict', **kwargs )
lcstat_dict_filename 	= lambda field, **kwargs : get_lcstat_filename( field, 'dict', **kwargs )

fetch_hatids 			= lambda field, **kwargs : fetch_processed_file(field, get_hatid_list_filename       , NO_HATID_LIST_FILE_FOUND        , **kwargs)
fetch_keylist 			= lambda field, **kwargs : fetch_processed_file(field, keylist_dict_filename         , NO_KEYLIST_DICTIONARY_FILE_FOUND, **kwargs)
fetch_twomass 			= lambda field, **kwargs : fetch_processed_file(field, twomass_dict_filename         , NO_TWOMASS_DICTIONARY_FILE_FOUND, **kwargs)
fetch_lcstat 			= lambda field, **kwargs : fetch_processed_file(field, lcstat_dict_filename          , NO_LCSTAT_DICTIONARY_FILE_FOUND , **kwargs)
fetch_missing_exposures = lambda field, **kwargs : fetch_processed_file(field, get_missing_exposures_filename, NO_MISSING_EXPOSURES_FILE_FOUND , **kwargs)

get_eras_dir = lambda **kwargs : os.path.join(kwargs['SCRATCH'], 'LCCACHE', 'ERAS')
get_fields_dir = lambda **kwargs : os.path.join(kwargs['SCRATCH'], 'LCCACHE', 'FIELDS')

get_era_contents_filename = lambda **kwargs : os.path.join(kwargs['SCRATCH'], 'LCCACHE', 'era_contents.dict')
get_field_contents_filename = lambda **kwargs : os.path.join(kwargs['SCRATCH'], 'LCCACHE', 'field_contents.dict')

make_hatid_from_filename  =lambda fname, **kwargs : fname.replace('.tfalc.gz', '')
is_a_lightcurve = lambda fname, **kwargs : '.tfalc.gz' in fname and 'HAT' in fname and len(fname) == ''
field_of = lambda hatid : hatid[4:7]

get_lc_filename_in_era = lambda hatid, era, **kwargs : os.path.join(get_eras_dir(**kwargs), era, os.path.basename(get_raw_lc_filename(hatid, **kwargs)))
get_lc_locations_filename = lambda **kwargs : os.path.join(kwargs['SCRATCH'], 'LCCACHE', 'lc_locations.dict')
get_lc_filename_in_field = lambda hatid, era, **kwargs : os.path.join(kwargs['SCRATCH'], 'LCCACHE', 'FIELDS', "%s.%s.tfalc.gz"%(hatid, era))

get_hatids_in_directory = lambda directory, **kwargs : set([ make_hatid_from_fname(f) for f in directory if is_a_lighcurve(f) in f ])


#'time_coordinate' : ('RJD', 'Column to use for exposure dates. Should always be RJD as this is common across HATNet/HATSouth LC formats'),
#	'flux_column' : ('TF1', 'Aperature to use for flux'