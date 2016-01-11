import cPickle as pickle
import os, sys
import numpy as np
import logging
from config import *

logger = logging.getLogger(__name__)


def generate_contents_file(parent_dir, outfilename, **kwargs):
	""" generates dictionary of { folder : hatids } (used by 
		generate_era_contents_file and generate_field_contents_file)
		
		* saves dictionary to outfilename
		* returns said dictionary and a status message
	"""
	logger.info("looking for folders in parent_dir=%s"%(parent_dir))
	if not os.path.isdir(parent_dir):
		logger.error("cannot find parent directory %s"%(parent_dir))
		return None, CANNOT_FIND_PARENT_DIRECTORY

	logger.info("making a list of folders in %s"%(parent_dir))
	folders = [ folder for files in os.listdir(parent_dir) if os.path.isdir(folder) ]

	logger.info("%d folders found in %s"%(len(folders), parent_dir))
	logger.info("compiling contents of each folder in %s"%(parent_dir))
	contents = { folder : get_hatids_in_directory(os.path.join(parent_dir, folder)) for folder in folders }

	logger.info("saving contents to %s"%(outfilename))
	pickle.dump(contents, open(outfilename, 'wb'))

	return contents, OK

def generate_era_contents_file(**kwargs):
	logger.info("generating contents of era folders")
	return generate_contents_file(get_eras_dir(**kwargs), get_era_contents_filename(**kwargs), **kwargs)

def generate_field_contents_file(**kwargs):
	logger.info("generating contents of field folders")
	return generate_contents_file(get_fields_dir(**kwargs), get_fields_contents_filename(**kwargs), **kwargs)

def generate_hatid_lc_locations(era_contents, **kwargs):
	locations = { }
	for era, contents in era_contents.iteritems():
		logger.info("adding hatids from era %s"%(era))
		for hatid in contents:
			if not hatid in locations: locations[hatid] = set([ era ])
			else: locations[hatid].add(era)

	fname = get_lc_locations_filename(**kwargs)
	logger.info("saving lc_locations contents to %s"%(fname))
	pickle.dump(locations, open(fname, 'wb'))

	return locations, OK

def field_contents_from_hatid_lc_locations(hatid_lc_locations, **kwargs):
	field_contents = {}
	for hatid, locations in hatid_lc_locations.iteritems():
		field = field_of(hatid)
		if not field in field_contents: field_contents[field] = set([ hatid ])
		else: field_contents[field].add(hatid)

	return field_contents

def symlink_era_lc_to_field(hatid, era, **kwargs):
	field = field_of(hatid)

	# Filename of the lightcurve
	fname_era = get_lc_filename_in_era(hatid, era, **kwargs)

	# Filename of the symlink to the lightcurve
	fname_field = get_lc_filename_in_field(hatid, era, **kwargs)

	# If there's already a symlink and we don't want to redo things, continue
	if os.path.exists(fname_field) and not any([ kwargs[r] in kwargs for r in ['redo', 'redo_all' ] if r in kwargs ]): return OK

	# Make the symlink
	os.symlink(fname_era, fname_field)

	return OK

def symlink_era_lcs_to_field(hatid, locations, **kwargs):
	for era in locations:
		msg = symlink_era_lc_to_field(hatid, era, **kwargs)

	return OK

def create_symlinks_for_field(field, hatid_lc_locations, **kwargs):
	field_directory = get_lc_dir(field, **kwargs)
	if not os.path.isdir(field_directory): 
		os.makedirs(field_directory)

	for hatid, locations in hatid_lc_locations:

		if not field_of(hatid) == field: continue

		msg = symlink_era_lcs_to_field(hatid, locations, **kwargs)

		# If we're planning on combining all of the lightcurves, 
		# then we're done 
		if kwargs['deal_with_multiple_lcs']: continue
		
		# Otherwise (simple mode) we just use the lightcurve
		# with the largest filesize
		filenames = [ get_lc_filename_in_era(hatid, era, **kwargs) for era in locations ]
		sizes = [ os.path.sizeof(f) for f in filenames]

		biggest_lc_fname = filenames[sizes.index(max(sizes))]

		final_lc_fname = get_raw_lc_filename(field_directory, hatid)
		os.symlink(biggest_lc_fname, final_lc_fname)

	return OK

def create_symlinks_for_all_fields(field_contents, hatid_lc_locations, **kwargs):
	for field, contents in field_contents.iteritems():
		hidlocs = { hatid : hatid_lc_locations[hatid] for hatid in contents }
		msg = create_symlinks_for_field(field, hidlocs, **kwargs)
		if not msg is OK: return msg

	return OK

if __name__ == '__main__':
	
















