import settings
import numpy as np
import utils.readhatlc as rhlc
import os, sys, gzip
import utils.miscutils as mutils
import utils.featureutils as futils
import cPickle as pickle
import argparse
from astroquery.simbad import Simbad
from math import *
import astropy.units as u
import astropy.coordinates as coord


SMALL = 1E-4
model = "rrab_v4"
iterations = 5
directory = settings.parent_dir
fnames = [ 
			( os.path.join(directory, 'candidate_results_%s_iter%04d.dat'%(model, i)), 
			  os.path.join(directory, 'candidates_%s_iter%04d.tar'%(model, i)) ) 
			for i in range(1, iterations + 1) 
		]

def unpack(tarfile):
	tarfile_rootname = tarfile.split("/")[-1]
	os.system("tar xf %s -C %s"%(tarfile, settings.candidates_dir))
	os.system("cp %s %s/"%(tarfile, settings.candidates_dir))

	cparent_dir = "%s/%s"%(settings.candidates_dir,tarfile_rootname.split('.')[-2])
	lcdir = "%s/lc"%(cparent_dir)
	features_dir = "%s/features"%(cparent_dir)
	scores_dir = "%s/scores"%(cparent_dir)
	return cparent_dir, lcdir, features_dir, scores_dir

def info(fts):
	return "%s %f %f %.4f"%(fts['hatid'], fts['ra'], fts['dec'], fts['V'])

nfound = 0
All_Features = {}
labels = {}
for results, tarfile in fnames:

	# Now read in and interpret the results
	mutils.logprint(" collate_discoveries : Now reading in the results")
	if not os.path.exists(results): 
		mutils.logprint(" collate_discoveries : Can't find results file '%s'"%(results))
		continue
	if not os.path.exists(tarfile):
		mutils.logprint(" collage_discoveries : Can't find tarfile '%s'"%(tarfile))
	

	f = open(results, 'r')
	
	while True:
		line = f.readline()
		if line == '': break
		line = line.replace('\n', '')

		splits = line.split(' ')
		assert(len(splits) <= 2)
		ID = splits[0]
		if len(splits) > 1: 
			CLASS = splits[1]
			print CLASS
			if CLASS == 'Possible-RRab': continue # Don't label things unless you're sure!
			if CLASS != 'RRab': CLASS = "none"
			else: CLASS = "RRAB"
		else: CLASS = None

		labels[ID] = CLASS
	f.close()

	mutils.logprint(" collate_discoveries : done.")

	# Unpack the tarfile
	pdir, lcdir, fdir, sdir = unpack(tarfile)

	# Load candidate ID's
	candidate_ids = pickle.load(open("%s/candidate_ids.pkl"%(pdir), 'rb'))

	# Fetch lightcurve filenames
	candidate_filenames = [ "%s/%s-full.tfalc.gz"%(lcdir, ID) for ID in candidate_ids ]
	candidate_feature_filenames = [ "%s/%s-feats.pkl"%(fdir, ID) for ID in candidate_ids ]
	candidate_score_filenames = [ "%s/%s.scores"%(sdir, ID) for ID in candidate_ids ]

	for ID, feat_fname in zip(candidate_ids, candidate_feature_filenames): All_Features[ID] = pickle.load(open(feat_fname, 'rb')) 

for ID in labels:
	if labels[ID] == 'RRAB':
		nfound += 1
		if not 'hatid' in All_Features[ID]:
			All_Features[ID]['hatid'] = ID
		print info(All_Features[ID])
print "----------------------------------------"
for ID in labels:
	if labels[ID] == 'RRAB':
		if not 'hatid' in All_Features[ID]:
			All_Features[ID]['hatid'] = ID
		f = All_Features[ID]
		print f['ra'], f['dec']
		c = coord.FK5(ra=f['ra'], dec=f['dec'], unit=(u.deg, u.deg))
		results = Simbad.query_region(c,radius=2 * u.arcminute,epoch='J2000', equinox=2000)
		print results.colnames

print " Done: %d found"%(nfound)
