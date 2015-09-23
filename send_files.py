import utils.readhatlc as rhlc
import numpy as np
import os, sys, gzip
import utils.miscutils as mutils
import utils.featureutils as futils
import settings
import cPickle as pickle
import argparse

parser = argparse.ArgumentParser(description='Send candidates')

parser.add_argument('--iteration', 
				dest 	= 'iteration', 
				default	= mutils.get_iteration_number(),
				type    = int,
				help	= 'Iteration to visualize/label')

parser.add_argument('--candidate-file',
				help	= 'File containing list of candidates')

args = parser.parse_args()

iteration = args.iteration

# Get list of candidate HAT ID's:
if args.candidate_file is None:
	candidate_fname = settings.get_candidate_fname(iteration)
else:
	candidate_fname = args.candidate_file

SMALL = 1E-5
show_fit = True

# Load candidate ID's
candidate_ids = pickle.load(open(candidate_fname, 'rb'))

# Prune out the ones that are already labeled
labeled_ids, categories = futils.LoadLabeledHatIDs()
candidate_ids = [ ID for ID in candidate_ids if not ID in labeled_ids ]

# Fetch lightcurve filenames
candidate_filenames = [ mutils.get_lc_fname(ID) for ID in candidate_ids ]

# Fetch feature filenames
candidate_feature_filenames = [ settings.hat_features_fname(ID, settings.model_prefix) for ID in candidate_ids ]

# Make directories
candidate_dir = "%s/candidates_%s"%(settings.SCRATCH, iteration)
lcdir = "%s/lc"%(candidate_dir)
fdir = "%s/features/"%(candidate_dir)
for dname in [ candidate_dir, lcdir, fdir ]:
	os.makedirs(dname)

print "Copying lightcurves."
for fname in candidate_filenames:
	os.system("cp %s %s/"%(fname, lcdir))
print "Copying features."
for fname in candidate_feature_filenames:
	os.system("cp %s %s/"%(fname, fdir))
print "Archiving."
os.system("tar cvf %s.tar %s/"%(candidate_dir, candidate_dir))

print "Done: %s.tar"%(candidate_dir)
