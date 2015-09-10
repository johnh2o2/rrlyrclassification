from pyvislc.vislc import Visualizer
import utils.readhatlc as rhlc
from Tkinter import Tk
import numpy as np
import os, sys, gzip
import utils.miscutils as mutils
import utils.featureutils as futils
import settings
import cPickle as pickle

# What iteration are we on?
iteration = mutils.get_iteration_number()

# Get list of candidate HAT ID's:
candidate_fname = settings.get_candidate_fname(iteration)
candidate_results_fname = settings.get_candidate_results_fname(iteration)

# Load candidate ID's
candidate_ids = pickle.load(open(candidate_fname, 'rb'))

# Prune out the ones that are already labeled
labeled_ids, categories = futils.LoadLabeledHatIDs()
candidate_ids = [ ID for ID in candidate_ids if not ID in labeled_ids ]

# Quit if there's nothing to do!
if len(candidate_ids) == 0:
	print "No ID's need to be labeled!"
	sys.exit()

# Fetch filenames of the gzipped csv lightcurves
candidate_filenames = [ mutils.get_lc_fname(ID) for ID in candidate_ids ]

# Add a 'hatid' if there isn't already one.
for fname,ID in zip(candidate_filenames, candidate_ids):
	lc = pickle.load(gzip.open(fname, 'rb'))
	if not 'hatid' in lc:
		lc['hatid'] = ID
	pickle.dump(lc, gzip.open(fname, 'wb'))


# Find the files that aren't in the data directory yet:
missing_ids = [ ID for fname,ID in zip(candidate_filenames, candidate_ids) if not os.path.exists(fname) ]
if len(missing_ids) > 0:
	print "%d missing ids:"%(len(missing_ids))
	print missing_ids

# Fetch the LC's from the HAT (HTTP) server
#mutils.fetch_lcs(missing_ids)

# Labels that the user can give to each lightcurve
flags = [ 'RRab', 'Not-variable', 'Variable-not-RRab' ] 
shortcuts = [ 'q','w','e' ] # hit these keys on the keyboard to set each of the labels...

# Start the visualizer
root = Tk()
root.geometry('+250+80') 
root.wm_title("Label RR Lyrae candidates")
app = Visualizer(root, candidate_filenames, 
					logfile=candidate_results_fname, flag_shortcuts=shortcuts, flags=flags, jah=True)
root.mainloop()


# Now read in and interpret the results
mutils.logprint(" label_candidates : Now reading in the results")
f = open(candidate_results_fname, 'r')
classification_results = {}
while True:
	line = f.readline()
	if line == '': break
	line.replace('\n', '')

	splits = line.split(' ')
	assert(len(splits) <= 2)
	ID = splits[0]
	if len(splits) > 1: 
		CLASS = splits[1]
		if CLASS != 'RRab': CLASS = "none"
		else: CLASS = "RRAB"
	else: CLASS = None

	classification_results[ID] = CLASS
f.close()
mutils.logprint(" label_candidates : done.")

# Now save these into the global "labeled HatIDs" file.
mutils.logprint(" label_candidates : Now saving the results!")
for ID in classification_results:
	CLASS = classification_results[ID]
	if CLASS is None: continue

	if ID in labeled_ids:
		if CLASS != categories[ID]:
			print "WARNING: ID '%s' is already labeled as '%s', but the labeler says it's a '%s'. The labeler will take precedence here!"%(ID, categories[ID], CLASS)
			categories[ID] = CLASS

futils.SaveLabeledHatIDs(categories, iteration)
