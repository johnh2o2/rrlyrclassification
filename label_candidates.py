from pyvislc.vislc import Visualizer
import utils.readhatlc as rhlc
from Tkinter import Tk
import numpy as np
import os, sys, gzip
import utils.miscutils as mutils
import utils.featureutils as futils
import settings
import cPickle as pickle
import argparse

parser = argparse.ArgumentParser(description='Visualize/label candidates')

parser.add_argument('--iteration', 
				dest 	= 'iteration', 
				default	= mutils.get_iteration_number(),
				type    = int,
				help	= 'Iteration to visualize/label')

parser.add_argument('--candidate-file',
				action 	= 'store_true',
				help	= 'File containing list of candidates')
	)

parser.add_argument('--include-labeled', 
				action  = 'store_true',
				help	= 'Dont skip already-labeled candidates')
	)

parser.add_argument('--output-file', 
				dest 	= 'outputfile', 
				default = None,
				help	= 'File to save the results of labeling')

parser.add_argument('--visualize', 
				dest 	= 'visualize', 
				action 	= 'store_true',
				help 	= 'Use pyvislc to visualize/label')

parser.add_argument('--dont-save', 
				dest 	= 'save', 
				action 	= 'store_false', 
				help 	= 'Dont save the results!')

parser.add_argument('--save', 
				dest 	= 'save', 
				action 	= 'store_true', 
				help 	= 'Save the results.')

args = parser.parse_args()

visualize = args.visualize
save = args.save
iteration = args.iteration
include_labeled = args.include_labeled

# Where do we store the results of the labeling?
if args.outputfile is None and save:
	candidate_results_fname = settings.get_candidate_results_fname(iteration)
elif save:
	candidate_results_fname = args.outputfile
else:
	candidate_results_fname = "label_candidates.results"

# Get list of candidate HAT ID's:
if args.candidatefile is None:
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

# Quit if there's nothing to do!
if len(candidate_ids) == 0 and not include_labeled:
	print "No ID's need to be labeled!"
	sys.exit()

# Fetch lightcurve filenames
candidate_filenames = [ mutils.get_lc_fname(ID) for ID in candidate_ids ]

# Find the files that aren't in the data directory yet:
missing_ids = [ ID for fname,ID in zip(candidate_filenames, candidate_ids) if not os.path.exists(fname) ]
if len(missing_ids) > 0:
	print "%d missing ids:"%(len(missing_ids))
	print missing_ids

All_Features = LoadAllFeatures([ ID for ID in candidate_ids if not ID in missing_ids ])

def ax_candidate_pf(ax, t, y, opts):
    """ Adds a phase-folded plot to a matplotlib.Axis instance

        inputs:
            ax          matplotlib.Axis instance to which phase-folded plot will be added
            t           Central values of phase bins
            y           Average magnitude in each phase bin
            opts        Options (required, should contain 'ylabel', 'period', 'yerr' )

    """
    ax.scatter(t , y, label="P=%.4f days"%(opts['period']),c='b'
        ,alpha=0.05, marker=',')
    V = opts['visualizer']
    phase = V.phase_offset
    HATID = V.lc['hatid']
    #score = scores[HATID]
    features = All_Features[HATID]
    period = features['p1']

    ax.scatter(t , y, label="P=%.4f days (Pbest=%.4f)"%(opts['period'], period),c='b'
        ,alpha=0.05, marker=',')

    amplitudes, phases, C = translate_features_to_fit_params(features)
    fitfunc = lambda T : C + sum( [ A * cos(2*np.pi*N*t - PHI) for A,N,PHI in zip(amplitudes, np.arange(1, len(amplitudes) + 1) , phases) ] )
    

    if show_fit and 2 * abs(period - opts['period'])/(period + opts['period']) < SMALL:
    	tfit = np.linspace(0, 1)
    	yfit = np.array( [ fitfunc( T  - phase) for T in tfit ])
    	ax.plot( tfit, yfit, lw=2, color='r')

    ax.set_ylabel(opts['ylabel'])
    ax.set_xlabel("Phase")
    ax.set_title("Phase-folded")
    ax.set_xlim(0,1)
    ax.legend(loc='upper right', fontsize=12)
    ax.invert_yaxis()

    ax.format_coord = lambda x, y : ''

if visualize:
	# Labels that the user can give to each lightcurve
	flags = [ 'RRab', 'Not-variable', 'Variable-not-RRab' ] 
	shortcuts = { 'q' : 'RRab', 'w' : 'Not-variable', 'e' : 'Variable-not-RRab' } # hit these keys on the keyboard to set each of the labels...

	# Start the visualizer
	root = Tk()
	root.geometry('+250+80') 
	root.wm_title("Label RR Lyrae candidates")
	app = Visualizer(root, candidate_filenames, 
						logfile=candidate_results_fname, flag_shortcuts=shortcuts, flags=flags, jah=True)

	app.phase_plot_axis_function = ax_candidate_pf
	app.set_phase_folded()
	root.mainloop()

# Now read in and interpret the results
mutils.logprint(" label_candidates : Now reading in the results")
f = open(candidate_results_fname, 'r')
classification_results = {}
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
		if CLASS != 'RRab': CLASS = "none"
		else: CLASS = "RRAB"
	else: CLASS = None

	classification_results[ID] = CLASS
f.close()
mutils.logprint(" label_candidates : done.")

# Now save these into the global "labeled HatIDs" file.
if save:
	mutils.logprint(" label_candidates : Now saving the results!")
	for ID in classification_results:
		CLASS = classification_results[ID]

		if CLASS is None: 
			continue
		if not ID in categories: 
			categories[ID] = CLASS
			continue
		if CLASS != categories[ID]:
			print "WARNING: ID '%s' is already labeled as '%s', but the labeler says it's a '%s'. The labeler will take precedence here!"%(ID, categories[ID], CLASS)
			categories[ID] = CLASS

	futils.SaveLabeledHatIDs(categories, iteration)
