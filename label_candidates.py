import settings
if not settings.RUNNING_ON_DELLA:
	from pyvislc.vislc import Visualizer
	from Tkinter import Tk
import numpy as np
import utils.readhatlc as rhlc
import os, sys, gzip
import utils.miscutils as mutils
import utils.featureutils as futils
import cPickle as pickle
import argparse
from math import *


SMALL = 1E-5
show_fit = True

parser = argparse.ArgumentParser(description='Visualize/label candidates')

parser.add_argument('--iteration',
				help	= 'Current iteration')

parser.add_argument('--include-labeled', 
				action  = 'store_true',
				help	= 'Dont skip already-labeled candidates')
	

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

parser.add_argument('--results', 
				help 	= 'Results file; this will be used to add to labeled data.')

parser.add_argument('--tarfile',
				help    = 'Tar file containing features, scores and lightcurves from candidates')

args = parser.parse_args()

tarfile = args.tarfile
visualize = args.visualize
save = args.save
include_labeled = args.include_labeled
results = args.results
iteration = args.iteration

if results is None and tarfile is None:
	raise Exception(" Need to specify either results or tarfile")

# Where do we store the results of the labeling?
if args.outputfile is None:
	candidate_results_fname = "candidate_results.dat"
else:
	candidate_results_fname = args.outputfile

if settings.RUNNING_ON_DELLA and not results is None:
	candidate_results_fname = results


def unpack(tarfile):
	tarfile_rootname = tarfile.split("/")[-1]
	os.system("tar xvf %s -C %s"%(tarfile, settings.candidates_dir))
	os.system("cp %s %s/"%(tarfile, settings.candidates_dir))

	cparent_dir = "%s/%s"%(settings.candidates_dir,tarfile_rootname.split('.')[-2])
	lcdir = "%s/lc"%(cparent_dir)
	features_dir = "%s/features"%(cparent_dir)
	scores_dir = "%s/scores"%(cparent_dir)
	return cparent_dir, lcdir, features_dir, scores_dir

if not tarfile is None:
	# Unpack the tarfile
	pdir, lcdir, fdir, sdir = unpack(tarfile)

	# Load candidate ID's
	candidate_ids = pickle.load(open("%s/candidate_ids.pkl"%(pdir), 'rb'))

	# Prune out the ones that are already labeled (if we want to)
	try:
		labeled_ids, categories = futils.LoadLabeledHatIDs()
		candidate_ids = [ ID for ID in candidate_ids if include_labeled or not ID in labeled_ids ]
	except:
		print "Warning: Could not open labeled hatids!"

	# Quit if there's nothing to do!
	if len(candidate_ids) == 0:
		print "No ID's need to be labeled!"
		sys.exit()

	# Fetch lightcurve filenames
	candidate_filenames = [ "%s/%s-full.tfalc.gz"%(lcdir, ID) for ID in candidate_ids ]
	candidate_feature_filenames = [ "%s/%s-feats.pkl"%(fdir, ID) for ID in candidate_ids ]
	candidate_score_filenames = [ "%s/%s.scores"%(sdir, ID) for ID in candidate_ids ]

	All_Features = { ID: pickle.load(open(feat_fname, 'rb')) for ID, feat_fname in zip(candidate_ids, candidate_feature_filenames) }

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

    amplitudes, phases, C = mutils.translate_features_to_fit_params(features)

    fitfunc = lambda T : C + sum( [ A * cos(2*np.pi*N*T - PHI) for A,N,PHI in zip(amplitudes, np.arange(1, len(amplitudes) + 1) , phases) ] )
   
    if show_fit and 2 * abs(period - opts['period'])/(period + opts['period']) < SMALL:
    	tfit = np.linspace(0, 1)
    	yfit = np.array( [ fitfunc( T  - phase*2*np.pi) for T in tfit ])
    	ax.plot( tfit, yfit, lw=2, color='r')

    ax.set_ylabel(opts['ylabel'])
    ax.set_xlabel("Phase")
    ax.set_title("Phase-folded")
    ax.set_xlim(0,1)
    ax.legend(loc='upper right', fontsize=12)
    ax.invert_yaxis()

    ax.format_coord = lambda x, y : ''


if visualize:
	
	class CustomVisualizer(Visualizer):
		def __init__(self, *args, **kwargs):
			super( CustomVisualizer, self).__init__( *args, **kwargs)
			self.phase_plot_axis_function = ax_candidate_pf
			self.set_phase_folded()
	# Labels that the user can give to each lightcurve
	flags = [ 'RRab', 'Not-variable', 'Variable-not-RRab', 'Possible-RRab' ] 
	shortcuts = { 'q' : 'RRab', 'w' : 'Possible-RRab', 'e' : 'Variable-not-RRab', 'r' : 'Not-variable' } # hit these keys on the keyboard to set each of the labels...

	# Start the visualizer
	root = Tk()
	root.geometry('+250+80') 
	root.wm_title("Label RR Lyrae candidates")
	app = CustomVisualizer(root, candidate_filenames, 
						logfile=candidate_results_fname, flag_shortcuts=shortcuts, flags=flags, jah=True)

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
		if Class == 'Possible-RRab': continue # Don't label things unless you're sure!
		if CLASS != 'RRab': CLASS = "none"
		else: CLASS = "RRAB"
	else: CLASS = None

	classification_results[ID] = CLASS
f.close()
mutils.logprint(" label_candidates : done.")

# Now save these into the global "labeled HatIDs" file.
if save:
	if iteration is None:
		iteration = mutils.get_iteration_number()

	labeled_ids, categories = futils.LoadLabeledHatIDs()
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
