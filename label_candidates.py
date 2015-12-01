import settings
if not settings.RUNNING_ON_DELLA:
	from pylabeler.labeler import Labeler 
	from pyvislc.vislc import Visualizer
	from pyvislc.utils import phase_fold
	from Tkinter import Tk

import numpy as np
import utils.readhatlc as rhlc
import os, sys, gzip
import utils.miscutils as mutils
import matplotlib.pyplot as plt
import utils.featureutils as futils
import cPickle as pickle
import argparse
from math import *


SMALL = 1E-4
show_fit = True

parser = argparse.ArgumentParser(description='Visualize/label candidates')

parser.add_argument('--iteration',
				help	= 'Current iteration')

parser.add_argument('--remove-labeled', 
				action  = 'store_true',
				help	= 'Skip already-labeled candidates')
	
parser.add_argument('--fast',
				action  = 'store_true',
				help    = 'Uses lo-res png figures; much faster but no interactivity')

parser.add_argument('--remake',
				action  = 'store_true',
				help    = 'regenerates the .png figures for each hatid')

parser.add_argument('--label-file', 
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

parser.add_argument('--directory', 
				help 	= 'parent directory of files.')

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
fast_visualize = args.fast
include_labeled = (not args.remove_labeled)
results = args.results
iteration = args.iteration
REMAKE_FIGURES = args.remake

if results is None and tarfile is None and args.directory is None:
	raise Exception(" Need to specify either results or tarfile")

# Where do we store the results of the labeling?
if args.outputfile is None:
	candidate_results_fname = "candidate_results.dat"
else:
	candidate_results_fname = args.outputfile

if settings.RUNNING_ON_DELLA and not results is None:
	candidate_results_fname = results

def set_dirs(cparent_dir):
	lcdir = "%s/lc"%(cparent_dir)
	features_dir = "%s/features"%(cparent_dir)
	scores_dir = "%s/scores"%(cparent_dir)
	return cparent_dir, lcdir, features_dir, scores_dir


def unpack(tarfile):
	tarfile_rootname = tarfile.split("/")[-1]

	os.system("tar xvf %s -C %s"%(tarfile, settings.candidates_dir))
	os.system("cp %s %s/"%(tarfile, settings.candidates_dir))

	cparent_dir = "%s/%s"%(settings.candidates_dir,tarfile_rootname.split('.')[-2])
	return set_dirs(cparent_dir)

if not tarfile is None or not args.directory is None:
	if not tarfile is None:
		# Unpack the tarfile
		pdir, lcdir, fdir, sdir = unpack(tarfile)
	else:
		pdir, lcdir, fdir, sdir = set_dirs(args.directory)

	# Load candidate ID's
	candidate_ids = pickle.load(open("%s/candidate_ids.pkl"%(pdir), 'rb'))
	print "%d hatids found"%(len(candidate_ids))

	# Prune out the ones that are already labeled (if we want to)
	try:
		labeled_ids, categories = futils.LoadLabeledHatIDs()
		if not include_labeled: 
			candidate_ids = [ ID for ID in candidate_ids if not ID in labeled_ids ]
			print "  -> %d IDs (pruned out labeled IDs)"%(len(candidate_ids))

	except:
		print "Warning: Could not open labeled hatids!"

	# Quit if there's nothing to do!
	if len(candidate_ids) == 0:
		print "No ID's need to be labeled!"
		sys.exit()

	# Fetch lightcurve filenames
	candidate_filenames = [ "%s/%s-full.tfalc.gz"%(lcdir, ID) for ID in candidate_ids ]
	candidate_score_filenames = [ "%s/%s.scores"%(sdir, ID) for ID in candidate_ids ]
	candidate_feature_filenames = [ "%s/%s-feats.pkl"%(fdir, ID) for ID in candidate_ids ]

	# Prune out IDs that have a None lightcurve
	candidate_ids, candidate_filenames, candidate_feature_filenames = \
		zip(*[ (hatid,lcfname, ftfname) for hatid, lcfname, ftfname in zip(candidate_ids, candidate_filenames, candidate_feature_filenames) \
													if os.path.exists(lcfname) and os.path.exists(ftfname)\
													and not pickle.load(gzip.open(lcfname,'rb')) is None \
													and not pickle.load(open(ftfname, 'rb')) is None ])
	print "  -> %d IDs (pruned out None lcs and features)"%(len(candidate_ids))
	candidate_ids = list(candidate_ids)
	candidate_filenames = list(candidate_filenames)
	candidate_feature_filenames = list(candidate_feature_filenames)
	
	All_Features = { ID: pickle.load(open(feat_fname, 'rb')) for ID, feat_fname in zip(candidate_ids, candidate_feature_filenames) if os.path.exists(feat_fname) }




def ax_candidate_pf(ax, t, y, opts):
	""" Adds a phase-folded plot to a matplotlib.Axis instance

	    inputs:
	        ax          matplotlib.Axis instance to which phase-folded plot will be added
	        t           Central values of phase bins
	        y           Average magnitude in each phase bin
	        opts        Options (required, should contain 'ylabel', 'period', 'yerr' )

	"""
	inds = np.arange(len(t))
	np.random.shuffle(inds)
	add_ones = inds[:len(inds)/2]
	T2 = [ T + 1. if i in add_ones else T for i,T in enumerate(t)  ]

	#ax.scatter(T2 , y, label="P=%.4f days"%(opts['period']),c='b'
	#    ,alpha=0.05, marker=',')

	V = opts['visualizer']
	phase = V.phase_offset
	HATID = V.lc['hatid']
	#score = scores[HATID]
	features = All_Features[HATID]
	period = features['p1']

	if period != opts['period']:
		V.set_phase_folded(period=period)
		return

	ax.scatter(T2 , y, label="P=%.4f days (Pbest=%.4f)"%(opts['period'], period),c='b',alpha=0.05, marker=',')

	amplitudes, phases, C = mutils.translate_features_to_fit_params(features)

	fitfunc = lambda T : C + sum( [ A * cos(2*np.pi*N*T - PHI) for A,N,PHI in zip(amplitudes, np.arange(1, len(amplitudes) + 1) , phases) ] )

	if show_fit and 2 * abs(period - opts['period'])/(period + opts['period']) < SMALL:
		tfit = np.linspace(0, 2)
		yfit = np.array( [ fitfunc( T - phase ) for T in tfit ])
		ax.plot( tfit, yfit, lw=2, color='r')

	ax.set_ylabel(opts['ylabel'])
	ax.set_xlabel("Phase")
	ax.set_title("Phase-folded")
	ax.set_xlim(0,2)
	ax.legend(loc='upper right', fontsize=9)
	

	ax.format_coord = lambda x, y : ''


def fast_ax_candidate_pf(ax,hatid,lc_fname, ft_fname, show_fit=True, phase=0.0):
	""" Adds a phase-folded plot to a matplotlib.Axis instance

	    inputs:
	        ax          matplotlib.Axis instance to which phase-folded plot will be added
	        t           Central values of phase bins
	        y           Average magnitude in each phase bin
	        opts        Options (required, should contain 'ylabel', 'period', 'yerr' )

	"""
	
	# Open lightcurve and get features.
	#lc_fname = candidate_filenames[candidate_ids.index(hatid)]
	if not os.path.exists(lc_fname):
		ax.text(0.5, 0.5, '%s not found'%(lc_fname), transform=ax.transAxes, ha='center')
		return False
	else:
		lc = pickle.load(gzip.open(lc_fname, 'rb'))
		if lc is None:
			ax.text(0.5, 0.5, 'Lightcurve is None'%(lc_fname), transform=ax.transAxes, ha='center')
			return False
	# Raw time/mags
	ytype = 'TF2'
	t_, y_ = lc['BJD'], lc[ytype]
	try:
		t_ = [ float(T) for T in t_ ]
		y_ = [ float(Y) for Y in y_ ]
	except ValueError:
		print "Warning: %s cannot convert t, y"%(hatid)
		return False

	# Get features/period
	features = pickle.load(open(ft_fname, 'rb'))
	period = features['p1']

	# Phase fold
	t, y, junk = phase_fold(t_, y_, period)

	# edit that to get 0, 2 phase
	inds = np.arange(len(t))
	np.random.shuffle(inds)
	add_ones = inds[:len(inds)/2]
	T2 = [ T + 1. if i in add_ones else T for i,T in enumerate(t)  ]
	
	# Scatter plot
	ax.scatter(T2 , y, label="P=%.4f days"%(period),c='b',alpha=0.05, marker=',')

	# Get info about the fit
	amplitudes, phases, C = mutils.translate_features_to_fit_params(features)

	fitfunc = lambda T : C + sum( [ A * cos(2*np.pi*N*T - PHI) for A,N,PHI in zip(amplitudes, np.arange(1, len(amplitudes) + 1) , phases) ] )

	if show_fit:
		tfit = np.linspace(0, 2)
		yfit = np.array( [ fitfunc( T - phase ) for T in tfit ])
		ax.plot( tfit, yfit, lw=2, color='r')

	ax.set_ylabel(ytype)
	ax.set_xlabel("Phase")
	ax.set_title("Phase-folded")
	ax.text(0.05, 0.95, '%s'%(hatid), transform=ax.transAxes, ha='left', va='top')
	ax.set_xlim(0,2)
	ax.legend(loc='upper right', fontsize=9)
	ax.invert_yaxis()
	

	ax.format_coord = lambda x, y : ''
	return True

fig_folder = os.path.join(settings.parent_dir, 'plots')
pf_fig_name = lambda hatid : '%s/%s_pf.png'%(fig_folder, hatid)
key_mapping = {
		 'q' : 'RRab', 'w' : 'RRab-RRc', 'e' : 'RRc', 'r' : 'Variable-not-RR', 't' : 'Unsure' , 'y' : 'Not-variable', 'p' : 'Problem'
	}
def fast_visualize(hatids):
	

	def label(hatid):
		f = plt.figure(1)
		ax = f.add_subplot(111)
		i = candidate_ids.index(hatid)
		lc_fname = candidate_filenames[i]
		feat_fname = candidate_feature_filenames[i]
		fast_ax_candidate_pf(ax, hatid, lc_fname, feat_fname)
		f.canvas.draw()
		f.savefig(pf_fig_name(hatid), dpi=100)
		plt.close()

	
	# generate lores images of all of them
	for hatid in hatids:
		if not os.path.exists(pf_fig_name(hatid)) or REMAKE_FIGURES:
			print "Making plots for %s..."%(hatid)
			label(hatid)
	
	image_files = { hatid : [ pf_fig_name(hatid) ] for hatid in hatids }
	label_file = candidate_results_fname

	# Now label away!!
	labeler = Labeler(image_files, label_file, key_mapping.values(), key_mapping)
	labeler.connect()
	plt.show(block=True)
		
if visualize and fast_visualize:
	fast_visualize(candidate_ids)

elif visualize:
	class CustomVisualizer(Visualizer):
		def __init__(self, *args, **kwargs):
			super( CustomVisualizer, self).__init__( *args, **kwargs)
			self.phase_plot_axis_function = ax_candidate_pf
			self.set_phase_folded()
	# Labels that the user can give to each lightcurve
	flags = key_mapping.values()
	shortcuts = key_mapping
	
	# Start the visualizer
	root = Tk()
	root.geometry('+250+80') 
	root.wm_title("Label RR Lyrae candidates")
	app = CustomVisualizer(root, candidate_filenames, 
						logfile=candidate_results_fname, flag_shortcuts=shortcuts, flags=flags, jah=True)

	root.mainloop()

ignore_types = [ 'RRab-RRc', 'RRc', 'Unsure', 'Problem' ]
# Now read in and interpret the results
mutils.logprint(" label_candidates : Now reading in the results")
if os.path.exists(candidate_results_fname):

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
			if CLASS in ignore_types: continue # Don't label things unless you're sure!
			if CLASS != 'RRab': CLASS = "none"
			else: CLASS = "RRAB"
		else: CLASS = None

		classification_results[ID] = CLASS
	f.close()
	mutils.logprint(" label_candidates : done.")
else:
	mutils.logprint(" label_candidates : No file %s (could just mean that you didn't label anything...)"%(candidate_results_fname))

# Now save these into the global "labeled HatIDs" file.
if save:
	if iteration is None:
		iteration = mutils.get_iteration_number()

	labeled_ids, categories = futils.LoadLabeledHatIDs()
	mutils.logprint(" label_candidates : Now saving the results!")
	for ID in classification_results:
		CLASS = classification_results[ID]

		if CLASS is None:# and ID in categories:
			#del categories[ID] 
			continue
		if not ID in categories: 
			categories[ID] = CLASS
			continue
		if CLASS != categories[ID]:
			print "WARNING: ID '%s' is already labeled as '%s', but the labeler says it's a '%s'. The labeler will take precedence here!"%(ID, categories[ID], CLASS)
			categories[ID] = CLASS

	futils.SaveLabeledHatIDs(categories, iteration)
