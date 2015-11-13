import numpy as np
from math import *
import os, sys
import miscutils as mutils
from settings import *
import cPickle as pickle 
from simulated_rrlyr import *

lsp_cutoff = 15.
Nperbin = 100
min_mag = 6.
max_mag = 16.
nbins = 20

m5lcdir = "/Users/jah5/Documents/Fall2014_Gaspar/rrlyr_classification/information/J_MNRAS_411_1744"

fake_lc_dir = "/Users/jah5/Documents/Fall2014_Gaspar/fake_lc"

# Read in M5 lightcurves and fit fourier components
print "Reading M5 RR Lyrae lightcurves .."
m5lcs = ParseM5Data(m5lcdir)
print "Calculating LSP .."
for lc in m5lcs:
	print "   > ", lc.info['VName']
	freq, lspower, nfreq, imax, prob = lsp(lc.times, lc.mags)
	lc.info['lsp'] = { 'periods' : np.power(freq[::-1], -1), 'powers' : lspower[::-1],  }

	best_period = 1./freq[imax]
	lc.best_period = best_period

	phases = [ t/best_period - floor(t/best_period) for t in lc.times ]
	lc.phases = phases

	period, amplitudes, phases, C =get_fourier_components(lc)
	args = convert_to_args(period, amplitudes, phases, C)

	lc.fourier_components = ( period, amplitudes, phases, C )

all_hatids = pickle.load(open(os.join(information_dir, 'hatid_field_list.pkl'), 'rb'))
np.random.shuffle(all_hatids)

magnitude_bin_edges = np.linspace(6, 16, nbins)
dm = magnitude_bin_edges[1] - magnitude_bin_edges[0]
magnitude_bins = [ [] for edge in magnitude_bin_edges ]

def update_bins(content, M, bins, bin_edges):
	pass
def not_enough_lcs_in_bin(M):
	pass
def get_magnitude(features):
	pass
def check(hatid):
	pass
def inject(lc, source_lc):

	# Need to figure out how to get the right colors into this...
	# How to get 2mass colors??
	# should we just ignore colors altogether??
	# AAAAAAAAAAAH
	
	pass

for hatid in all_hatids:
	if not check(hatid): continue
	features = get_features(hatid)
	M = get_magnitude(features)
	if not_enough_lcs_in_bin(M): 
		magnitude_bins = update_bins(get_lightcurve(hatid), M, magnitude_bins, magnitude_bin_edges)


for BIN in magnitude_bins:
	for LC in BIN:
		fake_lc = inject(LC, np.random.choice(m5lcs))


