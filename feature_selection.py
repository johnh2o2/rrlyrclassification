import numpy as np 
from multiprocessing import Pool
from scipy.optimize import curve_fit
from scipy.stats import skew, kurtosis
from fastperiod import specwindow, lombscargle, dworetsky, stetson
import lsp, sys, os, re
from os.path import exists
import matplotlib.pyplot as plt
from math import *
from utils import *
from time import time
import readhatlc as rhlc
import cPickle as pickle
VERBOSE=False
feat_dir = "/Users/jah5/Documents/Fall2014_Gaspar/work/saved_features"
small = eps
FORCE_REDO = True
FORCE_IDLIST_REDO = False
NPROC=4
min_nobs = 20
# features: 
# ===========
#  	rms
#   chi2
#	kurtosis
#	skew
#	frac_beyond1std
#	frac_beyond2std
#	frac_beyond3std
#	max_lsp_power_all		 -- this should not be used in classification
#	best_ls_per_all
#
#	rms_resid
#	chi2_resid
#	kurtosis_resid
#	skew_resid
#	frac_beyond1std_resid
#	frac_beyond2std_resid
#	frac_beyond3std_resid
#	max_lsp_power_resid		 -- this should not be used in classification
#	best_ls_per_resid
#
#   COLORS: V,(R-V,I-V,J-V,H-V,K-V) <-- do LDA on flux ratios to get a composite color: A
#	gal_lat					> (we have to be careful with these -- strong sample bias?)
#	gal_long				|
#	closest_neighbor_dist	> To implement later. Shouldn't be extremely important
#	closest_neighbor_mag    |
#	closest_neighbor_flag	|
#	
#   for each period:
#   ================
#    max_lsp_power			
#    best_lsp_period		
#    min_stetson_len			using n_peaks from LSP (and searching around them) find minimum stetson len
#    stetson_period				corresponding period
#    max_dworetsky_len			
#    (dworetsky_period)
#    (min_resid_std)
#	 amplitude
#	 constant_offset
#    pct_of_rms_explained
#    (rms_of_residual)
#    chi2_pf
#
#   for each harmonic:
#   ------------------
#		amplitude 				relative to the total amplitude
#		phase 					relative to the 0th harmonic
#
def BinnedDistro(x, nbins=20, zrange=(-5, 5)):
	zmin, zmax = zrange
	Z = (x - np.mean(x))/np.std(x)
	dz = float(zmax - zmin)/float(nbins)
	edges = [ zmin + i*dz for i in range(nbins) ]
	edges.append(zmax)
	distro = np.zeros(nbins)
	for z in Z:
		bin = float(z - zmin)/dz
		if bin < 0: bin = 0
		if bin > nbins - 1: bin = nbins-1

		distro[bin] += 1.
	return edges, distro/sum(distro)
def ScaledMoment(x, N):
	sig = np.std(x)
	Z = (x - np.mean(x))/sig
	return np.mean(np.power(Z, N))/(sig**N)
def FeatureVectorFileName(hatid):
	return "%s/%s-features-v3.pkl"%(feat_dir, hatid)
def Stetson(t, x, p, err=None):
	n_obs = len(x)
	j_range = n_obs-1

	fold_t = np.min(t) # fold at the first time element

	phase = (t - fold_t)/p - np.floor((t - fold_t)/p)
	phase_sort_ind = np.argsort(phase)

	phase_sorted = phase[phase_sort_ind]
	x_sorted = x[phase_sort_ind]
	if err is None:
		err = get_errs(x_sorted)
	
	err_sorted = err[phase_sort_ind]

	strlen_numer = 0.0
	strlen_denom = 0.0

	# build up the string                 
	for j in range(j_range):
	    
	    weights = ( ( (err_sorted[j])*(err_sorted[j]) + 
	                  (err_sorted[j+1])*(err_sorted[j+1]) ) *
	                (phase_sorted[j+1] - phase_sorted[j] + 1.0/n_obs) )
	    weights = 1.0/weights

	    strlen_numer += weights*np.fabs(x_sorted[j] - x_sorted[j+1])
	    strlen_denom += weights

	# add the last element
	weights = ( ( (err_sorted[-1])*(err_sorted[-1]) + 
	              (err_sorted[0])*(err_sorted[0]) ) *
	            (phase_sorted[0] - phase_sorted[-1] + 1.0/n_obs) )
	weights = 1.0/weights

	strlen_numer += weights*np.fabs(x_sorted[-1] - x_sorted[0])
	strlen_denom += weights

	strlen = strlen_numer/strlen_denom

	return strlen
def Dworetsky(t, x, p):
	j_range = len(x)-1
	fold_t = np.min(t) # fold at the first time element

	mod_x = (x - np.min(x))/(2.0*(np.max(x) - np.min(x))) - 0.25

	phase = (t - fold_t)/p - np.floor((t - fold_t)/p)

	phase_sort_ind = np.argsort(phase)
	phase_sorted = phase[phase_sort_ind]
	mod_x_sorted = mod_x[phase_sort_ind]

	strlen = 0.0

	# now calculate the string length
	for j in range(j_range):
	    strlen += np.sqrt( (mod_x_sorted[j+1] - mod_x_sorted[j]) * 
	                       (mod_x_sorted[j+1] - mod_x_sorted[j]) +
	                       (phase_sorted[j+1] - phase_sorted[j]) * 
	                       (phase_sorted[j+1] - phase_sorted[j]))

	strlen += np.sqrt( (mod_x_sorted[0] - mod_x_sorted[-1]) * 
	                   (mod_x_sorted[0] - mod_x_sorted[-1]) +
	                   (phase_sorted[0] - phase_sorted[-1] + 1) *
	                   (phase_sorted[0] - phase_sorted[-1] + 1))

	return strlen   
def BinWidths(X):
	dX_lefts = {}
	dX_rights = {}
	dX_lefts[X[0]] = X[1] - X[0]
	dX_rights[X[0]] = X[1] - X[0]
	dX_lefts[X[-1]] = X[-1] - X[-2]
	dX_rights[X[-1]] = X[-1] - X[-2]
	for i in range(1, len(X)-1): 
		dX_lefts[X[i]] = X[i] - X[i-1]
		dX_rights[X[i]] = X[i+1] - X[i]
	return dX_lefts, dX_rights
def LombScargle(t, x ):
	wk1,wk2,nout,jmax,prob = lsp.fasper(t, x, ofac, hifac, MACC)
	P = np.power(wk1[::-1],-1)
	LSP = wk2[::-1]
	return P, LSP
def ClipPeriods(P, LSP, minperiod, maxperiod):
	okinds = [ i for i in range(len(P)) if P[i] < maxperiod and P[i] > minperiod ] 
	return P[okinds], LSP[okinds]
def WeedOutHarmonicPeaks(peaks, periods, threshold=0.01):
	new_peaks = []
	for per in peaks:
		is_ok = True
		for Per in periods:
			if abs(per/Per - int(round(per/Per))) < small or abs(Per/per - int(round(Per/per))) < small: 
				is_ok = False
				break
		if is_ok: new_peaks.append(per)
	return new_peaks
def GetSearchPeriods(t, x, periods):
	# Get L-S; weed out long periods
	P, LSP = LombScargle(t, x)
	P, LSP = ClipPeriods( P, LSP, 0., max_per )
	# get a list of the top N peaks in the LSP
	peaks,pows = find_n_peaks(P, LSP, n_peaks)
	# weed out peaks that are harmonics
	peaks = WeedOutHarmonicPeaks(peaks, periods)
	# get the bin widths for each period (equal spacing in w ~ 1/p)
	dP_lefts, dP_rights = BinWidths(P)
	# Get periods surrounding the LSP peaks
	search_pers = [ ]
	for per in peaks:
		for E in np.linspace(-delta_P*dP_lefts[per], delta_P*dP_rights[per], NSEARCH):
			search_pers.append( per+E )
	return search_pers
def GetModelLightcurve(t, ws, amps, phs, c):
	arr = np.array([ A*np.cos(W*t - P) for A, W, P in zip(amps, ws, phs) ])
	mod = np.zeros(len(arr))
	for i in range(len(mod)):
		mod[i] = sum(arr[:,i]) + c
	return mod
def GetAmplitude(*args):
	pers = 2*np.pi*np.power(args[0],-1)
	T = np.linspace(0,max(pers))
	X = GetModelLightcurve(T, *args)
	
	return 0.5*(max(X) - min(X))
def Whiten(t, x, ws, amps, phs, c):
	return t, get_resid( t, x, ws, amps, phs, c )
def FitPeriods( t, x, verbose=VERBOSE ):
	resid = np.zeros(len(x))
	resid[:] = x[:]
	std_total = np.std(x)
	PERIODS, chi2s, stds = [], [], []
	while len(PERIODS) < npers:
		# start timer
		if verbose: t0 = time()
		# Get subset of periods to search
		search_pers = GetSearchPeriods(t, resid, PERIODS)
		if verbose: print "     seaching..."
		# Find the period that minimizes the scatter in the residual
		best_pers, sig_pers = find_best_of_n_periods(t,resid,search_pers,other_periods=PERIODS)
		# End timer.
		if verbose: dt = time() - t0
		if verbose: print "     done (%f s)"%(dt)
		
		PERIODS.append(best_pers[0])

		# Now fit the period(s) to the data and update the residual.
		fit_params = fit_periods(t, x, PERIODS, use_dwdt=False)
		
		# Analyze the residual
		chi2s.append(get_chi2_pf(t,resid,PERIODS[-1]))

		std_old = np.std(resid)	
		resid[:] = get_resid(t, x, *fit_params)[:]
		std_new = np.std(resid)
		stds.append(std_new)
		
		

		# Print relevant information
		if verbose:
			# Now phase fold the residual around the highest LSP peak, get the chi2
			P, LSP = LombScargle(t, resid)
			Pbest = P[np.argsort(LSP)[::-1][0]]

			chi2_new = get_chi2_pf(t,resid,Pbest)
			print "   PERIOD %d = %.5fd   chi2: %.3e --> %.3e [-dchi2 = %.3e]"%( len(PERIODS),  PERIODS[-1], chi2s[-1], chi2_new, chi2s[-1] - chi2_new)
			print "                         std : %.3e --> %.3e [-dstd  = %.3e] = %.3f%% (%.3f%% of total rms)"%(std_old, std_new, std_old - std_new, 100*(std_old - std_new)/std_old, 100*(std_old - std_new)/std_total)

	ws, amps, phs, c = fit_params
	return ws, amps, phs, c, chi2s, stds
	#plot_fit(times,mags,ws,amps,phs,c,show=True)
def SplitFitParam( par ):
	return [ [ par[i*npers + j] for i in range(nharmonics) ] for j in range(npers) ]
def SplitFitParams( ws, amps, phs, c ):
	WS = SplitFitParam(ws)
	AMPS =  SplitFitParam(amps)
	PHS = SplitFitParam(phs)
	C = [ c for i in range(npers) ] 

	return WS, AMPS, PHS, C
def GetLombScargleFeatures(t, x, prefix='', npeaks = NPEAKS_TO_SAVE, add_stetson=True, add_dworetsky=True):
	p, lsp = LombScargle(t, x)
	# Find peaks
	ppers, ppows = find_n_peaks(p, lsp, npeaks)
	inds = np.argsort(ppows)[::-1]
	ppers_sorted = ppers[inds]
	ppows_sorted = ppows[inds]
	features = {}
	# Save the (period, power) pairs for the top few peaks
	for j in range(npeaks):
		P = ppers_sorted[j]
		features['%slsp_peak%d_power'%(prefix, j+1)] = ppows_sorted[j]
		features['%slsp_peak%d_period'%(prefix, j+1)] = ppers_sorted[j]

		if add_stetson: features['%slsp_peak%d_stetson_strlen'%(prefix, j+1)] = Stetson(t, x, P)
		if add_dworetsky: features['%slsp_peak%d_dworetsky_strlen'%(prefix, j+1)] = Dworetsky(t, x, P)
	return features
def MergeDicts(d1, d2):
	d3 = {}
	for key1 in d1:
		if key1 in d2: raise Exception('MergeDicts error, key "%s" is in both dicts'%(str(key1)))
		d3[key1] = d1[key1]
	for key2 in d2: d3[key2] = d2[key2]
	return d3
def GetPeriodicFeatures(t, x):
	features = {}

	# Save Lomb-Scargle results for raw lightcurve
	lsp_raw_features = GetLombScargleFeatures(t, x, prefix='raw_')
	features = MergeDicts(lsp_raw_features, features)

	# Fit (multi-)periodic model to data
	ws, amps, phs, c, chi2s, stds = FitPeriods(t, x)
	
	std = np.std(x)
	# Reorganize arrays Arr = [ [ arr_p1h1, arr_p1h2, ... ], [ arr_p2h1, ... ] ]
	WS, AMPS, PHS, C = SplitFitParams(ws, amps, phs, c)

	# List of periods
	Pers = [ 2*np.pi/w[0] for w in WS ]

	# The chi2 for the phase-folded lightcurve around the best period
	features['chi2_pf'] = chi2s[0]
	features = {}
	for i in range(len(Pers)):

		features['p%d'%(i+1)] = Pers[i]

		# resid -- lightcurve whitened for all P != Pers[i]
		# v-3 ignores this: resid = np.zeros(len(x))
		# v-3 ignores this: resid[:] = x[:]
		# v-3 ignores this: for j in range(len(Pers)):
		# v-3 ignores this: 	if j == i: continue
		# v-3 ignores this: 	resid[:] = Whiten(t, resid, WS[j], AMPS[j], PHS[j], C[j])[1][:]


		# v-3 ignores this: lsp_features = GetLombScargleFeatures(t, resid, prefix='p%d_'%(i+1))
		# v-3 ignores this: features = MergeDicts(lsp_features, features)

		# Total amplitude of this component
		features['p%d_total_amplitude'%(i+1)] = GetAmplitude(WS[i], AMPS[i], PHS[i], C[i])

		# Harmonic breakdown of this component
		for j in range(nharmonics):
			features['p%dh%d_amplitude'%(i+1, j+1)] = AMPS[i][j]
			features['p%dh%d_phase'%(i+1, j+1)] = PHS[i][j]

		features['p%d_frac_of_rms_explained'%(i+1)] = ( std - stds[i] )/std
		# v-3 ignores this: features['p%d_chi2_pf'%(i+1)] = get_chi2_pf(t, resid, Pers[i]) 

	# get red. chi2 for residual, folded around the best LSP period
	# v-3 ignores this: resid = get_resid(t, x, ws, amps, phs, c)
	# v-3 ignores this: pr, powr = LombScargle(t, resid)
	# v-3 ignores this: pbest = pr[np.argsort(powr)[::-1][0]]

	# v-3 ignores this: features['resid_chi2_pf'] = get_chi2_pf(t, resid, pbest)

	# Save Lomb-Scargle results for residual

	# v-3 ignores this: features = MergeDicts(features, GetLombScargleFeatures(t, resid, prefix='resid_'))

	return features
def GetMagnitudeFeatures(lc):
	return {
		'V' 			: lc['mags'][0],
		'R-V'			: lc['mags'][1] - lc['mags'][0],
		'I-V'			: lc['mags'][2] - lc['mags'][0],
		'J-V'			: lc['mags'][3] - lc['mags'][0],
		'H-V'			: lc['mags'][4] - lc['mags'][0],
		'K-V'			: lc['mags'][5] - lc['mags'][0]
	}
def GetCoordinateFeatures(lc):
	return {
		'ra'		: lc['ra'],
		'dec'		: lc['dec']
	}
def FracBeyondNStd(t, x, N):
	sig = np.std(x)
	mu = np.mean(x)
	n = 0
	for X in x: 
		if abs((X - mu)/(N*sig)) > 1.: n+=1
	return float(n)/float(len(x))
def GetVariabilityFeatures(t, x ):
	var_feats =  {
		'std' 							: np.std(x),
		'median_absolute_deviation' 	: np.median(np.abs(x - np.median(x))),
		'skew'							: skew(x),
		'kurtosis' 						: kurtosis(x),
		#'beyond1std' 					: FracBeyondNStd(t, x, 1.0),
		#'beyond2std' 					: FracBeyondNStd(t, x, 2.0),
		#'beyond3std' 					: FracBeyondNStd(t, x, 3.0),
		#'beyond4std' 					: FracBeyondNStd(t, x, 4.0),
		#'beyond5std' 					: FracBeyondNStd(t, x, 5.0),
	}	
	edges, distro = BinnedDistro(x)
	for i in range(len(distro)):
		var_feats['binned_distro_z%.2f_%.2f'%(edges[i], edges[i+1])] = distro[i]
	return var_feats
def get_features(lc, nfreqs = npers, nharmonics=nharmonics, loud=True):
	Features = {}
	if loud: print "detrending..."
	if loud: t0 = time()
	detrend(lc)
	if loud: dt = time() - t0
	if loud: print "  done (%.3f s)"%(dt)

	if loud: print "selecting filters..."
	if loud: t0 = time()
	t, x = get_t_x(lc)
	if t is None: return None
	if len(t) < min_nobs: return None
	if loud: dt = time() - t0
	if loud: print "  done (%.3f s)"%(dt)
	if loud: print "getting periodic features..."
	if loud: t0 = time()
	Features = MergeDicts(Features, GetPeriodicFeatures(t, x))
	if loud: dt = time() - t0
	if loud: print "  done (%.3f s)"%(dt)
	if loud: print "getting magnitude features..."
	if loud: t0 = time()
	Features = MergeDicts(Features, GetMagnitudeFeatures(lc))
	if loud: dt = time() - t0
	if loud: print "  done (%.3f s)"%(dt)
	if loud: print "getting variability features..."
	if loud: t0 = time()
	Features = MergeDicts(Features, GetVariabilityFeatures(t, x))
	if loud: dt = time() - t0
	if loud: print "  done (%.3f s)"%(dt)

	if loud: print "getting coordinate features..."
	if loud: t0 = time()
	Features = MergeDicts(Features, GetCoordinateFeatures(lc))
	if loud: dt = time() - t0
	if loud: print "  done (%.3f s)"%(dt)
	#features['linear_slope'] = fit_out_linear_trend(t, x)['slope'] # This is risky. You should subtract off a linear trend from each (FLT, ...) subset,
																   # THEN detrend
	
	#features['ra'] = lc['ra']
	#features['dec'] = lc['dec']
	
	#if loud: print "     binned distro"
	#distro_bins = get_binned_distro(t, x)
	#for i, d in enumerate(distro_bins):
	#	features['mag_pdf_bin%d'%(i)] = d 
	
	return Features

if __name__ == '__main__':
	def process(ID):
		feat_fname = FeatureVectorFileName(ID)
		if os.path.exists(feat_fname) and not FORCE_REDO:
			return

		print "PID: %s; started"%(str(os.getpid()))
		lc = rhlc.read_hatlc(hat_fname(ID))
		if lc is None:
			return
		t0 = time()
		feat = get_features(lc)
		dt = time() - t0
		print "PID: %s; %s DONE (total time: %.3f s)"%(str(os.getpid()), ID, dt)
		
		pickle.dump(feat, open(feat_fname, 'wb'))
	p = Pool(NPROC)
	ids_fname = "ids_to_use.list"
	if os.path.exists(ids_fname) and not FORCE_IDLIST_REDO:
		ids_to_use = pickle.load(open(ids_fname, 'rb'))
	else:
		ids_to_use = []
		for ID in gcvs_m:
			if is_available(ID): ids_to_use.append(ID)
		pickle.dump(ids_to_use, open(ids_fname, 'wb'))
	print len(ids_to_use)," out of ", len(gcvs_m), " ID's available."
	#print "Fetching %d ids"%(len(ids_to_fetch))
	#fetch_lcs(ids_to_fetch)
	p.map(process, ids_to_use)
	