import pandas as pd 
import numpy as np
from math import *
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import zscore as ZSCORE
from fastperiod import specwindow, lombscargle
from paramiko import SSHConfig, SSHClient
from sklearn.lda import LDA
from lcutils.lcutils_config import *
from settings import *
import sys, os, re, gzip
import fastlombscargle as lsp
from os.path import exists
import dill as pickle
import readhatlc as rhlc
from mpi4py import MPI
import masterslave as msl

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
ROOT = (rank == 0)


HP = pickle.HIGHEST_PROTOCOL
def get_logfile(rank):
	return "%s/logfile.%04d"%(SCRATCH, rank)
save_log_of_each_node = True
terminal_printing = True
nodes_to_print = [ 0, 1 ]
def logprint(m, all_nodes=False):
	if VERBOSE and save_log_of_each_node:
		msg = "node %d: %s"%(comm.rank, m)
		f = open(get_logfile(rank), 'a')
		f.write("%s\n"%(m))
		f.close()
	if VERBOSE and all_nodes and terminal_printing and comm.rank in nodes_to_print: 
		print "node %d: %s "%(comm.rank, m)
	elif VERBOSE and ROOT and terminal_printing: print m


field_info = pickle.load(open(field_info_fname, 'rb'))
all_fields = [ F for F in field_info ]
hatid_field_list = {}
if not os.path.exists(gcvs_info_file):
	print "making gcvs info file..."
	
	hatid_field_list = pickle.load(open(hatid_field_list_fname, 'rb'))

	print "loading gcvs ids..."
	gcvs_ids = pickle.load(open(good_gcvs_hatids_fname, 'rb'))

	print "making info..."
	gcvs_info = {}
	for hatid in gcvs_ids:
		if not hatid in hatid_field_list: 
			gcvs_info[hatid] = None
		else:
			gcvs_info[hatid] = hatid_field_list[hatid]

	print "dumping.."
	pickle.dump(gcvs_info, open(gcvs_info_file, 'wb'))

else:
	gcvs_info = pickle.load(open(gcvs_info_file, 'rb'))


# Load hatids for each field
for field in fields_to_analyze:
	if field == 'gcvs': 
		for hatid in gcvs_info:
			hatid_field_list[hatid] = gcvs_info[hatid]
	else: 
		
		fname = os.path.join(hatids_in_fields_dir, "hatids_in_field_%s.list"%(field))
		'''
		if not os.path.exists(fname):
			field_ids = []
			if not os.path.exists(hatid_field_list_fname):
				get_hatid_field_list(fields_to_analyze)

				print "miscutils: loading hatid_field_list.."
				hatid_field_list = pickle.load(open(hatid_field_list_fname, 'rb'))

			for hatid in hatid_field_list:
				if hatid_field_list[hatid] == field: field_ids.append(hatid)

			pickle.save(field_ids, open(fname, 'wb'))
		'''
		field_ids = pickle.load(open(fname, 'rb'))
		if field_ids is None: 
			hatid_field_list[hatid] = []
			continue
		
		for hatid in field_ids:
			hatid_field_list[hatid] = field

def get_field_of(hatid):
	if not hatid in hatid_field_list: return None
	else: return hatid_field_list[hatid]

get_raw_lc_fname = lambda hatid :  "%s/%s/%s.tfalc"%(LCCACHE, get_field_of(hatid), hatid)
get_lc_fname = lambda hatid :  "%s/%s-full.tfalc.gz"%(LCCACHE, hatid)
#hatid_field_list = {}

def get_iteration_number():
	iteration = 1
	while os.path.exists(get_classifier_fname(iteration)): iteration += 1
	return iteration - 1
def nancleaned(arr):
	#print arr
	return np.array([ a for a in arr if not np.isnan(a) ])
def hat_fname(hid):
	return "%s/%s-hatlc.csv.gz"%(data_dir, hid)
def fit_out_linear_trend(t, x):
	params, pcov = curve_fit(lambda t, y0, m : m*t + y0, t, x)
	x -= params[0] + params[1]*t
	return {
		'slope' : params[1],
		'mag0' : params[0]
	}
def fetch_lcs(hatids, verbose=True):
	i = 0
	lh = len(hatids)
	while i < lh:
		j = i+grpsize
		if j > lh: j = lh
		IDgrp = hatids[i:j]
		filtered_grp = IDgrp
		#filtered_grp = []
		#for ID in IDgrp:
		#	if not os.path.exists(hat_fname(ID)): 
		#		filtered_grp.append(ID)
		has_all = (len(filtered_grp) == 0)
		if has_all: continue
		idlist=""
		for i,ID in enumerate(filtered_grp):
			idlist = "%s%s"%(idlist,ID)
			if i < len(filtered_grp) - 1: idlist = "%s,"%(idlist)
		url = 'https://hatsurveys.org/lightcurves/lc/direct?hatid=%s&apikey=%s'%(idlist, API_KEY['lcdirect'])
		if verbose: print url
		fetch_command = "curl -J -O '%s'"%(url)
		unzip_command = "unzip *zip; rm *zip"
		os.system(fetch_command)
		if len(filtered_grp) > 5: os.system(unzip_command)
		os.system("mv *hatlc*gz %s"%(data_dir))
		i += j
def get_t_x(lc, coltype='TF',selection='locally-brightest', ttype='BJD'):
	t = []
	x = []
	ytypes = [ "%s%d"%(coltype, j+1) for j in range(0,3) ]


	if selection == 'locally-brightest':
		# use the aperture with the brightest flux at each observation
		ndatas = 0
		for i in range(len(lc[ttype])):
			obs = nancleaned([ lc[y][i] for y in ytypes  ])
			if len(obs) == 0: continue
			ndatas += 1
			t.append(float(lc[ttype][i]) - float(lc[ttype][0]))
			x.append(min(obs))

	elif selection == 'globally-brightest':
		# find the aperture with the brightest mean flux
		brightest_ytype = ytype[0]
		brightest_mag = np.mean(nancleaned(lc[brightest_ytype]))
		for ytype in ytypes:
			if ytype == ytype[0]: continue
			mag = np.mean(lc[ytype])
			if mag < brightest_mag:
				brightest_ytype = ytype
				brightest_mag = mag
		ndatas = 0
		for i in range(len(lc[ttype])):
			X = float(lc[brightest_ytype][i])
			if np.isnan(X): continue
			ndatas += 1
			x.append(X)

			t.append(float(lc[ttype][i]) - float(lc[ttype][0]))

	elif selection in [1, 2, 3]:
		# use a specified aperture
		ytype = "%s%d"%(coltype, selection)
		ndatas = 0
		for i in range(len(lc[ttype])):
			X = float(lc[ytype][i])
			if np.isnan(X): continue
			ndatas += 1
			x.append(X)
			t.append(float(lc[ttype][i]) - float(lc[ttype][0]))

	else:
		raise ValueError("I dont understand {0}".format(selection))

	if ndatas < 3: 
		return None, None
	else:
		return np.array(t), np.array(x)
def is_peak(powers, i, imin, imax, per, peak_pers):
    p = powers[i]
    for I in range(imin,imax+1):
        if I == i: continue
        if powers[I] > p: return False
    #for pper in peak_pers:
    #    if abs(per - pper)/pper < 0.01: return False
    return True
def find_n_peaks(periods, powers, n_peaks, dn=7):
    peak_periods, peak_powers = [], []
    inds = np.argsort(powers)[::-1]
    J=0
    n = len(periods)
    while J < len(inds):
        if len(peak_periods) == n_peaks: break
        I = inds[J]
        imin = max([ I - dn, 0 ] )
        imax = min([ I + dn, n - 1])
        if is_peak(powers, I, imin, imax, periods[I], peak_periods): 
            peak_periods.append(periods[I])
            peak_powers.append(powers[I])
            J+= (imax - J)
        else: J+=1
    return np.array(peak_periods), np.array(peak_powers)
def detrend(lc, detrend_vars=[ 'FLT', 'FLD', 'STF', 'NET', 'CAM', 'TEL' ], magtypes = [ 'EP', 'TF' ]):
	for magtype in magtypes:
		for i in range(1, 4):
			magcol = "%s%d"%(magtype, i)
			for v in detrend_vars:
				if v not in lc: continue
				categories = np.unique(lc[v])
				if len(categories) == 1: continue

				nmembers = {}
				mus = {}
				rmses = {}
				most_populous_cat = None
				for cat in categories:
					MAGS = nancleaned(lc[magcol][np.where(lc[v] == cat)])
					nmembers[cat] = len(MAGS)
					if len(MAGS) > 2:
						mus[cat] = np.mean(MAGS)
						rmses[cat] = np.std(MAGS)
					else:
						mus[cat] = None
						rmses[cat] = None
				
				if len([ rmses[cat] for cat in categories if rmses[cat] is not None ]) == 1: continue

				ntot = sum([ nmembers[cat] for cat in categories if mus[cat] is not None ])
				mean_mu = sum([ mus[c]*float(nmembers[c])/float(ntot) for c in categories if mus[c] is not None ])
				mean_rms = sum([ rmses[c]*float(nmembers[c])/float(ntot) for c in categories if rmses[c] is not None ])
				for cat in categories:
					if mus[cat] is None: continue
					lc[magcol][np.where(lc[v] == cat)] -= mus[cat] - mean_mu
					if rmses[cat] > 0:
						lc[magcol][np.where(lc[v] == cat)] = (lc[magcol][np.where(lc[v] == cat)] - mean_mu)*(mean_rms/rmses[cat]) + mean_mu
def false_alarm_prob_to_sig(fap):
	MAX_SIG = 1E15

	if fap < pow(MAX_SIG + 1,-1): return MAX_SIG
	else: return 1./fap - 1.
def string_length(t, x , gamma=2.):
	dxdt = 0.
	N = 0
	for i in range(1,len(t)):
		if t[i] - t[i-1] == 0.: continue
		N+=1
		dxdt += abs((x[i] - x[i-1])/(t[i] - t[i-1]))
	dxdt/=N
	A = 1./dxdt

	return sum([ 
				pow(
					pow(abs(t[i] - t[i-1]),gamma) + pow(A*abs(x[i] - x[i-1]),gamma), 
				1./gamma) for i in range(1,len(t)) ]) #+ pow( 
					#pow(abs(t[0] - t[-1] + 1),gamma) + pow(A*abs(x[0] - x[-1]),gamma), 
				#1./gamma) ]) 
def bootstrap(x, nboots, bootsize, independent_trials=True):
	if independent_trials and nboots*bootsize > len(x): 
		if nboots == 1: 
			return x
		else:
			raise Exception("bootstrap: independent_trials is true, but nboots*bootsize > len(x)")
	inds = np.arange(len(x))
	np.random.shuffle(inds)
	bx = []
	if independent_trials:
		for i in range(nboots):
			bx.append([ x[j] for j in inds[i*bootsize: (i+1)*bootsize] ])
	else:
		for i in range(nboots):
			bx.append([ x[j] for j in inds[:bootsize] ])
			np.random.shuffle(inds)
	return bx
def bootstrapped_lc(t, x, nboots=n_boots, bootsize=boot_size, independent_trials = True):

	bslc = bootstrap([ [ T, X ] for T, X in zip(t, x) ] , nboots, bootsize, independent_trials=independent_trials )

	#print bslc
	bst = [ [ data[0] for data in sample ] for sample in bslc ]
	bsx = [ [ data[1] for data in sample ] for sample in bslc ]
	
	return bst, bsx
def get_dt(t):
	return np.mean([ t[i] - t[i-1] for i in range(1, len(t)) ])
def bootstrap_lsp(t, x, nboots, bootsize, independent_trials = True, data_range = None):
	if data_range is None:
		bst, bsx = bootstrapped_lc(t, x, nboots, bootsize)
		
	else:
		assert( data_range >= boot_size and data_range <= len(t) )

		start_inds = np.arange(len(t) - data_range)
		np.random.shuffle(start_inds)
		start_inds = start_inds[:n_boots]
		
		data_inds = np.arange(len(t))

		bst = [ ]
		bsx = [ ]
		for i,si in enumerate(start_inds):
			inds = data_inds[si:(si + data_range)]
			np.random.shuffle(inds)
			dinds = inds[:n_boots]

			print max(t[dinds]) - min(t[dinds]), data_range*get_dt(t)
			bst.append(t[dinds])
			bsx.append(x[dinds])
	

	return [ lsp.fasper(np.array(T), np.array(X), ofac, hifac, MACC) for T, X in zip(bst, bsx) ]
def get_errs(x, DI=DELTA_I):
	# Assign error values to magnitude measurements by 
	# calculating the standard deviation of 10 nearby measurements
	I=0
	stds = []
	while I < len(x):
		if I + 2*DI >= len(x):
			stds.extend(np.ones(len(x[I:]))*np.std(x[I:]))
			DI = len(x[I:])
		elif I + DI < len(x):
			stds.extend(np.ones(DI)*np.std(x[I:I+DI]))
			
		I += DI
	# If any stds are zero, give them a sensible value.
	for i in range(len(stds)):
		if stds[i] == 0: stds[i] = np.median(stds)
	return np.array(stds)
def avg_subtracted(tpf, xpf, NBINS=NBINS):
	I=0
	tavgs, xavgs = [],[]
	DI = len(xpf)/NBINS
	while I < len(tpf):
		tavgs.append(np.mean(tpf[I:min([I+DI, len(tpf)])]))
		xavgs.append(np.mean(xpf[I:min([I+DI, len(tpf)])]))
		I+=DI
	avgint = interp1d(tavgs, xavgs)
	def avg_func(t):

		if t < min(tavgs):
			i = 1
			while tavgs[i] == tavgs[0]: i+=1
			m = (xavgs[i] - xavgs[0])/(tavgs[i] - tavgs[0])
			b = xavgs[0] - m*tavgs[0]
			return m*t + b
		elif t > max(tavgs):
			i = len(tavgs) - 2
			while tavgs[i] == tavgs[-1]: i-=1
			m = (xavgs[-1] - xavgs[i])/(tavgs[-1] - tavgs[i])
			b = xavgs[-1] - m*tavgs[-1]
			return m*t + b
			
		else: return avgint(t)
	return tpf, np.array([ x - avg_func(t) for t,x in zip(tpf, xpf) ])
def bin_fixed_width(tpf, xpf, nbins=NBINS, mint=None, maxt=None):
	if mint is None: mint = min(tpf)
	if maxt is None: maxt = max(tpf)

	dT = (maxt-mint)/NBINS
	nums = np.zeros(NBINS)
	avgs = np.zeros(NBINS)
	stds = np.zeros(NBINS)
	all_ = [ [ ] for i in range(NBINS) ]

	for T, X in zip(tpf, xpf):
		i = max([ min([ int((T - mint)/dT), NBINS-1 ]), 0 ])
		all_[i].append(X)
	for i in range(NBINS):
		nums[i] = len(all_[i])
		if nums[i] == 0:
			avgs[i] = np.nan
			stds[i] = np.nan
			continue
		avgs[i] = np.mean(all_[i])
		stds[i] = np.std(all_[i])
	
	return avgs, stds, nums
def bin_phase_folded(t,x,stds,nbins=NBINS):
	# If there's a gap in the phase, this binning function might need to be changed
	nperbin = len(t)/nbins

	# Deal with small amounts of data here
	if nperbin == 0:
		print "Error, not enough data to bin"
		if len(t) < 4: return [ 0. ], [ 1. ]
		NB = 2
		nperbin = len(t)/NB
	else: NB = nbins
	
	bmus = []
	bts = []
	bstds = []
	for i in range(NB):
		I = i*nperbin
		J = (i+1)*nperbin
		# If we're at the last bin, include any stragglers
		if i == NB - 1: 
			sel_stds = stds[I:]
			sel_vals = x[I:]
			sel_ts = t[I:]
		else:
			sel_stds = stds[I:J]
			sel_vals = x[I:J]
			sel_ts = t[I:J]
		bts.append(np.mean(sel_ts))
		bstds.append(sqrt(sum([ s**2 for s in sel_stds ])/len(sel_stds)))
		bmus.append(np.mean(sel_vals))
	return np.array(bts), np.array(bmus), np.array(bstds)
def get_binned_distro(t, x, nbins=8):
	sorted_x = np.zeros(len(x))
	sorted_x[:] = x[:]
	
	np.sort(sorted_x)
	nperbin = len(sorted_x)/nbins
	leftover = len(sorted_x) - nperbin*nbins
	
	mu = np.mean(sorted_x)
	std = np.std(sorted_x)
	dev = np.zeros(nbins)

	I = 0
	J= 0 
	for i in range(nbins):
		if i == 0: dJ = nperbin + leftover/2
		elif i == nbins - 1: dJ = len(sorted_x) - J
		else: dJ = nperbin
		J += dJ 
		
		dev[i] = (np.mean(sorted_x[I:J]) - mu)/std
		I += dJ

	return dev
def beyond_nstd(t,x, n):
	std = np.std(x)
	mean = np.mean(x)
	return float(len([ X for X in x if abs(X - mean)/std > n ]))/float(len(x))
def fitting_function(*args):
	# Fits multiple periods!
	t = args[0]
	nharms = args[1]
	nprs = args[2]
	angular_frequencies = [ args[3 + i] for i in range(nprs) ]
	dwdt = args[3 + nprs]
	c = args[4 + nprs]
	o = 5 + nprs
	val = c*np.ones(len(t))
	T = max(t) - min(t)
	for i,w in enumerate(angular_frequencies):
		for h in range(nharms):
			A = o + 2*(i*nharms + h)
			B = o + 2*(i*nharms + h) + 1
			n = h + 1
			W = (w + dwdt*(t - T/2))
			val += args[A]*np.cos(n*W*t) + args[B]*np.sin(n*W*t)
	return val
def phase_fold(t,x,p,dp=DPHASE, dwdt=0.0):
	#if dwdt == 0.: dpdt = 0.
	#else: p = 2pi/w; w = 2pi/p dw = -2pi/p^2 dp

	nparr = np.zeros(len(t), dtype=np.dtype([('tpf', np.float_),( 'x', np.float_)]))
	T0 = min(t)
	TOBS = max(t) - T0
	if dwdt == 0.:
		for i, T in enumerate(t):
			nparr[i]['tpf'] = dp*((T - T0)/(dp*p) - int((T - T0)/(dp*p)))
			nparr[i]['x'] = x[i]
	else:
		W = 2*np.pi/p

		for i, T in enumerate(t):
			P = 2*np.pi/( W + dwdt*(T - TOBS/2) )
			nparr[i]['tpf'] = dp*((T - T0)/(dp*P) - int((T - T0)/(dp*P)))
			nparr[i]['x'] = x[i]
	nparr.sort(order='tpf')
	return nparr['tpf'], nparr['x']
def symmetry_measure(t, x, p, nbins=NBINS):
	tpf, xpf = phase_fold(t, x, p, dp=2.0)

	xpf_avgs_1, xpf_err_1, xpf_nums_1 = bin_fixed_width(tpf, xpf, nbins=nbins, mint=0.0, maxt=1.0)
	xpf_avgs_2, xpf_err_2, xpf_nums_2 = bin_fixed_width(tpf, xpf, nbins=nbins, mint=1.0, maxt=2.0)

	vals = []
	for x1,x2,e1,e2,n1,n2 in zip(xpf_avgs_1, xpf_err_1, xpf_nums_1, xpf_avgs_2, xpf_err_2, xpf_nums_2):
		if np.isnan(x1) or np.isnan(x2): continue
		vals.append( (x1 - x2)**2/(e1**2/n1 + e2**2/n2) )
	if len(vals) < 3: return None
	else: return -sum(vals)/(len(vals) - 2)
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
def SplitFitParam( par ):
	return [ [ par[i*npers + j] for i in range(nharmonics) ] for j in range(npers) ]
def SplitFitParams( ws, amps, phs, c ):
	WS = SplitFitParam(ws)
	AMPS =  SplitFitParam(amps)
	PHS = SplitFitParam(phs)
	C = [ c for i in range(npers) ] 

	return WS, AMPS, PHS, C
def translate_features_to_popt(features, use_dwdt=False):
	popt = []
	popt.append(features['constant_offset'])
	for p in range(npers):
		for h in range(nharmonics):
			phi = features['p%dh%d_phase'%(p+1,h+1)]
			amp = features['p%dh%d_amplitude'%(p+1, h+1)]
			popt.append(amp*cos(phi))
			popt.append(amp*sin(phi))
	return popt
def translate_popt_to_fit_params(popt, use_dwdt=False):
	if use_dwdt: ADD = 2
	else: ADD = 1
	As = [ popt[ ADD + 2*(freq_no*nharmonics + harm_no)  ]
				for harm_no in range(nharmonics) 
			for freq_no in range(npers) ]
	Bs = [ popt[ ADD + 2*(freq_no*nharmonics + harm_no) + 1  ]
				for harm_no in range(nharmonics) 
			for freq_no in range(npers) ]				  
	


	def phase_shift(A, B): 
		if A == 0 and B > 0: return np.pi/2
		elif A == 0 and B < 0: return -np.pi/2
		elif A < 0:
			return np.arctan(B/A) + np.pi
		else: return np.arctan(B/A)

	# popt[0] is the constant offset.
	# then it's (A,B) pairs for each harmonic, for each frequency: 
	#		popt = [ C, A_f1h1, B_f1h1, A_f1h2, B_f1h2, ..., A_f1hNH, B_f1NH, A_f2h1, B_f2h1, ... ]

	amplitudes = 	[ 	sqrt(A**2 + B**2) 
						for A,B in zip(As, Bs) ]

	phases =     	[ 	phase_shift(float(A),float(B)) 
						for A,B in zip(As, Bs) ]

	if use_dwdt:
		return amplitudes, phases, popt[0], popt[1]
	else:
		return amplitudes, phases, popt[0]
def translate_features_to_fit_params( features, use_dwdt = False):
	return translate_popt_to_fit_params(translate_features_to_popt(features, use_dwdt=use_dwdt) ,use_dwdt=use_dwdt)
def translate_fit_params_to_features(amps,phs,c):
	features = {}
	features['constant_offset'] = c
	WS, AMPS, PHS, C = SplitFitParams(np.ones(len(amps)), amps, phs, c)
	for i in range(npers):
		features['p%d_total_amplitude'%(i+1)] = GetAmplitude(WS[i], AMPS[i], PHS[i], C[i])

		# Harmonic breakdown of this component
		for j in range(nharmonics):
			features['p%dh%d_amplitude'%(i+1, j+1)] = AMPS[i][j]
			features['p%dh%d_phase'%(i+1, j+1)] = PHS[i][j]
	return features
def translate_popt_to_features(popt, use_dwdt=False):
	fpars = translate_popt_to_fit_params(popt, use_dwdt=use_dwdt)
	#print fpars
	return translate_fit_params_to_features(*fpars)
def fit_periods(t,x,ps, use_dwdt=False, return_errors=False):
	ws = [ 2*np.pi/float(p) for p in ps ]
	frqs = [ ws[freq_no] * (harm_no + 1)
				for harm_no in range(nharmonics) 
			for freq_no in range(npers) ]

	if use_dwdt:
		ADD = 2
		ff = lambda T, *args : fitting_function(T, nharmonics, len(ps), *(ws + args) )
	else:
		ADD = 1
		def ff( T, *args ):
			new_args = [ 0. ]
			new_args.extend(args)
			new_args = tuple(new_args)
			return fitting_function(T, nharmonics, len(ps), *(tuple(ws) + new_args) )
		
	# p0 is the initial guess for each parameter.
	p0 = np.ones(2*nharmonics*len(ws) + ADD)
	p0[0] = np.mean(x) # guess a constant offset of <magnitude> (this isn't crucial)
	if use_dwdt: p0[1] = 0.
	# fit!
	popt, pcov = curve_fit(ff, t, x, p0=p0, absolute_sigma=True)
	fit_pars = translate_popt_to_fit_params(popt, use_dwdt=use_dwdt)
	if return_errors:

		return tuple([frqs]) + fit_pars + tuple([pcov])
	else:
		return tuple([frqs]) + fit_pars
def bs_fit_periods(t, x, ps, nsub=boot_size):
	bst, bsx = bootstrapped_lc(t, x, nboots=1, bootsize=nsub)

	bst = bst[0]
	bsx = bsx[0]

	return fit_periods(bst, bsx, ps)
def correct_for_window(LSP0, WLSP, NSIG=5):
	# Find outlier-free std and mean of the window
	sig = np.std(WLSP)
	mu = np.mean(WLSP)
	LSPF = np.zeros(len(LSP0))
	WLSP_no = [ L for L in WLSP ]
	while True:
		N = len(WLSP_no)
		WLSP_no = [ L for L in WLSP_no if L - mu < NSIG*sig ]
		if N - len(WLSP_no) == 0: break
		
		mu = np.mean(WLSP_no)
		sig = np.std(WLSP_no)

	# Set all LSP power to zero in places where the window
	# function has a high value
	for i in range(len(LSP0)):
		if WLSP[i] - mu > NSIG*sig:
			LSPF[i] = mu
		else:
			LSPF[i] = LSP0[i]
	return LSPF
def get_resid(t, x, ws, amps, phs, c, dwdt=0.0):
	tobs = max(t) - min(t)
	#DW=  -dwdt*(max(t) - min(t))/2.0
	return x - np.array([ c + sum([  A*cos((W+dwdt*(T - tobs/2.))*T - phi) for A,W,phi in zip(amps, ws, phs) ]) for T in t ])
def find_outliers(t,x,sigma_threshold):
	assert(len(t) == len(x))
	zscore = np.zeros(len(t))
	zscore[:] = avg_subtracted(t,x)[1][:]
	sigma, mu = np.std(zscore), np.mean(zscore)
	zscore[:] = (zscore[:] - mu)/sigma
	return [ i for i,z in enumerate(zscore) if abs(z) > sigma_threshold ]
def remove_outliers(t,x,o, return_outliers = False):
	tp, xp, tout, xout = [], [], [], []
	I=0
	outliers = np.sort(o)
	for i in range(len(t)):
		if I == len(outliers): 
			for j in range(i,len(t)):
				tp.append(t[j])
				xp.append(x[j])
			break
		if i == outliers[I]: 
			tout.append(t[i])
			xout.append(x[i])
			I+=1
			continue
		tp.append(t[i])
		xp.append(x[i])
	if return_outliers: return tp, xp, tout, xout
	else: return tp, xp
def prune_outliers(t, x, sigma_threshold=5):
	assert(len(t) == len(x))
	outliers = []
	tp, xp = [],[]
	tp.extend(t)
	xp.extend(x)
	tpru, xpru = [],[] #
	tout, xout = [],[] #
	tmod, xmod = [],[] #
	while True:
		#print "Iter %d"%(len(tpru)+1)
		new_outliers = find_outliers(tp, xp, sigma_threshold)
		tmod.append(tp) #
		xmod.append(np.subtract(xp,avg_subtracted(tp,xp)[1])) #
		if len(new_outliers) == 0: break
		tp, xp, to, xo = remove_outliers(tp,xp,new_outliers, return_outliers=True)
		tpru.append(tp) #
		xpru.append(xp) #
		tout.append(to) #
		xout.append(xo) #
		outliers.extend(new_outliers)


	outliers = np.sort(outliers)
	'''
	NR = 3
	FW = 4
	nrows = len(tpru)/NR + len(tpru)%NR
	ncols = 0
	while nrows*ncols < len(tpru): ncols+=1
	f, axes = plt.subplots(nrows,ncols,figsize=(FW*nrows, FW*ncols))

	for i in range(len(tpru)):
		if ncols > 1: ax = axes[i/ncols][i%ncols]
		else: ax = axes[i]
		ax.set_title("Iter %d"%(i+1))
		ax.scatter(tpru[i],xpru[i],color='k',marker=',',alpha=0.1)
		ax.scatter(tout[i],xout[i],color='r',marker=',',alpha=1.0)
		ax.plot(tmod[i],xmod[i], color='r')
		ax.set_xlim(0,DPHASE)
		ax.invert_yaxis()

	plt.show()

	'''
	return np.array(tp), np.array(xp)
def red_chi2(t,y,yfit,sigma,ndof):
	assert(all(sigma>0))
	assert(len(t) == len(y))
	nu = len(t) - ndof - 1
	return sum(((y - yfit)/sigma)**2)/nu
def get_chi2_pf(t,x,p, DI=DELTA_I,nbins=NBINS):

	tpf, xpf = phase_fold(t, x, p)
	
	# Error estimates for phase-folded lightcurve
	stds0 = get_errs(xpf,DI=DI)

	# Now bin!
	tbins, mus, stds = bin_phase_folded(tpf,xpf,stds0,nbins=nbins)

	# Make sure all of the stds > 0
	std = []
	for s in stds:
		if s < eps: std.append( eps )
		else: std.append(s)
	std = np.array(std)

	# Fit a horizontal line to the results
	popt, pcov = curve_fit(lambda T, A : A*np.ones(len(T)),tbins , mus,  p0=[ np.mean(mus) ], sigma=std)
	A = popt[0]

	nu = len(tbins) - 2

	return red_chi2(tbins, mus, A*np.ones(len(tbins)), std, 1)
	#return sum(np.power(np.divide((mus - A),std),2))/float(nu)
def LN_norm(x, n):
	mu = np.mean(x)
	return pow(np.mean( abs(x - mu)**n ), 1./n)
def find_best_of_n_periods(t,x,pers,other_periods=[],fraction_of_smallest_stds_to_use=fraction_of_smallest_stds_to_use):
	sigs = []
	scores = []
	#print pers
	if not use_bootstrap:
		BOOTSIZE = len(t)
		NBOOTS   = 1
	else:
		BOOTSIZE = boot_size
		NBOOTS = n_boots
	for p in pers:
		tot_pers = [ p ]
		tot_pers.extend(other_periods)
		bst, bsx = bootstrapped_lc(t, x, nboots=NBOOTS, bootsize=min([ len(x)/NBOOTS, BOOTSIZE ]))
		sses = []
		for bt, bx in zip(bst, bsx):
			ws, amps, phs, c = fit_periods(bt, bx, tot_pers, use_dwdt=False)
			
			resid = get_resid(bt, bx, ws, amps, phs, c)
			
			sses.append(np.std(resid))
		scores.append(np.mean(sses))
		
	pers_sorted = np.array(pers)[np.argsort(scores)]
	bp = pers_sorted[0]
	
	return [ bp ], []
def ZipAndMerge(x1, x2):
	X = [ ]
	assert(len(x1) == len(x2))

	for X1,X2 in zip(x1, x2):

		if not hasattr(X1, '__iter__'):  X1 = [ X1 ]
		if not hasattr(X2, '__iter__'):  X2 = [ X2 ]
		x = []
		x.extend(X1)
		x.extend(X2)
		X.append(x)
	return np.array(X)
def make_obs(Feats, keylist=None):

	Obs = []
	for F in Feats:
		obs = []
		if keylist is None:
			for k in Feats[F]:
				obs.append(Feats[F][k])
		else:
			for k in keylist:
				obs.append(Feats[F][k])
		Obs.append(obs)
	return np.array(Obs)


def make_keylist_dict(klist_data):
	keylist = {}
	for kl in klist_data:
		keylist[kl['TSTF']] = { }
		for k in keylist_dt.names:
			if k == 'TSTF': continue
			keylist[kl['TSTF']][k] = kl[k]
	return keylist
def get_keylist(field, ssh_host_name=ssh_host_name):
	
	local_fname = get_local_keylist_fname(field)

	os.system("scp %s:%s/keylist.txt %s"%(ssh_host_name, get_keylist_dir(field), local_fname))
	klist_data = np.loadtxt(local_fname, dtype=keylist_dt)
	return make_keylist_dict(klist_data)
def load_2mass(field, twomass_remote_dir='~', twomass_fname='colors_field%s.dat', ssh_host_name=ssh_host_name):
	twomass_dict = {}

	if not os.path.exists(twomass_local_fname):
		os.system("scp %s:%s %s"%(ssh_host_name, get_remote_2mass_fname(field), get_local_2mass_fname(field)))
	twomass_data = np.loadtxt(get_local_2mass_fname(field), dtype=twomass_dt)
	#print twomass_data['hatid']
	for tm in twomass_data:
		#print tm, tm['hatid']
		#sys.exit()
		twomass_dict[tm['hatid']] = {}
		for c in twomass_dt.names:
			#print c
			if c == 'hatid': continue
			twomass_dict[tm['hatid']][c] = tm[c]

	return twomass_dict
def fix_times(lc, tcols = [ 'BJD' ]):
	# just checks that the LC entries are in chronological order;
	for tcol in tcols: assert(all([ lc[tcol][i+1] > lc[tcol][i] for i in range(len(lc[tcol] - 1)) ]) ) 
	return lc
def add_2mass(lc, tmdat):
	lc['ra'] = tmdat['ra']
	lc['dec'] = tmdat['dec']
	lc['mags'] = [ tmdat[x] for x in ('Vmag', 'Rmag', 'Imag', 'jmag', 'hmag', 'kmag') ]
	lc['ugriz'] = [ tmdat[x] for x in ( 'umag', 'gmag', 'rmag', 'imag', 'zmag' )]
	lc['ndet'] = len(lc['TF1'])
	

	lc['twomassid'] = 0
	lc['hatstations'] = np.unique(lc['STF'])
	lc['filters'] = [ flt for flt in np.unique(lc['FLT']) ]
	return lc
def add_keylist_data(lc, kl):

	cols = [ 'FLT', 'EXP', 'BJD' ]

	# HATnorth sources have a slightly different format than HATsouth sources;
	#  north TSTFC = <TSTF>_<CCD>
	#  south TSTF  = <TSTF>

	if 'HAT' in lc: 
		is_hatnorth_source = True
	else:
		is_hatnorth_source = False
		
	# initialize
	for col in cols:
		lc[col] = []

	missing_exposures = []
	# For each measurement
	for i in range(len(lc['TF1'])):

		# Obtain CCD/exposure codes
		if is_hatnorth_source:
			EXP, CCD = lc['TSTFC'][i].split('_')
		else:
			EXP = lc['TSTF'][i]
		not_available = (EXP not in kl)
		if not_available: missing_exposures.append(EXP)
			#logprint(" ERROR: exposure %s is not in the keylist!!"%(EXP), all_nodes=True)
			
		# Add info from keylist!
		for col in cols:
			if not_available:
				lc[col].append('?')
			else:
				lc[col].append(kl[EXP][col])
	if len(missing_exposures) > 0: logprint("   add_keylist_data: %d out of %d exposures are missing from the keylist!!"%(len(missing_exposures), len(lc['TF1'])), all_nodes=True)
	# Convert to numpy arrays
	for col in cols:
		lc[col] = np.array(lc[col])

	return lc

def get_colnames(local_fname):
	logprint(" in get_colnames: local_fname = %s"%(local_fname), all_nodes=True)	
	with open(local_fname, 'r') as f:
		for line in f:
			if '#' in line: continue
			val1 = line.split()[0]
			if 'HAT' in val1: return HATLC_COL_DEFS['hs']['tfalc']
			else: return HATLC_COL_DEFS['hn']['tfalc']
	return None

def load_tfalc(local_fname):
	lc = {}
	colnames = get_colnames(local_fname)
	if colnames is None: return None
	for c in colnames: lc[c] = []

	with open(local_fname, 'r') as f:
		for line in f:
			data = line.split("#")
			vals = data[0].split()
			if len(vals) < len(colnames): continue

			for i,c in enumerate(colnames):
				try:
					lc[c].append(TEXTLC_OUTPUT_COLUMNS[c][3](vals[i]))
				except ValueError, e:
					raise Exception("{0} -- {1} not convertable... {2}".format(c, vals[i], local_fname))

	for c in colnames: lc[c] = np.array(lc[c])
	lc['frame'] = []
	lc['STF'] = []
	if 'HAT' in colnames:
		lc['CCD'] = []
		for tstf in lc['TSTFC']:
			station, frameccd = tstf.split('-')
			frame, ccd = frameccd.split('_')
			lc['STF'].append(station)
			lc['CCD'].append(ccd)
			lc['frame'].append(frame)
		lc['CCD'] = np.array([ int(CCD) for CCD in lc['CCD']])
	else:
		for tstf in lc['TSTF']:
			station, frame = tstf.split('-')
			lc['STF'].append(station)
			lc['frame'].append(frame)

	lc['frame'] = np.array([ int(f) for f in lc['frame']])
	lc['STF'] = np.array([ int(stf) for stf in lc['STF']])

	lc['cols'] = colnames
	return lc
def prune_out_bad_inds(lc):
	bad_inds = []
	for i,t in enumerate(lc['BJD']):
		try:
			T = float(t)
		except:
			bad_inds.append(i)

	cols_to_prune = [ c for c in lc if hasattr(lc[c], '__iter__') and len(lc[c]) == len(lc['TF1']) ]
	for c in cols_to_prune:
		lc[c] = np.array([ lc[c][i] for i in range(len(lc['TF1'])) if not i in bad_inds ])

	return lc
def load_full_tfalc(local_fname, keylist_dat, twomass_dat):

	lc = load_tfalc(local_fname)
	if lc is None: return None
	lc = add_keylist_data(lc, keylist_dat)
	if lc is None: return None
	lc = add_2mass(lc, twomass_dat)
	return lc

def get_remote_tfalc_fname(hatid, field = None):
	if field is None:
		field = get_field_of(hatid)
		if field is None: return None
	return "%s/%s.tfalc"%(field_info[field], hatid)

def download_lc(hatid, sftp=None):
	close_it = False
	if sftp is None:
		ssh, sftp = open_ssh_connection()
		close_it = True


	fname_r = get_remote_tfalc_fname(hatid)
	if fname_r is None: 
		logprint("                                   %s  !  fname is none!"%(hatid), all_nodes=True)
		return False

	fname_l = get_raw_lc_fname(hatid)

	if not rexists(sftp, fname_r):
		logprint("                                   %s  !  rexists returns false"%(hatid), all_nodes=True)
		return False

	sftp.get(fname_r, fname_l)
	if close_it: close_ssh_connection(ssh, sftp)

	return True
def save_and_return(lc, hatid, save_full_lc):
	# Save lightcurve
	if save_full_lc: 
		f = gzip.open(get_lc_fname(hatid), 'wb')
		logprint("   ==> saving %s lc to %s"%(hatid, get_lc_fname(hatid)), all_nodes=True)
		pickle.dump(lc, f)
		f.close()
	return lc

def load_full_tfalc_from_scratch(hatid, field=None, keylist_dat=None, twomass_dat=None, ssh=None, sftp=None, save_full_lc=True, min_observations=5, delete_raw_lc=True, force_redo = force_redo):
	logprint("  load_full_tfalc_from_scratch ** %s"%(hatid), all_nodes=True)

	# Load full lightcurve if one is available
	if os.path.exists(get_lc_fname(hatid)): 
		logprint("    (loading) %s"%(hatid), all_nodes=True)
		try:
			f = gzip.open(get_lc_fname(hatid), 'rb')
			full_tfalc = pickle.load(f)
			f.close()
			if not force_redo or ( force_redo and not full_tfalc is None ): return full_tfalc
		except:
			logprint("               !                  %s  ! can't load %s"%(hatid, get_raw_lc_fname(hatid)), all_nodes=True)
	
	# Find field of HATID
	if field is None:
		logprint("                                  %s  >  getting field..."%(hatid), all_nodes=True)

		field = get_field_of(hatid)

	# Make sure .tfalc lightcurve exists on the local system...
	if not os.path.exists(get_raw_lc_fname(hatid)): 
		logprint("                    !             %s  ! No tfalc lightcurve on the system."%(hatid), all_nodes=True)
		return save_and_return(None, hatid, save_full_lc)
	
	# Obtain twomass/color data for hatid
	if twomass_dat is None:

		# Try loading twomass_info_for_field if it isn't loaded already...
		if twomass_info_for_field is None: add_twomass_info_field(field)
		if twomass_info_for_field is None: raise Exception("line 1033 in miscutils: tried to load twomass_info_for_field, but it's STILL None.")
		if not field in twomass_info_for_field: add_twomass_info_field(field)
		
		if is_gcvs(hatid) and 'gcvs' in fields_to_analyze:
			if not hatid in twomass_info_for_field['gcvs']:
				logprint("                    !             %s  !  No twomass information available for hatid in `gcvs` field."%(hatid), all_nodes=True)
				return save_and_return(None, hatid, save_full_lc)
			twomass_dat = twomass_info_for_field['gcvs'][hatid]

		elif not field in twomass_info_for_field:
			
			logprint("                    !             %s  !  No twomass information loaded for field %s.\n\t\t\t You'll need to run make_2mass_info_file.py again for this field"%(hatid, field), all_nodes=True)
			return save_and_return(None, hatid, save_full_lc)
			#return
		elif not hatid in twomass_info_for_field[field]:
			logprint("                    !             %s  !  No twomass information loaded.\n\t\t\t You'll need to run make_2mass_info_file.py again for this hatid"%(hatid), all_nodes=True)
			return save_and_return(None, hatid,  save_full_lc)
		elif twomass_info_for_field[field][hatid] is None:
			logprint("                    !             %s  !  No twomass information is available!", all_nodes=True)
			return save_and_return(None, hatid,  save_full_lc)
		else:
			twomass_dat = twomass_info_for_field[field][hatid]
	# Load keylist for field
	if keylist_dat is None:
		logprint("                                  %s  >  getting keylist data..."%(hatid), all_nodes=True)
		keylist_dat = load_keylist(field)
		if keylist_dat is None:
			logprint("                    !             %s  !  KEYLIST is no good :("%(hatid), all_nodes=True)
			
			return save_and_return(None, hatid,  save_full_lc)

	logprint("                                  %s  >  getting twomass data..."%(hatid), all_nodes=True)

	

	# Load tfalc lightcurve; add twomass/keylist data to it.
	logprint("                                  %s  >  loading lc..."%(hatid), all_nodes=True)
	lc = load_full_tfalc(get_raw_lc_fname(hatid), keylist_dat, twomass_dat)
	if lc is None: 
		logprint("                    !             %s  !  lc was none!"%(hatid), all_nodes=True)
		return save_and_return(None, hatid, save_full_lc)

	logprint("                                  %s  >  pruning out bad inds..."%(hatid), all_nodes=True)
	lc['hatid'] = hatid
	# Prune out missing keylist datapoints:
	lc = prune_out_bad_inds(lc)

	# Ignore this lightcurve if there are too few points
	if len(lc['BJD']) < min_observations: 
		logprint("                    !             %s  !  too few observations!!"%(hatid), all_nodes=True)
		return save_and_return(None, hatid, save_full_lc)
	
	logprint(logprint("                                  %s  >  %d observations are OK!"%(hatid, len(lc['BJD'])), all_nodes=True))
	# To save space, once the pickled lightcurve is generated, the original tfalc is discarded.
	if delete_raw_lc and RUNNING_ON_DELLA:
		logprint("                                  %s  >  DELETING tfalc LC."%(hatid), all_nodes=True)
		os.remove(get_raw_lc_fname(hatid))

	return save_and_return(lc, hatid, save_full_lc)

def rexists(sftp, path):
    """os.path.exists for paramiko's SCP object
    """
    try:
        sftp.stat(path)
    except IOError, e:
        if e[0] == 2:
            return False
        raise
    else:
        return True
def get_field_ids(field, sftp):
	if not field in field_info: return None
	field_dir = field_info[field]
	if not rexists(sftp, field_dir): return None
	ids = []
	for fname in sftp.listdir(field_dir):
		if not 'HAT' in fname: continue
		if not '.tfalc' in fname: continue
		if len(fname) != 21: continue
		ids.append(fname[:15])
	return ids

def is_gcvs(hatid):
	return ( hatid in gcvs_info )
#def get_field_ids(field):
#	return [ hatid for hatid in hatid_field_list if hatid_field_list[hatid] == field ]
def get_available_field_ids(field):

	dirname = "%s/%s"%(LCCACHE, field)
	
	files = os.listdir(dirname)

	field_ids = []
	for f in files:
		# Skip non-tfalc lightcurves
		if not 'tfalc' in f: continue
		# if more than 1 period in filename, skip.
		if len(f.split('.')) > 2: continue
		
		field_ids.append(f.replace('.tfalc', ''))

	return field_ids
def print_list(l):
	pstr = ""
	for i,item in enumerate(l):
		if i == len(l) - 1: pstr = "%s%s"%(pstr, item)
		else: pstr = "%s%s,"%(pstr, item)

	return pstr

def get_remote_tfalc_filename(hatid):
	field = get_field_of(hatid)
	if field is None: return None
	return "%s/%s.tfalc"%(field_info[field], hatid)
def safe_open_keylist(fname):
	try:
		keylist_data = np.loadtxt(fname, dtype = keylist_dt)
		return keylist_data
	except ValueError, e:
		logprint(" WARNING: problem using loadtxt on keylist %s"%(fname), all_nodes=True)
		try:
			kl = open(fname, 'r')
		except:
			logprint(" %s problem summary: open() in 'r' mode does not work"%(fname), all_nodes=True)
			return None


		keys = []
		for key, dt in keylist_dt_arr:
			if 'unknown' in key: continue
			keys.append((key, dt))

		keylist_dict = {key : [] for key, dt in keys}

		lno = 0
		bad_inds = []
		for line in kl:
			if 'NULL' in line: 
				bad_inds.append(lno)
				lno+=1
				continue

			entries = line.split()
			if len(entries) == len(keys):
				for i in range(len(keys)):

					if dt == float: entry = float(entries[i])
					else: entry = entries[i]

					keylist_dict[keys[i][0]].append(entry)
			else:
				bad_inds.append(lno)
				lno+=1
				continue
			lno+=1

		kl.close()
		keylist_dat = np.empty(lno - len(bad_inds), dtype=keylist_dt)

		for i in range(lno - len(bad_inds)):
			for key in keylist_dict:
				keylist_dat[i][key] = keylist_dict[key][i]
		logprint(" %s: %d/%d lines are problematic"%(fname, len(bad_inds), lno), all_nodes=True)


		return keylist_dat
def load_keylist(field):
	klfname_l = get_local_keylist_dict_fname(field)

	if not os.path.exists(klfname_l):
		logprint("Warning: keylist data for field %s is not available!"%(field))
		return None

	keylist = pickle.load(gzip.open(klfname_l, 'rb'))

	return keylist
def make_hatid_list(field, sftp=None, ssh=None):
	# Get list of hatids that are in this field
	close_it = False
	if sftp is None: 
		close_it = True
		ssh, sftp = open_ssh_connection()
	
	field_ids = get_field_ids(field, sftp)
	fname = get_hatids_in_field_fname(field)
	pickle.dump(field_ids, open(fname, 'wb'))

	if close_it: close_ssh_connection(ssh, sftp)

	return field_ids
def make_hatid_lists(fields = None, sftp=None, ssh=None):
	# SSH
	close_it = False
	if ( not ROOT or size == 1 ) and sftp is None:
		ssh, sftp = open_ssh_connection()
		close_it = True

	# Get all missing fields if no list is specified
	logprint("Making hatid lists")
	if fields is None: 
		logprint("Getting missing fields.")
		fields = []
		for field in all_fields:
			fname = get_hatids_in_field_fname(field)
			if not overwrite and os.path.exists(fname): continue
			fields.append(field)
		if len(fields) == 0: 
			logprint("No fields are missing!")
			return None
		logprint("done!")

	# Get list of field ids.
	if size == 1:
		# Serial
		hatid_lists = {}
		for field in fields:
			hatid_lists[field] = make_hatid_list(field, sftp=sftp, ssh=ssh)
		return hatid_lists

	elif ROOT:
		# Master/slave
		all_field_ids = msl.master(fields)

		# save information.
		bad_fields = []
		for field, ids in all_field_ids:
			if ids is None: 
				bad_fields.append(field)
				
			fname = get_hatids_in_field_fname(field)
			pickle.dump(ids, open(fname, 'wb'))

		logprint("%d out of %d fields (%d total fields) have invalid paths!"%(len(bad_fields), len(fields), len(all_fields)))
		logprint(print_list(bad_fields))

	else:
		
		msl.slave(lambda field : ( field, get_field_ids(field, sftp)))

	if close_it: close_ssh_connection(ssh, sftp)
	
get_2mass_binary = lambda field : "/home/jhoffman/2massread"
get_2mass_data_dir = lambda field : "/H/CAT/2MASS/2MASS_JH_AP/data"
get_2mass_shell_script = lambda field : "/home/jhoffman/color-file.sh"
def open_ssh_connection():
	logprint("  open_ssh_connection: in function.", all_nodes=True)
	config = SSHConfig()
	with open(os.path.expanduser('~/.ssh/config')) as config_file:
	    config.parse(config_file)
	d = config.lookup(ssh_host_name)
	logprint("  open_ssh_connection: d['hostname'] = %s, d.get('user') = %s"%(d['hostname'], d.get('user')), all_nodes=True)
	ssh = SSHClient()
	ssh.load_system_host_keys() #NOTE: no AutoAddPolicy() 
	logprint("  open_ssh_connection: loaded system host keys!", all_nodes=True)
	ssh.connect(d['hostname'], username=d.get('user'))
	logprint("  open_ssh_connection: connected.", all_nodes=True)
	sftp = ssh.open_sftp()
	logprint("  open_ssh_connection: opened sftp.", all_nodes=True)

	return ssh, sftp
def close_ssh_connection(ssh, sftp):
	if sftp: sftp.close()
	if ssh: ssh.close()
def conv(r, t):
	if isinstance(t, str): return r
	else: return t(r)

def get_2mass_data_for_hatid_over_ssh(hatid, sftp):
	command = "/home/jhoffman/2massread --cat %s -g %s"%(get_2mass_data_dir(None), hatid)
	stdin, stdout, stderr = sftp.exec_command(command)
	ncols = len(twomass_dt_arr)

	line = stdout.read()
	try:
		cols = line.split('\n')[1].split()
	except:
		print "Can't get line..."
		return None
	
	#print len(cols)
	# Is the number of columns ok?
	if len(cols) != ncols: 
		print "Len cols is %d, but ncols is %d"%(len(cols), ncols)
		return None
		

	# Is the value of each column ok?	
	col_is_bad = False
	for j,col in enumerate(cols):
		name, dt = twomass_dt_arr[j]
		try:
			conv(col, dt)
		except:
			print "col %s is 'bad'"%(name), col
			col_is_bad = True
			break
	if col_is_bad:
		return None

	twomass_dict = {}
	for j,col in enumerate(cols):
		name, dt = twomass_dt_arr[j]
		twomass_dict[name] = conv(cols[j], dt) 

	return twomass_dict

def get_2mass_data_for_hatid(hatid):
	field = get_field_of(hatid)
	if field is None:
		logprint("Warning (get_2mass_data): %s does not have a field."%(field))
		return None
	if twomass_info_for_field is None: add_twomass_info_field(field)
	elif not field in twomass_info_for_field: add_twomass_info_field(field)

	if twomass_info_for_field is None:
		raise Exception("in get_2mass_data_for_hatid in miscutils. Twomass_info_for_field is None even after we tried to load it!!")
	return twomass_info_for_field[field][hatid]
	
def make_2mass_data_for_field(field, client, sftp):

	#if os.path.exists(get_local_2mass_fname(field)):
	#	return np.loadtxt(get_local_2mass_fname(field))

	fname = get_hatids_in_field_fname(field)
	if not os.path.exists(fname):
		logprint("List of HATID's for field %s does not yet exist. Making one!"%(field), all_nodes=True)
		res = make_hatid_lists([ field ])
		if res is None: 
			logprint("Oh wait, we don't know where that field is on the system yet!", all_nodes=True)
			return None
	hatids = pickle.load(open(fname, 'rb'))
	
	results = np.empty(len(hatids), dtype=twomass_dt)
	f = open(get_local_2mass_fname(field), 'w')
	for i, hatid in enumerate(hatids):
		print i, len(hatids), hatid
		res = get_2mass_data_for_hatid(hatid, client, field=field)
		for c in res:
			results[i][c] = res[c]
			f.write("{0} ".format(res[c]))
		f.write("\n")
		
	f.close()

	return np.array(results)
def is_candidate(scores, min_score, min_frac_above_min_score):
	nabove = sum([ 1 for s in scores if s > min_score ])
	if float(nabove)/float(len(scores)) > min_frac_above_min_score: return True
	return False
def load_bad_ids(iteration):
	BAD_IDS = []
	for ID in bad_ids: BAD_IDS.append(ID)
	if iteration == 0: return BAD_IDS
	for i in range(0,iteration):
		bad_ids_ = pickle.load(open(get_bad_ids_fname(i), 'rb'))
		for ID in bad_ids_: BAD_IDS.append(ID)
	return BAD_IDS
def get_hatid_field_list(fields):
	logprint("Getting hatid field list!", all_nodes=True)
	ssh, sftp = open_ssh_connection()
	for field in fields:
		fname = get_hatids_in_field_fname(field)
		if not os.path.exists(fname):
			logprint("List of HATID's for field %s does not yet exist. Making one!"%(field), all_nodes=True)
			res = make_hatid_list(field, sftp=sftp)
			if res is None: 
				logprint("Oh wait, we don't know where that field is on the system yet!", all_nodes=True)
				continue
		hatids = pickle.load(open(fname, 'rb'))
		if hatids is None: continue

		for hatid in hatids:
			hatid_field_list[hatid] = field
	pickle.dump(hatid_field_list, open(hatid_field_list_fname, 'wb'))
	close_ssh_connection(ssh, sftp)

load_lightcurve = load_full_tfalc_from_scratch

twomass_info_for_field = None
def load_2mass_info_for_field():
	print "miscutils: loading twomass_info_for_fields.."
	global twomass_info_for_field
	twomass_info_for_field = { field : twomass_info_file(field) for field in fields_to_analyze }


def add_twomass_info_field(field):
	print "miscutils: loading twomass_info_for_fields for field %s"%(field)
	global twomass_info_for_field
	if twomass_info_for_field is None:
		print "miscutils: twomass info is None."
		twomass_info_for_field = { field : twomass_info_file(field)}
	else:
		twomass_info_for_field[field] = twomass_info_file(field)

def get_bagged_samples(Categories, size, ftest=0.):
	Samples = {}
	for ID in Categories:
		if not Categories[ID] in Samples: Samples[Categories[ID]] = [ ID ]
		else: Samples[Categories[ID]].append(ID)

	bags = []
	test_bag = []
	num_test = { cat : int(ftest)*len(Samples[cat]) for cat in Samples }
	for i in range(size):
		bag = []
		for cat in Samples:
			di = int( 1./float(size) * (len(Samples[cat])-num_test[cat]) )
			i1, i2 = i*di, (i+1)*di
			if i == size - 1:
				bag.extend(Samples[cat][i1:])
			else:
				bag.extend(Samples[cat][i1:i2])
		

		np.random.shuffle(bag)
		bags.append(bag)
	if ftest > 0:
		for cat in Samples:
			test_bag.extend(Samples[cat][-num_test[cat]:])
		np.random.shuffle(test_bag)
	return bags, test_bag

def composite_prediction(ModelPreds, Class, clfrs=None):
	if clfrs is None: 
		#return np.mean([ M[Class] for M in ModelPreds ])
		return max(M[Class] for M in ModelPreds)
	return clfrs[Class].predict_proba(np.ravel(ModelPreds))[0][1] / sum( [ clfr.predict_proba(np.ravel(ModelPreds))[0][1] for clfr in clfrs ] )

class BaggedModel:
	def __init__(self, scalers=None, models=None, w=None):

		# Assumes that either scalers is either None or a list with the same length as models
		self.models = models
		self.scalers = scalers
		self.nclasses_ = None
		self.clfrs_ = None

		# weight the models equally by default
		#self.composite_prediction_ = lambda ModelPreds, Class : np.mean(ModelPreds)

	def get_model_preds_(self, X):
			
		Ys = []
		for s, m in zip(self.scalers, self.models):
			if not s is None: Ys.append(m.predict_proba(s.transform(X)))
			else: Ys.append(m.predict_proba(X))

		Yp = []

		nclasses = len(Ys[0][0])
		if not self.nclasses_ is None: assert(nclasses==self.nclasses_)
		else: self.nclasses_ = nclasses

		for i in range(len(X)):
			Yp.append([ [ Ys[j][i][k] for k in range(nclasses) ] for j in range(len(self.models)) ])

		return Yp

	def predict_proba(self, X):
		
		Ys = self.get_model_preds_(X)

		Y = [ ]
		for i in range(len(X)):
			Y.append([ composite_prediction(Ys[i], k) for k in range(self.nclasses_) ])
		return Y

	def fit(self, X, Y):
		
		ModelPredictions = self.get_model_preds_(X)

		mpreds = []
		for i in range(len(X)):
			mpreds.append( np.ravel(ModelPredictions[i]) )
			
		self.clfrs_ = [ SVC(**svm_params) for i in range(self.nclasses_) ]
		for Label in range(self.nclasses_):
			y = [ 1 if yval == Label else 0 for i,yval in enumerate(Y) ]
			self.clfrs_[Label].fit(mpreds, y)

		#print self.clfrs_[0].predict_proba(mpreds[0])
		self.composite_prediction_ = lambda ModelPreds, Class : \
			self.clfrs_[Class].predict_proba(ModelPreds)[0][1] \
					/ sum( [ clfr.predict_proba(ModelPreds)[0][1] for clfr in self.clfrs_ ] )

	def save(self, root_name):
		filenames = []
		if not self.models is None:
			# Save scalers and models
			for i, m in enumerate(self.models):
				if not self.scalers is None: s = self.scalers[i]
				else: s = None

				fname_s = "%s_s%d.pkl"%(root_name,i)
				fname_m = "%s_m%d.pkl"%(root_name,i)
				fname_sm = "%s_sm%d.pkl"%(root_name,i)
				pickle.dump(s, open(fname_s, 'wb'), HP )
				pickle.dump(m, open(fname_m, 'wb'), HP )
				pickle.dump((fname_s, fname_m), open(fname_sm, 'wb'), HP)

				filenames.append(fname_sm)
			bagged_model_fname = "%s.pkl"%(root_name)
			
			# Save additional classifiers if there are any
			if not self.clfrs_ is None:
				bagged_model_clfrs_fname = "%s_clfrs.pkl"%(root_name)
				clfr_fnames = []
				for clfr in self.clfrs_:
					fname = "%s_clfr%d.pkl"%(root_name,i)
					pickle.dump(clfr, open(fname, 'wb'), HP)
					clfr_fnames.append(fname)
				pickle.dump(clfr_fnames, open(bagged_model_clfrs_fname, 'wb'), HP)

			pickle.dump(filenames, open(bagged_model_fname, 'wb'), HP)

		else:
			raise Exception("Can't save an empty BaggedModel (models == None)")
		return True

	def load(self, root_name):
		bagged_model_fname = "%s.pkl"%(root_name)
		bagged_model_clfrs_fname = "%s_clfrs.pkl"%(root_name)

		filenames = pickle.load(open(bagged_model_fname, 'rb'))
		self.scalers = []
		self.models = []
		
		# Load scalers and models
		for f in filenames:
			fname_s, fname_m = pickle.load(open(f, 'rb'))

			s = pickle.load(open(fname_s,'rb'))
			m = pickle.load(open(fname_m,'rb'))

			self.scalers.append(s)
			self.models.append(m)

		# if there are also additional classifiers, load them too!
		if os.path.exists(bagged_model_clfrs_fname):
			clfr_fnames = pickle.load(open(bagged_model_clfrs_fname, 'rb'))
			self.clfrs_ = [ pickle.load(open(fname, 'rb')) for fname in clfr_fnames ]
			self.nclasses_ = len(self.clfrs_)
			#self.composite_prediction_ = self.composite_prediction

		return True

