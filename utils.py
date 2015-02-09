import pandas as pd 
import numpy as np
import defaults
from math import *
import defaults
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import zscore as ZSCORE
from fastperiod import specwindow, lombscargle
import sys, os, re
import fastlombscargle as lsp
from os.path import exists
import readhatlc as rhlc

data_dir = "/Users/jah5/Documents/Fall2014_Gaspar/data"
time_col = 'BJD'
grpsize = 10
bad_ids = [ 'HAT-079-0000101', 'HAT-128-0000156', 'HAT-141-0001285', 'HAT-141-0004548', 'HAT-142-0004019'
'HAT-150-0012878', 'HAT-168-0002894', 'HAT-189-0002780', 'HAT-196-0018339', 'HAT-207-0011053', 
'HAT-248-0000036', 'HAT-277-0004093', 'HAT-287-0017860', 'HAT-292-0028865', 'HAT-339-0136924',
'HAT-362-0002588', 'HAT-388-0000557', ' HAT-135-0007139 ', 'HAT-189-0006202', 'HAT-239-0006835','HAT-241-0014081',
'HAT-241-0018480']
look_at_again = [ 'HAT-223-0003186', 'HAT-237-0002943', 'HAT-242-0026174', 'HAT-242-0034689','HAT-256-0005695',
'HAT-292-0100671', 'HAT-332-0001158', 'HAT-339-0101490', 'HAT-363-0012214','HAT-431-0000070', 'HAT-437-0000456' ]
variable_star_classes = {
	'Eruptive' : ['FU', 'GCAS', 'I', 'IA', 'IB', 
			'IN', 'INA', 'INB', 'INT', 'IT', 'IN(YY)', 
			'IS', 'ISA', 'ISB', 'RCB', 'RS', 'SDOR', 'UV', 'UVN', 'WR'],

	'Pulsating' : ['ACYG', 'BCEP', 'BCEPS', 'CEP', 'CEP(B)', 'CW', 'CWA', 'CWB', 'DCEP', 'DCEPS',
           'DSCT', 'DSCTC', 'GDOR', 'L', 'LB', 'LC', 'M', 'PVTEL', 'RPHS', 'RR', 'RR(B)', 'RRAB',
           'RRC', 'RV', 'RVA', 'RVB', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'SXPHE', 'ZZ', 'ZZA', 'ZZB'],

    'Rotating' : [ 'ACV', 'ACVO', 'BY', 'ELL', 'FKCOM', 'PSR', 'SXARI' ],

    'Cataclysmic' : ['N', 'NA', 'NB', 'NC', 'NL', 'NR',
           'SN', 'SNI', 'SNII', 'UG', 'UGSS', 'UGSU', 'UGZ', 'ZAND'],

    'Eclipsing binary' : ['E', 'EA', 'EB', 'EW', 'GS', 'PN', 'RS', 'WD', 'WR', 'AR', 'D', 'DM',
           'DS', 'DW', 'K', 'KE', 'KW', 'SD'],
    'RR Lyrae' : [ 'RRAB', 'RRC', 'R', 'RR', 'RR(B)' ]
}
feat_dir = "/Users/jah5/Documents/Fall2014_Gaspar/features"
hat_features_fname = lambda hatid : "/Users/jah5/Documents/Fall2014_Gaspar/features/%s-features-v3.pkl"%(hatid)

COL_TYPE = 'TF'
COL_SELECTION = 'locally-brightest'

nharmonics = 8
npers      = 1

NPEAKS_TO_SAVE = 5

n_peaks    = 1
DELTA_I    = 10
delta_P    = .01

DPHASE     = 1.0
NSEARCH    = 1

ofac       = 3
hifac      = 2
MACC       = 3
NBINS      = 50
eps        = 10E-7
max_per    = 10.

use_bootstrap = True
boot_size  = 100
n_boots = 1
fraction_of_smallest_stds_to_use = 0.5
'''
# v2
COL_TYPE = 'TF'
COL_SELECTION = 'locally-brightest'

nharmonics = 6
npers      = 2

NPEAKS_TO_SAVE = 5

n_peaks    = 5
DELTA_I    = 10
delta_P    = 1.0

DPHASE     = 1.0
NSEARCH    = 7

ofac       = 4
hifac      = 4
MACC       = 4
NBINS      = 50
eps        = 10E-7
max_per    = 10.

use_bootstrap = True
boot_size  = 500
n_boots = 1
fraction_of_smallest_stds_to_use = 0.5
'''
API_KEY = {
    'lcdirect' : 'ZjNmZjQ0NzY4MTQxNzQ0Zjk1OTdlNzY1MTAxOTY1YTQyNDNlMzZlZmE2MWE3M2E3YTY0OWE1MDM5ZDU5NmRjYQ'
}

datafiles = os.listdir(data_dir)
regexp = re.compile("(HAT-...-.......)-hatlc\.csv\.gz")
available_hatids = []
for df in datafiles:
	result = regexp.match(df)
	if result is None: continue
	elif result.groups(1) is None: continue
	else: hatid = result.groups(1)[0]
	available_hatids.append(hatid)
available_hatids = np.array(available_hatids)
#print available_hatids[0]
#print len(available_hatids)

#def get_n_random_available_hatids(n):


phs13_lcdir = "/nfs/phs3/ar1/lcserver/lightcurves"
phs13_list_dt = np.dtype([
		('hatid', 'S15'),
		('ndet', np.int_),
		('relative_path','S200')
	])

phs13_list_file = 'lcs_on_phs13_gcvs.txt'
phs13_list = {}
with open(phs13_list_file,'r') as f:
	head = 0
	for line in f:
		if head < 2:
			head += 1
			continue
		cols = [ l.strip() for l in line.split('|') ]
		phs13_list[cols[0]] = {
			'ndet' : int(cols[1]),
			'relpath' : cols[2]
			}


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
def is_available(hid):
	fname = hat_fname(hid)
	if not os.path.exists(fname): return False
	try:
		lc = rhlc.read_hatlc(fname)
	except:
		return False
	if lc is None: return False
	return True
def get_n_random_available_hatids(n):
	inds = np.arange(0,len(available_hatids))
	np.random.shuffle(inds)
	return available_hatids[inds[:n]]
def fetch_lcs(hatids):
	i = 0
	lh = len(hatids)
	while i < lh:
		j = i+grpsize
		if j > lh: j = lh
		IDgrp = hatids[i:j]
		filtered_grp = []
		for ID in IDgrp:
			if not os.path.exists(hat_fname(ID)): 
				filtered_grp.append(ID)
		has_all = (len(filtered_grp) == 0)
		if has_all: continue
		idlist=""
		for i,ID in enumerate(filtered_grp):
			idlist = "%s%s"%(idlist,ID)
			if i < len(filtered_grp) - 1: idlist = "%s,"%(idlist)
		url = 'https://hatsurveys.org/lightcurves/lc/direct?hatid=%s&apikey=%s'%(idlist, API_KEY['lcdirect'])
		print url
		fetch_command = "curl -J -O '%s'"%(url)
		unzip_command = "unzip *zip; rm *zip"
		os.system(fetch_command)
		os.system(unzip_command)
		os.system("mv *hatlc*gz %s"%(data_dir))
		i += j
def make_batch_dl_file(hatids):
	script = "#!/bin/bash\n mkdir gcvs_source_lcs\n"
	for hid in hatids:
		if hid not in phs13_list: continue
		script= "%s cp %s/%s gcvs_source_lcs/%s-hatlc.csv.gz\n"%(script, phs13_lcdir, phs13_list[hid]['relpath'], hid)
	script = "%s tar cvf gcvs_source_lcs.tar gcvs_source_lcs/*"%(script)
	with open("dlscript.sh",'w') as sfile:
		sfile.write(script)
def fetch_lc_from_phs13(hatid):
	if hatid in phs13_list:
		remote_fname = "%s/%s"%(phs13_lcdir, phs13_list[hatid]['relpath'])
		command = "scp -i ~/.ssh/id_dsa hat:%s %s/%s-hatlc.csv.gz"%(remote_fname, data_dir, hatid)
		os.system( command )
		return True
	else: return False
def fetch_lcs_from_phs13(hatids):
	for hid in hatids: fetch_lc_from_phs13(hid)
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
			t.append(lc[ttype][i] - lc[ttype][0])
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
			X = lc[brightest_ytype][i]
			if np.isnan(X): continue
			ndatas += 1
			x.append(X)
			t.append(lc[ttype][i] - lc[ttype][0])

	elif selection in [1, 2, 3]:
		# use a specified aperture
		ytype = "%s%d"%(coltype, selection)
		ndatas = 0
		for i in range(len(lc[ttype])):
			X = lc[ytype][i]
			if np.isnan(X): continue
			ndatas += 1
			x.append(X)
			t.append(lc[ttype][i] - lc[ttype][0])

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


def find_n_peaks(periods, powers, n_peaks, dn=defaults.settings['peaks-dn']):
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
def fast_lsp(t, x):
	bslsp = bootstrap_lsp(t, x, n_boots, boot_size)

	#Freqs,Power,nout,jmax,prob 
	pers = [ l[0][::-1]**(-1) for l in bslsp ]
	lsps = [ l[1][::-1] for l in bslsp ]


	DT = get_dt(t)
	min_p = min([ p[0] for p in pers ])
	drange = int(min_p/(DT*ofac))

	while True:
		if drange < boot_size: drange = boot_size

		bslsp_hf = bootstrap_lsp(t, x, n_boots, boot_size, data_range=drange)
		
		pers_hf = [ l[0][::-1]**(-1) for l in bslsp_hf ]
		lsps_hf = [ l[1][::-1] for l in bslsp_hf ]

		pers_new = []
		lsps_new = []
		for i in range(n_boots):

			imax = len(pers_hf[i])-1
			while pers_hf[i][imax] > min_p and imax > 0: imax -= 1
			if imax <= 1: continue

			all_pers = [ p for p in pers_hf[i][:imax+1] ]
			all_lsps = [ l for l in lsps_hf[i][:imax+1] ]
			all_pers.extend(pers[i])
			all_lsps.extend(lsps[i])

			pers[i] = np.array(all_pers)[:]
			lsps[i] = np.array(all_lsps)[:]


		if drange == boot_size: break

		min_p = min([ p[0] for p in pers_hf ])
		drange = int(min_p/(DT*ofac))

	imin = min([ len(P) for P in pers ] )
	#imax = max([ len(P) for P in pers ] )
	#P0min = min([ P[0] for P in pers ])
	#P0max = max([ P[0] for P in pers ])
	#Pfmin = min([ P[-1] for P in pers ])
	#Pfmax = max([ P[-1] for P in pers ])
	#print P0min, P0max
	#print Pfmin, Pfmax
	#print imin, imax
	def eat_away(arr, num):
		if num == 0: return arr
		return arr[num/2:len(arr) - int(ceil(float(num)/2.))]

	#for p in pers:
	#	num =  len(p) - imin

	#	print num, p[:num/2 + 3], p[::-1][:num/2 + 3]
	#	print eat_away(p, num)[:num/2 + 3], eat_away(p, num)[::-1][:num/2 + 3]
	pers = [ eat_away(p, len(p) - imin) for p in pers ]
	lsps = [ eat_away(l, len(l) - imin) for l in lsps ]

	#imin = min([ len(P) for P in pers ] )
	#imax = max([ len(P) for P in pers ] )
	#lmin = min([ len(L) for L in lsps ] )
	#lmax = max([ len(L) for L in lsps ] )
	#print imin, lmin, imax, lmax


	

	P = [ np.mean([ p[i] for p in pers ]) for i in range(len(pers[0])) ]
	L = [ np.mean([ l[i] for l in lsps ]) for i in range(len(pers[0])) ]

	Lerr = [ np.std([ l[i] for l in lsps ]) for i in range(len(pers[0])) ]
	#print len(P), len(L), len(Lerr)
	return np.array(P), np.array(L), np.array(Lerr)
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
#def analyze_frequency evolution

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
def fit_periods(t,x,ps, use_dwdt=False):
	ws = tuple([ 2*np.pi/float(p) for p in ps ])

	if use_dwdt:
		ADD = 2
		ff = lambda T, *args : fitting_function(T, nharmonics, len(ps), *(ws + args) )
	else:
		ADD = 1
		def ff( T, *args ):
			new_args = [ 0. ]
			new_args.extend(args)
			new_args = tuple(new_args)
			return fitting_function(T, nharmonics, len(ps), *(ws + new_args) )
		
	# p0 is the initial guess for each parameter.
	p0 = np.ones(2*nharmonics*len(ws) + ADD)
	p0[0] = np.mean(x) # guess a constant offset of <magnitude> (this isn't crucial)
	if use_dwdt: p0[1] = 0.
	# fit!
	popt, pcov = curve_fit(ff, t, x, p0=p0)

	
	# popt[0] is the constant offset.
	# then it's (A,B) pairs for each harmonic, for each frequency: 
	#		popt = [ C, A_f1h1, B_f1h1, A_f1h2, B_f1h2, ..., A_f1hNH, B_f1NH, A_f2h1, B_f2h1, ... ]
	As = [ popt[ ADD + 2*(freq_no*nharmonics + harm_no)  ]
				for harm_no in range(nharmonics) 
			for freq_no in range(len(ws)) ]
	Bs = [ popt[ ADD + 2*(freq_no*nharmonics + harm_no) + 1  ]
				for harm_no in range(nharmonics) 
			for freq_no in range(len(ws)) ]				  
	


	def phase_shift(A, B): 
		if A == 0 and B > 0: return np.pi/2
		elif A == 0 and B < 0: return -np.pi/2
		elif A < 0:
			return np.arctan(B/A) + np.pi
		else: return np.arctan(B/A)

	angular_freqs = [  (harm_no+1)*ws[freq_no]
						for harm_no in range(nharmonics)
				 	for freq_no in range(len(ws)) ]

	amplitudes = 	[ 	sqrt(A**2 + B**2) 
						for A,B in zip(As, Bs) ]

	phases =     	[ 	phase_shift(float(A),float(B)) 
						for A,B in zip(As, Bs) ]
	if use_dwdt:
		return angular_freqs, amplitudes, phases, popt[0], popt[1]
	else:
		return angular_freqs, amplitudes, phases, popt[0]
def bs_fit_periods(t, x, ps, nsub=boot_size):
	bst, bsx = bootstrapped_lc(t, x, nboots=1, bootsize=nsub)

	bst = bst[0]
	bsx = bsx[0]

	return fit_periods(bst, bsx, ps)
#def binned_stds(t, x, nbins=NBINS):
#	return np.
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
			#tpf, xpf = phase_fold(bt, bx, p)
			#mpf = xpf - avg_subtracted(tpf, xpf)

			#resid = avg_subtracted(tpf, xpf)
			ws, amps, phs, c = fit_periods(bt, bx, tot_pers, use_dwdt=False)
			
			resid = get_resid(bt, bx, ws, amps, phs, c)
			#resid = avg_subtracted(*phase_fold(bst, bsx, p))
			sses.append(np.std(resid))
			

		scores.append(np.mean(sses))
		
	
	'''
	best_p = pers[np.argsort(scores)[0]]
	f = plt.figure()
	ax = f.add_subplot(111)
	ax.axhline(np.sort(scores)[0], color='r')
	ax.set_title(str(best_p))
	ax.axvline(best_p, color='r')
	ax.plot(np.sort(pers), [ scores[i] for i in np.argsort(pers) ])
	#plt.show()
	
	f2 = plt.figure()
	

	pers_sorted = np.array(pers)[np.argsort(scores)]
	
	Nplt = 5

	for i in range(Nplt):
		plt_fit = (Nplt, 2, 2*i + 1)
		plt_res = (Nplt, 2, 2*i + 2)

		ax_fit = f2.add_subplot(*plt_fit)
		ax_res = f2.add_subplot(*plt_res)
		p = pers_sorted[i]

		bst, bsx = bootstrapped_lc(t, x, nboots=1, bootsize=boot_size)

		bst = bst[0]
		bsx = bsx[0]

		ws, amps, phs, c = fit_periods(bst, bsx, [ p ])
		
		resid = get_resid(bst, bsx, ws, amps, phs, c)

		bstpf, bsxpf = phase_fold(bst, bsx, p)
		trpf, rpf = phase_fold(bst, resid, p)

		ax_fit.scatter(bstpf, bsxpf, color='b', marker=',', alpha=0.1)
		ax_res.scatter(trpf, rpf, color='b', marker=',', alpha=0.1)
		ax_res.set_title("P=%.3f, score=%.3f"%(p, np.sort(scores)[i]))
		ax_res.invert_yaxis()
		ax_fit.invert_yaxis()
		ax_res.set_xlim(0,1)
		ax_fit.set_xlim(0,1)

	f2.tight_layout()
	
	plt.show()
	'''
	#plt.show()


	pers_sorted = np.array(pers)[np.argsort(scores)]
	bp = pers_sorted[0]
	#print bp
	#if symmetry_measure(t, x, bp) < symmetry_measure(t, x, 2*bp): return [ bp ], []
	#else: return [ 2*bp ], []
	return [ bp ], []
	#return pers_sorted, np.sort(scores)
#def eclipsing_binary_fitting_function(wt1, wb1, wt1, wb2, d1, d2, phoffset):
def get_blazhko(t, x, p, nps=2):

	sets = []
	S = []
	t0 = t[0]
	for i in range(len(t)):
		if len(S) == 0:
			S.append([t[i], x[i]])
			t0 = t[i]
		elif t[i] - t0 > nps*p:
			sets.append(S)
			S = []
		else:
			S.append([t[i], x[i]])

	As = []
	tims = []
	sig = np.std(x)
	#minN = 100
	for i,S in enumerate(sets):
		ts = [ s[0] for s in S ]
		xs = [ s[1] for s in S ]
		#if len(ts) < minN: continue
		
		try:
			ws, amps, phs, c = fit_periods(ts, xs, [ p ])
		except:
			continue
		TS = np.linspace(i*nps*p, (i+1)*nps*p)
		ms = - get_resid(TS, np.zeros(len(TS)), ws, amps, phs, c)
		A = 0.5*(max(ms) - min(ms))
		if A > 3*sig: continue
		tims.append((i+0.5)*nps*p)
		As.append(A)
	tims = np.array(tims)
	As = np.array(As)

	#print len(tims), len(As)
	wk1,wk2,nout,jmax,prob = lsp.fasper(tims, As, ofac, hifac, MACC)
	LSPr = wk2[::-1]
	P = wk1[::-1]**(-1)
	p0 = 1./wk1[jmax]

	tpf, spf = phase_fold(tims, As, p0)
	f = plt.figure()
	axs = f.add_subplot(121)
	axs.scatter(tpf, spf, color='b', marker=',', alpha=0.1 )
	axs.set_ylabel("scatter")
	axs.set_xlabel("Phase")
	axs.set_xlim(0,1)
	axs.set_title("P=%0.4f (%.2f)"%(p0, (max(t) - min(t))/p0))

	axlsp = f.add_subplot(122)
	axlsp.plot(P, LSPr)
	axlsp.set_xscale('log')
	axlsp.set_ylabel("L-S Power")
	axlsp.set_xlabel("P [days]")

	plt.show()


#print "PHS13 files = %s"%(len(phs13_list))
gcvs_m = []
gcvs_m_types = {}
with open('Full_GCVS_Cross-matched_HATIDS.catalog', 'r') as f:
	for line in f:
		splits = line.split(' ')
		if 'HAT' in splits[0]: gcvs_m.append(splits[0])

		#if 
		else:
			continue
		found_type = False
		for i,s in enumerate(splits):
			if i < len(splits) - 1 and i > 0:
				if s != '': 
					gcvs_m_types[splits[0]] = s
					found_type = True
		if not found_type: gcvs_m_types[splits[0]] = "?"

#types_to_use = [ 'E', 'EW', 'EB', 'EA', 'R', 'RRAB', 'RRC', 'RR', 'RR(AB)' ]
types_to_use = [ 'RRAB', 'RRC', 'RR', 'R']
ids_to_use = [ hid for hid in gcvs_m if gcvs_m_types[hid] in types_to_use and os.path.exists(hat_fname(hid)) and not hid in bad_ids and not hid in look_at_again ]
