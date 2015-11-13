mutils = None
import numpy as np
from math import *
import os
from fastlombscargle import fasper as lspraw
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

nharmonics=8

lsp = lambda x, y : lspraw(np.array(x), np.array(y), 3, 200, 2)

def phase_shift(A, B): 
	if A == 0 and B > 0: return np.pi/2
	elif A == 0 and B < 0: return -np.pi/2
	if A == 0 and B == 0: return 0.
	elif A < 0:
		return np.arctan(B/A) + np.pi
	else: return np.arctan(B/A)

def convert_to_args(period, amplitudes, phases, C):
	args = [ period, C ]
	for a,p in zip(amplitudes, phases):
		args.append(a*cos(p))
		args.append(a*sin(p))

	return args

def convert_from_args(*args):
	period = args[0]
	C = args[1]
	assert(not np.isnan(period))
	assert(not np.isnan(C))
	amplitudes = []
	phases = []

	for i in range(2, len(args), 2):
		a = args[i]
		b = args[i+1]
		assert(not np.isnan(phase_shift(a, b)))
		phases.append(phase_shift(a, b))
		assert(not np.isnan(sqrt(a**2 + b**2)))
		amplitudes.append(sqrt(a**2 + b**2))


	assert(len(amplitudes) == len(phases))


	return period, amplitudes, phases, C

def fit_function(times, *args):
	period, amplitudes, phases, C = convert_from_args(*args)

	assert(not np.isnan(period))
	for a in amplitudes: assert(not np.isnan(a))
	for p in phases: assert(not np.isnan(p))

	ns = np.arange(1, len(amplitudes) + 1)
	y =  C + np.sum( [ A*np.cos(2*pi*times*N/period - PHI) for A, N, PHI in zip(amplitudes, ns, phases) ], axis=0 )

	for Y in y: assert(not np.isnan(Y))
	return np.array(y)

def get_fourier_components(lc):
	period = lc.best_period
	p0 = np.zeros(nharmonics*2 + 1)

	popt, pcov = curve_fit(lambda t, *args : fit_function(t, period, *args), lc.times, lc.mags, p0=p0)

	args = (period,) + tuple(popt)
	return convert_from_args(*args)



class Lightcurve:
	def __init__(self, *args, **kwargs):
		self.info = None
		pass

	def load(self, filename):
		pass

	def set_times(self, times):
		self.times = times

	def set_mags(self, mags):
		self.mags = mags

	def set_info(self, info):
		if self.info is None:
			self.info = { key : info[key] for key in info }

class HATLightcurve(Lightcurve):
	def __init__(self, *args, **kwargs):
		self.magcol = 'TF1'
		self.timecol = 'BJD'

	def load(self, hatid):
		lc = mutils.load_full_tfalc_from_scratch(hatid, twomass_dat=tmdat(hatid), save_full_lc=True, force_redo=False)
		self.info = { key : lc[key] for key in lc if not key in [ self.timecol, self.magtype ]}
		t0 = min(lc[self.timecol])

		self.times = [ t - t0 for t in lc[self.timecol] ]
		self.mags = [ m for m in lc[self.magtype] ]

def clean_stars_file(fname):
	stars_raw = open(fname, 'r')
	stars_cleaned = open("%s.cleaned"%(fname),'w')
	indices1 = np.arange(25, 30)
	indices2 = np.arange(31, 36)

	nanstr = ' nan '
	rawdat = stars_raw.read()

	newdat = rawdat.replace('     ', ' nan ')
	lines = newdat.split('\n')
	for line in lines:
		stars_cleaned.write("%s\n"%(line[:47]))

	stars_cleaned.close()
	stars_raw.close()

def ParseM5Data(dir_name, max_num = None):
		cols_table1 = [ 	('VName'	, 	'S3'), 
								('HJD'		,	np.float_), 
								('Bmag'		, 	np.float_), 
								('r_Bmag'	, 	'S3') 
							]
		cols_table3 = [ 	('HJDo'		, 	np.float_), 
								('mago'		,	np.float_), 
								('r_mago'	, 	'S3'), 
								('mag'		,	np.float_), 
								('HJD'		,	np.float_), 
								('VName'	, 	'S3') 
							]
		cols_stars =  [ 		('VName'	, 	'S3'), 
								('RAh'		, 	np.float_), 
								('RAm'		, 	np.float_), 
								('RAs'		, 	np.float_), 
								#('DE'		, 	'S1'), 
								('DEd'		, 	np.float_),
								('DEm'		, 	np.float_),
								('DEs'		, 	np.float_),
								('Bmag'		, 	np.float_),
								('Vmag'		, 	np.float_),
								('Pa'		, 	np.float_),
							]

		dt_table1 =  np.dtype(cols_table1)
		dt_table3 =  np.dtype(cols_table1)
		dt_stars = np.dtype(cols_stars)

		# Missing data in the stars file is just blank spaces >:(
		clean_stars_file(os.path.join(dir_name, "stars.dat"))

		fname_stars = os.path.join(dir_name, "stars.dat.cleaned")
		fname_table1 = os.path.join(dir_name, "table1.dat")
		#fname_table3 = os.path.join(dir_name, "table3.dat")

		stars_dat = np.loadtxt(fname_stars, dt_stars)
		table1_dat = np.loadtxt(fname_table1, dt_table1)
		#table3_dat = np.loadtxt(fname_table3, self.dt_table3)

		lcnames = [ s['VName'] for s in stars_dat ]
		lightcurves = [ ]
		i=0
		for name in lcnames:
			if not max_num is None and i > max_num: break
			i = [ j for j in range(len(stars_dat)) if stars_dat['VName'][j] == name ][0]

			star_info = { col : stars_dat[i][col] for col in stars_dat.dtype.names }
			times = [ dat['HJD'] for dat in table1_dat if dat['VName'] == name ]
			times = [ t - times[0] for t in times ]
			mags = [ dat['Bmag'] for dat in table1_dat if dat['VName'] == name ]

			lc = Lightcurve()
			lc.set_info(star_info)
			lc.set_times(times)
			lc.set_mags(mags)

			lightcurves.append(lc)
			i += 1
		return lightcurves

if __name__ == '__main__':
	print "Reading M5 RR Lyrae lightcurves .."
	m5lcs = ParseM5Data("/Users/jah5/Documents/Fall2014_Gaspar/rrlyr_classification/information/J_MNRAS_411_1744", 5)

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

		lc.best_fit_times = np.linspace(0, max(lc.times), 10000)

		lc.best_fit = fit_function(lc.best_fit_times, *args)

		lc.best_fit_phases = np.linspace(0, 1, 100)
		lc.best_fit_pf = fit_function(lc.best_fit_phases * period, *args)

	for lc in m5lcs:
		f = plt.figure()

		ax = f.add_subplot(111)
		ax.plot(lc.info['lsp']['periods'], lc.info['lsp']['powers'])

		ax.set_xscale('log')
		ax.set_ylabel("LSP power")
		ax.set_xlabel("period")
		f = plt.figure(figsize=(10,4))

		ax = f.add_subplot(111)
		ax.scatter(lc.phases, lc.mags, marker=',', alpha=0.1)
		
		ax.plot(lc.best_fit_phases, lc.best_fit_pf,lw=2)
		ax.set_xlim(0, 1)
		ax.set_ylabel("mag")
		ax.invert_yaxis()
		ax.set_xlabel("phase")
		plt.show(block=True)