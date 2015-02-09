import matplotlib.pyplot as plt
from fastperiod import specwindow, lombscargle
import lsp, sys, os
from os.path import exists
from utils import *
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from math import *
import readhatlc as rhlc
from multiprocessing import Pool
from time import time

#-HAT-079-0000101 RR   -- very strange magnitudes...
#HAT-087-0019610 RRAB -- period seems to be off by a few percent...
#HAT-095-0017322 RRAB -- maybe a double mode? 
#HAT-121-0043209 RRAB -- maybe also a double mode?
#HAT-125-0005338 RRAB -- dm
#HAT-125-0010856 RRAB -- dm
#-HAT-128-0000156 RRC  -- doesnt seem to be a detection, a mag8.5 star, too. Maybe a mismatch?
#HAT-135-0001326 RRC  -- dm??
#-HAT-141-0001285 RRC  -- doesn't seem to be a detection...
#-HAT-141-0004548 RRC  -- no match?
#-HAT-142-0004019 RRAB -- looks more like an RRC if anything -- bright burst at peaks...LOOK AT THIS
#HAT-143-0002021 RRAB -- dm
#HAT-147-0001851 RRAB -- dm
#HAT-147-0008064 RRAB -- looks like it's actually the wrong period
#HAT-149-0000203 RR   -- looks like an RRC, possibly dm
#HAT-149-0000895 RRC  -- dm? or wrong period
#HAT-149-0001139 RRAB -- also looks like the wrong period
#-HAT-150-0012878 RR   -- looks like no detection -- 14.25 mag
#HAT-150-0096348 RRAB -- offset period? (only slightly if so)
#HAT-151-0004188 RRC  -- ok...so this looks strange. 
#							one trough is lower than the other, but that only makes sense for EB/EW/EA
#HAT-151-0012955 RRAB -- dm
#HAT-151-0013885 RRAB -- dm
#HAT-151-0015029 RR   -- dm
#HAT-153-0014745 RRAB -- dm!!
#HAT-165-0015950 RRAB -- dm
#-HAT-168-0002894 RR   -- not sure what's going on here; not a very strong signal?
#HAT-186-0001972 RRAB -- maybe a *slightly* off period
#HAT-188-0001131 RRAB -- dm?
#-HAT-189-0002780 RRAB -- non-detection
#HAT-190-0005199 RRAB -- dm? ok, but low amplitude and strange bumps
#HAT-191-0000875 RRAB -- mayyybe double mode
#HAT-191-0001576 RRAB -- somewhat odd residual (same period)
#HAT-191-0003869 RRAB -- DM
#HAT-191-0008416 RRAB -- DM
#HAT-192-0005095 RRAB -- maybe sliiightly offset period
#HAT-192-0006225 RRAB -- DM
#HAT-192-0007720 RRAB -- possibly double mode
#HAT-192-0008922 RRAB -- possibly double mode
#HAT-192-0009467 RRAB -- DM
#HAT-193-0003992 RRAB -- DM
#HAT-193-0007586 RRAB -- dm?
#HAT-193-0009403 RRAB -- dm
#HAT-194-0003219 RRAB -- dm?
#HAT-194-0011239 RRAB -- OK, but some large outliers!
#-HAT-196-0018339 RRAB -- either non-detection or weak source or contamination (+dm?)
#HAT-197-0034214 RRC  -- dm?
#HAT-197-0038304 RRAB -- slightly offset period?
#HAT-199-1738692 RRAB -- dm?
#HAT-204-0011866 RRAB -- DM
#HAT-204-0018332 RRAB -- DM
#HAT-207-0003276 RRC  -- OK, there are several bright outliers. probably contamination
#-HAT-207-0011053 RRC  -- non-detection
#HAT-213-0001712 RRC  -- slightly offset period?
#HAT-215-0028041 RRAB -- looks like slightly wrong period?
#+HAT-223-0003186 RRAB -- something funky is going on here
#HAT-230-0003941 RRAB -- OK, lots of wild and crazy outliers
#HAT-231-0004186 RRAB -- DM
#HAT-234-0000486 RRAB -- (OK) there seems to be a uniform background? Or maybe DM?
#HAT-235-0001251 RRAB -- (OK) same deal (background)
#HAT-235-0004994 RRAB -- (OK) same deal (background)
#HAT-235-0005901 RRAB -- dm?
#+HAT-237-0002943 RRAB -- uhhh....really wierd shape; reasonably RR-like
#HAT-237-0002973 RR   -- slightly off? or dm?
#HAT-238-0008638 RRAB -- background or dm? (there are more like this that I haven't noted)
#HAT-238-0010246 RRAB -- dm?
#HAT-238-0011016 RRAB -- dm??
#HAT-238-0013573 RRAB -- dm
#+HAT-242-0026174 RRAB -- weak, but probably present, signal (probably wrong period too)
#+HAT-242-0034689 RRAB -- weak, but probably present, signal
#HAT-242-0038489 RRAB -- dm?
#-HAT-248-0000036 RRAB -- weird -- long period (95 d) variability; doesn't look like RR...
#HAT-250-0009555 RRAB -- DM
#+HAT-256-0005695 RRC  -- INVESTIGATE -- don't know what's going on here; LOTS of noise.
#HAT-260-0015176 RR   -- dm?
#HAT-264-0010203 RRC  -- seems like dm
#HAT-268-0003414 RRAB -- dm?
#HAT-268-0006729 RRAB -- dm?
#HAT-269-0004177 RRAB -- dm?
#HAT-270-0002241 RRAB -- dm?
#HAT-270-0004608 RRAB -- dm?
#HAT-276-0006129 RRAB -- dm?
#HAT-277-0003609 RRAB -- dm
#-HAT-277-0004093 RRC  -- non-det?
#HAT-277-0007376 RRAB -- dm? maybe?
#HAT-277-0008287 RRAB -- dm
#HAT-284-0007754 RRAB -- dm
#HAT-285-0011744 RRAB -- dm?
#HAT-285-0013544 RRAB -- dm?
#HAT-286-0008179 RRAB -- dm?
#HAT-286-0017559 RRAB -- dm
#HAT-286-0018039 RRAB -- dm?
#-HAT-287-0017860 RR   -- no detection?
#HAT-287-0018542 RRAB -- dm
#HAT-287-0025460 RRAB -- great double mode
#-HAT-292-0028865 RR   -- non-det
#+HAT-292-0100671 RR   -- weak signal or non-det (looks C-ish if anything)
#+HAT-332-0001158 RRC -- wrong period (short period below pmin)
#+HAT-339-0101490 RRAB -- hmm...either weak signal or no detection; definitely variable, probably RRAB
#-HAT-339-0136924 RR  -- non det
#-HAT-362-0002588 RR -- non det
#+HAT-363-0012214 RRAB -- OK, lots of spurious noise 
#-HAT-388-0000557 RR -- non det
#+HAT-431-0000070 RR -- fundamental mode and 2nd mode are switched
#+HAT-437-0000456 RR -- unfortunate period -- strange double bump at the peak

'''
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
def get_t_x(lc, coltype=COL_TYPE,selection=COL_SELECTION, ttype='BJD'):
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
		brightest_ytype = ytypes[0]
		brightest_mag = np.mean(nancleaned(lc[brightest_ytype]))
		for ytype in ytypes:
			if ytype == ytypes[0]: continue
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

def fitting_function(*args):
	# Fits multiple periods!
	t = args[0]
	nharms = args[1]
	nprs = args[2]
	angular_frequencies = [ args[3 + i] for i in range(nprs) ]
	c = args[3 + nprs]
	o = 4 + nprs
	val = c*np.ones(len(t))
	for i,w in enumerate(angular_frequencies):
		for h in range(nharms):
			A = o + 2*(i*nharms + h)
			B = o + 2*(i*nharms + h) + 1
			n = h + 1
			val += args[A]*np.cos(n*w*t) + args[B]*np.sin(n*w*t)
	return val
def phase_fold(t,x,p,dp=DPHASE):
	nparr = np.zeros(len(t), dtype=np.dtype([('tpf', np.float_),( 'x', np.float_)]))
	T0 = min(t)
	for i, T in enumerate(t):
		nparr[i]['tpf'] = dp*((T - T0)/(dp*p) - int((T - T0)/(dp*p)))
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
def fit_periods(t,x,ps):
	ws = tuple([ 2*np.pi/float(p) for p in ps ])

	ff = lambda T, *args : fitting_function(T, nharmonics, len(ps), *(ws + args) )

	# p0 is the initial guess for each parameter.
	p0 = np.ones(2*nharmonics*len(ws) + 1)
	p0[0] = np.mean(x) # guess a constant offset of <magnitude> (this isn't crucial)

	# fit!
	popt, pcov = curve_fit(ff, t, x, p0=p0)

	
	# popt[0] is the constant offset.
	# then it's (A,B) pairs for each harmonic, for each frequency: 
	#		popt = [ C, A_f1h1, B_f1h1, A_f1h2, B_f1h2, ..., A_f1hNH, B_f1NH, A_f2h1, B_f2h1, ... ]
	As = [ popt[ 1 + 2*(freq_no*nharmonics + harm_no)  ]
				for harm_no in range(nharmonics) 
			for freq_no in range(len(ws)) ]
	Bs = [ popt[ 1 + 2*(freq_no*nharmonics + harm_no) + 1  ]
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
def get_resid(t, x, ws, amps, phs, c):
	return x - np.array([ c + sum([  A*cos(W*T - phi) for A,W,phi in zip(amps, ws, phs) ]) for T in t ])
'''
def plot_fit(t,x,ws,amps,phs,c, show=False, dwdt = 0.):
	dp_plt = 2.0
	# Do Lomb-Scargle
	wk1,wk2,nout,jmax,prob = lsp.fasper(t, x, ofac, hifac, MACC)
	P = np.power(wk1[::-1],-1)
	LSP = wk2[::-1]

	NP = len(ws)/nharmonics
	# Make phase-folded plots of each mode
	WS = [ [ ws[j*NP + i] for j in range(nharmonics) ] for i in range(NP) ]
	AMPS = [ [ amps[j*NP + i] for j in range(nharmonics) ] for i in range(NP) ]
	PHASES = [ [ phs[j*NP + i] for j in range(nharmonics) ] for i in range(NP) ]
	#print 2*np.pi/np.array(ws)
	#for i in range(len(WS)):
	#	print 2*np.pi/np.array(WS[i])
	tpfs = []
	tpf_fit = np.linspace(0,dp_plt,100)
	xpfs = []
	Pers = []
	Pows = []
	xpf_fits = []
	for i in range(len(WS)):
		if len(WS) > 1:
			ws_out = [ ws[j*NP + k] for j in range(nharmonics) for k in range(NP) if k!=i ]
			amps_out = [ amps[j*NP + k] for j in range(nharmonics) for k in range(NP) if k!=i ]
			phs_out = [ phs[j*NP + k] for j in range(nharmonics) for k in range(NP) if k!=i ]
			
			R = get_resid(t, x, ws_out, amps_out, phs_out, c)
			Tpf, Xpf = phase_fold(t, R, 2*np.pi/WS[i][0], dp=dp_plt, dwdt=dwdt)
			wk1,wk2,nout,jmax,prob = lsp.fasper(t, R, ofac, hifac, MACC)
			p = wk1[::-1]**(-1)
			l = wk2[::-1]
			Pers.append(p)
			Pows.append(l)
			tpfs.append(Tpf)
			xpfs.append(Xpf + c)
		else:
			Tpf, Xpf = phase_fold(t, x, 2*np.pi/WS[i][0], dp=dp_plt, dwdt=dwdt)
			wk1,wk2,nout,jmax,prob = lsp.fasper(t, x, ofac, hifac, MACC)
			p = wk1[::-1]**(-1)
			l = wk2[::-1]
			Pers.append(p)
			Pows.append(l)
			tpfs.append(Tpf)
			xpfs.append(Xpf)

		# THIS ISN'T WORKING (it's giving a fitting function that looks reasonable but that is out of phase (but the residual is flat...))
		xpf_fits.append( np.array([ c + sum([  A*cos(2*np.pi*W*T/WS[i][0] - phi) for A,W,phi in zip(AMPS[i], WS[i], PHASES[i]) ]) for T in tpf_fit ]))

	
	#tpru, xpru = prune_outliers(tpf, xpf, [ c + sum([  A*cos(2*np.pi*W*T/ws[0] - phi) for A,W,phi in zip(amps, ws, phs) ]) for T in tpf  ])
	#tpru, xpru = prune_outliers(tpf, xpf)


	resid = get_resid(t, x, ws, amps, phs, c, dwdt = dwdt)
	#tf, xf = phase_fold(t, x - resid, 2*np.pi/ws[0])
	#tpf_fit, xpf_fit, sig_fit = bin_phase_folded(tf, xf, get_errs(xpf))

	# Get LSP of residual
	wk1,wk2,nout,jmax,prob = lsp.fasper(t, resid, ofac, hifac, MACC)
	Pr = wk1[::-1]**(-1)
	LSPr = wk2[::-1]
	best_per_res = 1./wk1[jmax]

	# Phase fold residual
	tpf_r, xpf_r = phase_fold(t, resid, best_per_res, dp=dp_plt, dwdt=dwdt)
	#tpf_r, xpf_r = phase_fold(t, resid,2*np.pi/ws[0] , dp=dp_plt)

	# Plot results

	f = plt.figure(figsize=(4*(2+len(WS)), 8), tight_layout=True)
	#axpfs = []
	#for i in range(len(WS)): axpfs.append(f.add_subplot(2, 1 + npers, i))


	for i in range(len(WS) + 1): 
		if i == len(WS):
			tpf = tpf_r
			xpf = xpf_r
			title = "Residual (P=%.4f d) %.3fP0"%(best_per_res, best_per_res*WS[0][0]/(2*np.pi) )
			ls_title = "Lomb-Scargle (resid)"
			lsperiods = Pr 
			lspowers = LSPr
			lsp_label = "P=%.4f"%(best_per_res)

		else:
			tpf = tpfs[i]
			xpf = xpfs[i]
			if i == 0:
				title = "Mode 1: P=%.4f d"%(2*np.pi/WS[0][0])
			else:
				title = "Mode %d: P=%.4f d (P=%.3fP0)"%(i+1, 2*np.pi/WS[i][0], WS[0][0]/WS[i][0])
			ls_title = "Lomb-Scargle"
			lsp_label = "P=%.4f"%(2*np.pi/WS[i][0])
			lsperiods = Pers[i]
			lspowers = Pows[i]

		ax_phase_folded = f.add_subplot(2, 2 + len(WS), i + 2)

		ax_phase_folded.scatter(tpf, xpf, marker=',', alpha=0.1)
		ax_phase_folded.axvline(1.0, color='k', ls=':')
		ax_phase_folded.set_title(title)
		ax_phase_folded.set_xlim(0,dp_plt)
		if i < len(WS): ax_phase_folded.plot(tpf_fit, xpf_fits[i], color='r', lw=2)
		ax_phase_folded.set_xlabel("Phase")
		ax_phase_folded.set_ylabel("Magnitude")
		ax_phase_folded.invert_yaxis()

		ax_lsp = f.add_subplot(2,2+len(WS), 2+len(WS) + i + 2)
		ax_lsp.plot(lsperiods, lspowers, label=lsp_label)
		ax_lsp.legend(loc='best')
		ax_lsp.set_title(ls_title)
		ax_lsp.set_xlabel("Period [days]")
		ax_lsp.set_ylabel("Power")
		ax_lsp.set_xscale('log')
	#ax_phase_folded.legend(loc='best')
	ax_lsp = f.add_subplot(2,2+len(WS), 1)
	ax_lsp.plot(P, LSP)
	#ax_lsp.legend(loc='best')
	ax_lsp.set_title("Lomb-scargle (FULL)")
	ax_lsp.set_xlabel("Period [days]")
	ax_lsp.set_ylabel("Power")
	ax_lsp.set_xscale('log')
	
	if show: plt.show()
def plot_fit_custom(t,x,tt,xt, label_data = "Data", label_fit="Fit"):
	f = plt.figure()
	
	ax = f.add_subplot(111)
	ax.set_ylabel('Mag')
	ax.set_xlabel('Phase')
	ax.scatter(t,x,marker=',',color='b',alpha=0.1, label=label_data)
	ax.plot(tt,xt,color='r', label=label_theory, lw=2)
	ax.legend(loc='best')
	ax.invert_yaxis()


# Get GCVS cross-matched HAT-IDs

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
ids_to_use = [ hid for hid in gcvs_m if gcvs_m_types[hid] in types_to_use and os.path.exists(hat_fname(hid))]
#last_index = ids_to_use.index("HAT-431-0017178")
#ids_to_use = [ ids_to_use[i] for i in range(last_index, len(ids_to_use)) ]
print len(ids_to_use)
#ids_to_use = [ 'HAT-062-0001360', 'HAT-062-0004746', 'HAT-063-0000158', 'HAT-063-0000487', 
#'HAT-063-0005403', 'HAT-063-0021269', 'HAT-064-0003647', 'HAT-079-0000101','HAT-080-0003698',
#'HAT-081-0003305', 'HAT-081-0010146', 'HAT-081-0012119', 'HAT-081-0012735' ]
#HAT-062-0001360 EA, HAT-062-0004746 EA, HAT-063-0000158 E (?), HAT-063-0000487 EA (right period i think, messed up fit)
#HAT-063-0005403 EW -- period overestimate i think
#HAT-063-0021269 EW -- period overestimate (need larger hifac)
#HAT-064-0003647 EW -- period overestimate
#HAT-079-0000101 RR -- period overestimate? or contamination from window function
#HAT-080-0003698 EB -- period *underestimate* (why is this happening? Not enough peaks?)
#HAT-081-0003305 EA -- period underestimate!! UGH!
#HAT-081-0010146 EA -- period underestimate
#HAT-081-0012119 EW -- period underestimate
#HAT-081-0012735 EA -- period underestimate
#HAT-082-0000154 EA -- I think period underestimate
#HAT-086-0007407 EB -- again, period underestimate
#HAT-086-0008141 EW -- period underestimate
#HAT-086-0008911 EW -- period underestimate
#HAT-087-0012338 EA -- period underestimate
#HAT-087-0018179 EA -- period underestimate and also not a good period.
#HAT-087-0020583 EB -- period underestimate
#HAT-087-0076119 EA -- looks like the completely wrong period
#HAT-088-0008417 EB -- period underestimate
#HAT-088-0015390 EW -- period underestimate
#HAT-088-0032993 EA -- just wrong period i think
#HAT-089-0029508 E  -- wrong period?
#HAT-092-0000033 EB -- ?
#HAT-092-0006183 EB -- period underestimate
#HAT-093-0007032 EW -- period underestimate?
#HAT-093-0048500 E -- period underestimate
#HAT-094-0012791 EW -- period underestimate
#HAT-094-0013710 EA -- period underestimate
#HAT-094-0034989 EA -- period underestimate
#HAT-095-0001631 EB -- period underestimate
#HAT-095-0005509 EW -- period underestimate
#HAT-095-0011797 EA -- period underestimate
#HAT-095-0218850 EW -- period underestimate
#HAT-107-0004403 EW -- period underestimate
#HAT-107-0005897 EW -- period underestimate
#HAT-108-0004206 RRAB -- either a very unfortunate period or wrong period
#HAT-113-0000380 E -- period underestimate
#HAT-113-0002034 EA -- period underestimate
#HAT-113-0009054 EW -- period underestimate
#HAT-113-0011722 EW -- period underestimate
#HAT-114-0002366 EW -- period underestimate
#HAT-115-0002112 EA -- period underestimate


#ids_to_use = [ 'HAT-150-0000235' ]
#ids_to_use = [ 'HAT-433-0037239' ]
#ids_to_use = [ 'HAT-056-0003869']
#ids_to_use = [ 'HAT-432-0034152', 'HAT-164-0004542', 'HAT-192-0001163', 
#'HAT-063-0004468', 'HAT-340-0035449', 'HAT-127-0004752', 'HAT-059-0000847', 
#'HAT-126-0032086', 'HAT-162-0003089', 'HAT-192-0002811']
#ids_to_use = [ 'HAT-192-0001163' ]
#np.random.shuffle(ids_to_use)
#print len(ids_to_use)
# HATID, rms, sj, sk, 
#

#ids_to_use = []
#ids_to_use.append("HAT-264-0010203")

#ids_to_use.extend(available_hatids[:100])
#ids_to_use.extend( get_n_random_available_hatids(100) )



FORCE_REDO = True
#OK_IDS = []
if not  exists("chi2values.dat") or FORCE_REDO:
	chi2_file = open("chi2values.dat",'w')
	for i,HATID in enumerate(ids_to_use):
		print HATID, gcvs_m_types[HATID]
		lc = rhlc.read_hatlc(hat_fname(HATID))

		if lc is None: continue
		if len(lc['TF1']) < 20: continue
		
		detrend(lc)
		
		times, mags = get_t_x(lc)
		if times is None: continue
		#OK_IDS.append(HATID)

		#continue

		resid = np.zeros(len(mags))
		resid[:] = mags[:]

		PERIODS = [ ]
		std_total = np.std(mags)
		
		while len(PERIODS) < npers:
			
			
			t0 = time()

			wk1,wk2,nout,jmax,prob = lsp.fasper(times, resid, ofac, hifac, MACC)
			P = np.power(wk1[::-1],-1)
			LSP = wk2[::-1]

			okinds = [ i for i in range(len(P)) if P[i] < max_per ] 
			P = P[okinds]
			LSP = LSP[okinds]
			Pbest = P[np.argsort(LSP)[::-1][0]]

			# Now, get a list of the top N peaks (and periods within delta*dP of those peaks), 
			# and look for which period gives you the LEAST (local) scatter around the mean
			
			peaks,pows = find_n_peaks(P, LSP, n_peaks)
			dP_lefts = {}
			dP_rights = {}
			dP_lefts[P[0]] = P[1] - P[0]
			dP_rights[P[0]] = P[1] - P[0]
			dP_lefts[P[-1]] = P[-1] - P[-2]
			dP_rights[P[-1]] = P[-1] - P[-2]
	
			for i in range(1, len(P)-1): 
				dP_lefts[P[i]] = P[i] - P[i-1]
				dP_rights[P[i]] = P[i+1] - P[i]

			
			search_pers = [ ]

			
			def is_kind_of_in(val, arr, eps=0.05):

				for a in arr:
					if abs(a-val)/a < eps: return True, a
				return False, None
			'''
			new_peaks = []
			for s in peaks:
				nharms = 1
				while True:
					is_in, val = is_kind_of_in(s*2. , peaks)
					if not is_in: break
					nharms+=1
					s=val
				if nharms == 1: new_peaks.append(s)
			#	print s, nharms
			peaks=new_peaks
			'''
			new_peaks = []
			small = 0.01
			for per in peaks:
				is_ok = True
				for Per in PERIODS:
					if abs(per/Per - int(round(per/Per))) < small or abs(Per/per - int(round(Per/per))) < small: 
						is_ok = False
						break
				if is_ok: new_peaks.append(per)
			peaks = new_peaks
			#print len(peaks)
			for per in peaks:
				for E in np.linspace(-delta_P*dP_lefts[per], delta_P*dP_rights[per], NSEARCH):
					search_pers.append( per+E )
			#print min(peaks), max(peaks)
			#print min(search_pers), max(search_pers)
			print "     seaching..."
			'''
			for p in peaks: 
				I = 0
				while P[I] < p: I+=1
				if I == 0: dP = P[I+1] - P[I]
				else: dP = P[I] - P[I-1]
				search_pers.extend(np.linspace( p - delta_P*dP, p + delta_P*dP, NSEARCH))
			'''
			best_pers, sig_pers = find_best_of_n_periods(times,resid,search_pers,other_periods=PERIODS)
			
			dt = time() - t0
			pinds = np.argsort(np.abs(peaks - best_pers[0]))
			pk = peaks[pinds[0]]
			dPleft = dP_lefts[pk]
			dPright = dP_rights[pk]
			#print np.sort(np.abs(peaks - best_pers[0]))[0]
			dP = pk - best_pers[0]
			if dP < 0: 
				offset = abs(dP)/dPleft
			else:
				offset = abs(dP)/dPright
			
			print "     done (%f s), %.3e fractions of a grid point away!"%(dt, offset)
			

			PERIODS.append(best_pers[0])
			#PERIODS.append(Pbest)
			#print PERIODS
			#get_blazhko(times, mags, PERIODS[0])

			# Now fit the period(s) to the data and update the residual.
			
			ws, amps, phs, c = fit_periods(times, mags, PERIODS, use_dwdt=False)
			
			chi2_old = get_chi2_pf(times,resid,PERIODS[-1])
			std_old = np.std(resid)
			resid[:] = get_resid(times, mags, ws, amps, phs, c)[:]
			std_new = np.std(resid)
			#bps, sps = find_best_of_n_periods(times,resid,search_pers)
			#print best_pers[0], bps[0]
			wk1,wk2,nout,jmax,prob = lsp.fasper(times, resid, ofac, hifac, MACC)
			P = np.power(wk1[::-1],-1)
			LSP = wk2[::-1]
			Pbest = wk1[jmax]**(-1)
			chi2_new = get_chi2_pf(times,resid,Pbest)

			
			print "   PERIOD %d = %.5fd   chi2: %.3e --> %.3e [-dchi2 = %.3e]"%( len(PERIODS),  PERIODS[-1], chi2_old, chi2_new, chi2_old - chi2_new)
			print "                         std : %.3e --> %.3e [-dstd  = %.3e] = %.3f%% (%.3f%% of total rms)"%(std_old, std_new, std_old - std_new, 100*(std_old - std_new)/std_old, 100*(std_old - std_new)/std_total)
		plot_fit(times,mags,ws,amps,phs,c,show=True)


'''
skipped_ids = []
for HID in ids_to_use:
	if not HID in OK_IDS: skipped_ids.append(HID)


print "SKIPPED: ",len(skipped_ids)
for ID in skipped_ids: print "'",ID,"',"
'''




