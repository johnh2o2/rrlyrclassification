# Tests that mpi_pool.py behaves as expected
# when multiple nodes are involved.
# John Hoffman Feb. 9 2015
from sys import exit
import mpi_pool as mpip
import numpy as np
from math import *
import mpi4py.MPI as mpi
import feature_selection as fs
import os
from settings import *
os.environ['DYLD_LIBRARY_PATH'] = DYLD_LIBRARY_PATH
from lcutils.lcutils_config import *
from lcutils.lcutils_processing import get_hatnet_lc_locations, scp_fetch_hatnet_lcs, collect_fetched_lightcurves, consolidate_hatnet_lightcurves
from utils import *
from featureutils import *
from hatbinlc import read_lc as read_binary_lc
import cPickle as pickle
comm = mpip.Comm()
rank = comm.rank
#HATLC_COL_DEFS

colnames = HATLC_COL_DEFS['hn']['tfalc']
dts = {}
for c in colnames:
	dt = TEXTLC_OUTPUT_COLUMNS[c][3]
	dts[c] = dt

dt_hn = np.dtype([ (c, dts[c]) for c in colnames ])

#for c, dt in [ (c, dts[c]) for c in colnames ]:
#	print c, dt
#print dt_hn.names
field_dirs = {
	219 :  "/nfs/phn15/ar1/H/4KRED/4KAP_LC/219"
}
default_lctype = 'tfalc'
rem_dir = "/nfs/phn15/ar1/H/4KRED/4KAP_LC/219"
keylist_dir = "~/2007_hatnet_phot/G219/BASE"
ssh_host_name = "phn1"
#<hatid> <ra> <dec> <jmag> <jerr> <hmag> <herr> <kmag> <kerr> <flags> 
#<Bmag> <Vmag> <Rmag> <Imag> <umag> <gmag> <rmag> <imag> <zmag> <hatfield> <hatfield objid>


#[objectinfo[x] for x in ('vmag','rmag','imag', 'jmag','hmag','kmag')]
#p = mpip.startPool()

def get_keylist(field, keylist_dir = keylist_dir, ssh_host_name=ssh_host_name):
	keylist = {}
	local_fname = "keylist_field%s.txt"%(field)
	os.system("scp %s:%s/keylist.txt %s"%(ssh_host_name, keylist_dir, local_fname))
	klist_data = np.loadtxt(local_fname, dtype=keylist_dt)
	for kl in klist_data:
		keylist[kl['TSTF']] = { }
		for k in keylist_dt.names:
			if k == 'TSTF': continue
			keylist[kl['TSTF']][k] = kl[k]
	return keylist
def get_all_lc_filenames(field_number=219, ext="tfalc"):
	hatid_remote_filenames = {}
	# TODO
def get_remote_lc_fname(hatid, lctype=default_lctype, field_number = None):
	if field_number is None:
		field_number, ID = hatid.split('-')[1:]
	fname = "%s/%s.%s"%(field_dirs[int(field_number)], hatid, lctype)
	return fname
def transfer_lc(remote_fname, local_fname=None, ssh_host_name=ssh_host_name):
	fname = remote_fname.split("/")[-1]
	move=True
	#print fname
	if local_fname == None:
		local_fname = "%s/%s"%(LCCACHE, fname)
		move=False
	command = "scp %s:%s %s"%(ssh_host_name, remote_fname, LCCACHE)
	os.system(command)
	if move: 
		command = "mv %s/%s %s"%(LCCACHE, fname, local_fname)
		os.system(command)
	print local_fname
	assert(os.path.exists(local_fname))#, "scp of %s:%s didn't appear to work..."%(ssh_host_name, remote_fname))

	return local_fname
def load_2mass(field, twomass_remote_dir='~', twomass_fname='colors_field%s.dat', ssh_host_name=ssh_host_name):
	twomass_dict = {}
	twomass_local_fname = twomass_fname%(field)
	if not os.path.exists(twomass_local_fname):
		os.system("scp %s:%s/%s %s"%(ssh_host_name, twomass_remote_dir, twomass_local_fname, twomass_local_fname))
	twomass_data = np.loadtxt(twomass_local_fname, dtype=twomass_dt)
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
	

	lc['cols'] = colnames

	lc['twomassid'] = 0
	lc['hatstations'] = np.unique(lc['STF'])
	lc['filters'] = [ flt for flt in np.unique(lc['FLT']) ]
	return lc
def add_keylist_data(lc, kl):

	cols = [ 'FLT', 'EXP', 'BJD' ]
	for col in cols:
		lc[col] = []

	for i in range(len(lc['TF1'])):
		for col in cols:
			lc[col].append(kl[lc['TSTF'][i]][col])
	for col in cols:
		lc[col] = np.array(lc[col])
	return lc
def load_tfalc(local_fname):
	lc = {}
	for c in colnames: lc[c] = []
	
	with open(local_fname, 'r') as f:
		for line in f:
			data = line.split("#")
			vals = data[0].split()
			if len(vals) < len(colnames): continue
			for i,c in enumerate(colnames):
				lc[c].append(TEXTLC_OUTPUT_COLUMNS[c][3](vals[i]))
	for c in colnames: lc[c] = np.array(lc[c])
	lc['STF'] = []
	lc['frame'] = []
	for tstf in lc['TSTF']:
		station, frame = tstf.split('-')
		lc['STF'].append(station)
		lc['frame'].append(frame)

	lc['frame'] = np.array([ int(f) for f in lc['frame']])
	lc['STF'] = np.array([ int(stf) for stf in lc['STF']])
	return lc
def load_lightcurve(local_fname):
	#header, data = read_binary_lc(local_fname, 100, False)
	#return header, data 
	#data = np.loadtxt(local_fname, dtype=dt_hn)
	return
	#return data


def get_mc_fit_features(features, pcov_file, N=100 ):
	
	assert(os.path.exists(pcov_file))
	pcov = pickle.load(open(pcov_file,'rb'))
	popt = translate_features_to_popt(features)

	X = np.random.multivariate_normal(popt, pcov, N)
	#print X[0]
	##print X.shape, len(popt)
	feats = []
	for x in X:
		fts = translate_popt_to_features(x)
		for c in features:
			if c in fts: continue
			fts[c] = features[c]
		feats.append(fts)
	#print len(feats), feats[0]
	return feats
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
def score_features(features, pcov_file, model_name='v4', N=5000, kind="other"):

	mag_scaler, mag_model, other_scaler, other_model, \
	skip_features, mag_features, vartypes_to_classify, \
	other_keylist, mag_keylist = pickle.load(open("%s/classifier_%s.clfr"%(model_dir, model_name), 'rb'))

	scalers = { 'mag' : mag_scaler, 'other' : other_scaler}
	models  = { 'mag' : mag_model, 'other' : other_model}
	observs = { 'mag' : None, 'other' : None}
	scaler  =   scalers[kind]
	model   =   models[kind]

	feats = get_mc_fit_features(features,pcov_file,N=N)
	Feats = { i : f for i, f in enumerate(feats)  }

	F = CleanFeatures(Feats)
	F = AddCustomFeatures(F)
	
	observs['mag'], observs['other'] = SplitFeatures(F, mag_features)

	observations = make_obs(observs[kind], keylist=other_keylist)

	if not scaler is None:
		observations = scaler.transform(observations)
	scores = model.predict_proba(observations)

	return np.array([ s[1] for s in scores ])
def test_lcload(hatid, tm, keylist):
	fname = transfer_lc("%s/%s.tfalc"%(rem_dir, hatid), local_fname="%s/%s.tfalc"%(data_dir, hatid))
	lc = load_tfalc(fname)
	lc = add_2mass(lc, tm)
	lc = add_keylist_data(lc, keylist)
	lc['hatid'] = hatid
	return lc
	#print dat.columns
	#print dat['TSTF']
	#sys.exit()
#HATIDs = [ "HAT-199-1738692" ]
#def retrieve_unlabeled_hatids
HATIDs = [ 'HAT-161-0020286', 'HAT-285-0013544', 'HAT-235-0001251', 'HAT-388-0060110', 'HAT-079-0000101', 'HAT-173-0000563', 'HAT-199-1738692' ]
#HATIDs = [ 'HAT-079-0000101', 'HAT-173-0000563', 'HAT-199-1738692' ]
#HATIDs = [ 'HAT-161-0020286' ]


if rank == 0:
	#pcov_file = get_pcov_file(HATID)
	#lc = rhlc.read_hatlc(hat_fname(HATID))
	#detrend(lc)
	#t, x = get_t_x(lc)
	#fts = fs.GetPeriodicFeatures(t, x, save_pcov=True, pcov_file=pcov_file)
	'''
	twomass_dict = load_2mass("219")
	keylist = get_keylist("219")
	

	remote_fname =  get_remote_lc_fname(HATID, field_number=219)
	#print remote_fname
	local_fname = transfer_lc(remote_fname)
	lc = load_tfalc(local_fname)
	lc = add_keylist_data(lc, keylist)
	lc = add_2mass(lc, twomass_dict[HATID])
	

	#collected_lc_dict = collect_fetched_lightcurves(HATID, lcdir=None, outdir=None,removefetched=False,saveinfo=False, ignorecollected=True)

	#lcdata, stfs = consolidate_hatnet_lightcurves(collected_lc_dict, removecollected=False)
	#lcdict = process_consolidated_lightcurve(HATID, lcdata, colnames, database=None)
	#for key in lcdict: print key
	F = fs.get_features(lc, detrend_vars = ['STF'], loud=True, save_pcov=True, pcov_file=pcov_file)	
	'''
	Feats = LoadAllFeatures(HATIDs)

	for ID in HATIDs:
		#if Feats[ID] is None: 
		#	print "Skipping", ID
		#	continue
		lc = rhlc.read_hatlc(hat_fname(ID))
		if lc is None:
			print "Skipping ", ID 
			continue
		F = fs.get_features(lc, loud=True, save_pcov=True, pcov_file=get_pcov_file(ID))
		if F is None:
			print "Skipping ", ID 
			continue
		scores = score_features(Feats[ID], pcov_file=get_pcov_file(ID))
		#continue

		#print np.mean(scores), np.std(scores)
		PlotProbHistogram(scores, np.ones(len(scores)), title="Histogram of classifier scores for %s"%(ID), bins=50, logy=True)

		# Maybe save an interpolation function for p(score)
		# + Flag possible RRAB candidates
		# + 

		#print len(scores)
		#print np.mean(scores), np.std(scores)
		#plt.hist(scores, bins=10)
		#plt.show()
	
#p.exit()
