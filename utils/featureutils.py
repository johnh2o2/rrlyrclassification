from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import StratifiedKFold
from sklearn.decomposition import RandomizedPCA
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.lda import LDA
from sklearn.qda import QDA
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.grid_search import GridSearchCV
import matplotlib.pyplot as plt
import readhatlc as rhlc
import numpy as np
import os, sys
from settings import *
from scipy.interpolate import interp1d
import feature_selection as fs
from miscutils import *
import cPickle as pickle
from sklearn.metrics import confusion_matrix

def LoadLabeledHatIDs():
	dat = np.loadtxt(get_labeled_hatids_fname(), dtype=dt_labeled_hatids)
	categories = {}
	for d in dat:
		categories[d['ID']] = d['label']
	return dat['ID'], categories

#def LoadUnlabeledHatIDs(): 

def SaveLabeledHatIDs(categories, iteration):
	f = open(get_labeled_hatids_fname(), 'w')
	for ID in categories:
		f.write("%-20s%-10i%-20s\n"%(ID, iteration, categories[ID]))
	f.close()

def GCVS_GetVartypeClasses(vt):
	classes = []
	for vtclass in variable_star_classes:
		if vt in variable_star_classes[vtclass]:  classes.append(vtclass)
	return classes
def GCVS_GetVartypeNameToUse(vt):
	vtclasses = GCVS_GetVartypeClasses(vt)

	if vt in skip_vartypes: return None
	for vtclass in skip_vartypes:
		if vtclass in vtclasses: return None
	
	if classify_categories:
		classes = []
		for vtclass in vtclasses:
			if vtclass in vartypes_to_classify: classes.append(vtclass)
		if len(classes) > 1: raise Exception("Variable type %s falls into several classes; this program assumes binary classification"%(vt))
		elif len(classes) == 0: return 'none'
		else: return classes[0]
	if vt in vartypes_to_classify: return vt
	return "none"
def GCVS_GetRandomSample(num, cat):
	inds = np.arange(len(cat))
	np.random.shuffle(inds)
	num_per_cat = {}
	for c in vartypes_to_classify: num_per_cat[c] = 0
	if num is None:
		full = lambda : False
	else:
		full = lambda : all([ num_per_cat[c] >= num for c in vartypes_to_classify ])
	selected_inds = []
	for i in inds:
		if full(): break

		ID = cat[i]['id']

		vt = cat[i]['vartype']
		#vtclasses = GCVS_GetVartypeClass(vt)

		vtname = GCVS_GetVartypeNameToUse(vt)
		if vtname is None: continue
		if num_per_cat[vtname] < num or num is None:
			num_per_cat[vtname] += 1
			selected_inds.append(i)
	return cat[selected_inds]

def GetHATIDsToFetch(IDs):
	ids_to_fetch = []
	for ID in IDs:
		if not os.path.exists(get_lc_fname(ID)): ids_to_fetch.append(ID)
	return ids_to_fetch
def GetCategoriesForEachHATID(cat):
	categories = {}
	for source in cat:
		categories[source['id']] = GCVS_GetVartypeNameToUse(source['vartype'])
	return categories



def GetGoodHATIDs(IDs, categs ):
	good_ids = []
	i=0
	for ID in IDs:
		i += 1
		#print ID, i, len(IDs), 
		#if not ID in phs13_list: 
		#	print ""
		#	continue
		if ID in bad_ids or ID in look_at_again:				continue
		if categs[ID] not in vartypes_to_classify: 				continue
		if not os.path.exists(get_lc_fname(ID)): 				continue
		if LoadFeatures(ID) is None: 							continue
		
		good_ids.append(ID)
		

		
	return good_ids
def SelectVmagRange(features, Vmin = None, Vmax = None):
	new_features = {}
	Vmags = [ features[ID]['V'] for ID in features ]
	if Vmin is None: Vmin = min(Vmags)
	if Vmax is None: Vmax = max(Vmags)

	for ID in features:
		if features[ID]['V'] < Vmin or features[ID]['V'] > Vmax: continue
		new_features[ID] = {}
		for f in features[ID]:
			new_features[ID][f] = features[ID][f]

	return new_features
def PlotProbHistogram(probs, labels, title="Histogram of classifier scores", bins=50, logy=False):
	#print probs, labels
	fpr, tpr, _ = roc_curve( labels, probs )

	probs_per_class = {}
	for c in label_names: 
		probs_per_class[label_names[c]] = [] 

	for p,l in zip(probs, labels):
		probs_per_class[label_names[l]].append(p)

	f = plt.figure()
	ax = f.add_subplot(111)
	ax.set_title(title)
	colors = [ 'r', 'k', 'b', 'g', 'c' ]

	for i, l in enumerate(np.unique(labels)):
		ax.hist(probs_per_class[label_names[l]], color=colors[i], alpha=0.3, normed=True, label=label_names[l], range=(0,1), bins=bins)
	if logy: ax.set_yscale('log')
	ax.set_ylabel("pdf")
	ax.set_xlabel("Prob score")
	plt.show()
def GetBestColor(x, y):
	xy = zip(x,y)
	np.random.shuffle(xy)
	X, Y = zip(*xy)

	X = StandardScaler().fit_transform(X)
	stds = [  np.std([ f[i] for f in X ]) for i in range(len(X[0]))  ]


	lda = LDA()
	lda.fit(X, Y)

	W = [ l[0] for l in  lda.scalings_]

	components = np.divide(W, stds)
	components /= sqrt(np.dot(components, components))

	return components
def MakeMagnitudeModel(xtrain, ytrain, CLFR=default_clfr ):
	#Scale the data to z-scores.
	scaler = StandardScaler().fit(xtrain)
	xtrain_scaled = scaler.transform(xtrain)
	clfr = CLFR()
	clfr.fit(xtrain_scaled, ytrain)
	
	return scaler , clfr
def MakeOtherModel(xtrain, ytrain):
	clfr =  RandomForestClassifier(**rfc_params)
	clfr.fit(xtrain, ytrain)
	return None, clfr
	#clfr = SVC(**svm_params)
	#scaler = StandardScaler().fit(xtrain)
	#xtrain_scaled = scaler.transform(xtrain)
	#clfr.fit(xtrain_scaled, ytrain)
	#return scaler, clfr
def PlotRandomForestImportances( RFC, Keylist ):
	imps = RFC.feature_importances_
	inds = np.argsort(imps)
	imps = np.sort(imps)
	labels = [ Keylist[i] for i in inds]

	x = np.arange(len(imps)) + 1
	f = plt.figure(figsize=(10,9))
	ax = f.add_subplot(111)

	ax.barh(x - 0.5,imps, height=1.0, alpha=0.5 )
	ax.set_yticks(x)
	ax.set_ylim(min(x) - 0.5, max(x) + 0.5)
	ax.set_yticklabels(labels, fontsize=10)
	f.subplots_adjust(left=0.3, top=0.95, bottom=0.05)
	plt.show()
def interp1d_improved(x, y):
	func = interp1d(x, y)
	def func2(X):
		if X < min(x): return 0
		if X > max(x): return 1
		else: return func(X)
	return func2
def AddROCtoAxis(ax, tprs, fprs, color='r', tpr_range=[0,1], fpr_range=[0,1]):
	FPR = np.linspace(0, 1)
	tpr_funcs = [ interp1d_improved(fpr, tpr) for fpr, tpr in zip(fprs, tprs) ]
	
	tpr_mean = np.array([ np.mean([ t(fpr) for t in tpr_funcs ]) for fpr in FPR ])
	tpr_sig = np.array([ np.std([ t(fpr) for t in tpr_funcs ]) for fpr in FPR ])
	
	ax.set_ylabel("True Positive Rate")
	ax.set_xlabel("False Positive Rate")
	ax.plot(FPR, tpr_mean, color=color, lw=2)
	ax.fill_between(FPR, tpr_mean - 0.5*tpr_sig, tpr_mean + 0.5*tpr_sig, facecolor=color, alpha=0.3)
	ax.set_xlim(fpr_range[0],fpr_range[1])
	ax.set_ylim(tpr_range[0],tpr_range[1])
def PlotMagnitudeDist(features):
	RR_Lyrae_Vmags = [ features[ID]['V'] for ID in features if categories[ID] == 'RR Lyrae' ]
	All_Vmags = [ features[ID]['V'] for ID in features ]
	f = plt.figure()
	ax = f.add_subplot(111)
	ax.hist(All_Vmags, bins=20, normed=True, color='k', alpha=0.3, range=(0,16))
	ax.hist(RR_Lyrae_Vmags, bins=20, normed=True, color='r', alpha=0.3, range=(0, 16))
	ax.set_xlabel("V mag")
	ax.set_ylabel("pdf")
	#plt.show()
def VmagSplit(features, nbins=3):
	RR_Lyrae_Vmags = [ features[ID]['V'] for ID in features if categories[ID] == 'RR Lyrae' ]
	RR_Lyrae_Vmags = np.sort(RR_Lyrae_Vmags)
	bin_width = len(RR_Lyrae_Vmags)/nbins

	Vmins, Vmaxs, Vavgs = [], [], []
	for n in range(nbins):
		Vmins.append(RR_Lyrae_Vmags[n*bin_width])
		if n == nbins - 1: 
			Vmaxs.append(RR_Lyrae_Vmags[-1])
			Vavgs.append(np.mean(RR_Lyrae_Vmags[n*bin_width:]))
		else: 
			Vmaxs.append(RR_Lyrae_Vmags[(n+1)*bin_width])
			Vavgs.append(np.mean(RR_Lyrae_Vmags[n*bin_width:(n+1)*bin_width]))
	return Vmins, Vmaxs, Vavgs
def EvalModel(Xtrain, Xtest, scaler, clfr):
	if scaler is None:
		PPtrain, PPtest = clfr.predict_proba(Xtrain), clfr.predict_proba(Xtest)
	else:
		PPtrain, PPtest = clfr.predict_proba(scaler.transform(Xtrain)), clfr.predict_proba(scaler.transform(Xtest))
	
	TestProbs = np.array([ p[1] for p in PPtest ])
	TrainProbs = np.array([ p[1] for p in PPtrain ])
	
	return TrainProbs, TestProbs
def AddCustomFeature(features):
	new_features = {}
	for f in features:
		new_features[f] = features[f]
	for p in range(1,npers+1):
		for h in range(1, nharmonics+1):
			new_features['p%d R%d/R1'%(p, h)] = features['p%dh%d_amplitude'%(p, h)]/features['p%dh1_amplitude'%(p)]
			#new_features['normed_p%dh%d_phase'%(p, h)] = features['p%dh%d_phase'%(p, h)] = features['p1h1_phase']
	#for p in range(1, 6):
		#new_features['raw_lsp_peak%d_power'%(p)]/=features['raw_lsp_peak1_power']
		#new_features['resid_lsp_peak%d_power'%(p)]/=features['resid_lsp_peak1_power']
		#new_features['p1_lsp_peak%d_power'%(p)]/=features['p1_lsp_peak1_power']
		#new_features['p2_lsp_peak%d_power'%(p)]/=features['p2_lsp_peak1_power']
	return new_features
def AddCustomFeatures(features):
	new_features = {}
	for ID in features:
		new_features[ID] = AddCustomFeature(features[ID])
	return new_features
def LoadFeatures(ID):
	#print "IN LOAD FEATURES", comm.size, comm.rank
	logprint("  featureutils (LoadFeatures): loading features for hatid: %s"%(ID), all_nodes=True)

	feat_fname = hat_features_fname(ID, model_prefix)
	if os.path.exists(feat_fname) and not overwrite:
		try:
			return pickle.load(open(feat_fname, 'rb'))

		except: 
			raise Exception("Cannot load features for HAT-ID %s (filename: '%s')"%(ID, feat_fname))
	else:
		logprint("                             : no features found; generating them!")
		LC = load_lightcurve(ID)
		if LC is None: 
			logprint("                             : LC is NONE :(")
			return None
		features = fs.get_features(LC, save_pcov=True, pcov_file=get_pcov_file(ID))
		pickle.dump(features, open(feat_fname, 'wb'))
		return features
def LoadAllFeatures(IDs):
	Features = {}
	for ID in IDs:
		Features[ID] = LoadFeatures(ID)
	return Features
def CleanFeatures(feats):
	cleaned_feats = {}
	for ID in feats:
		if feats[ID] is None: continue
		cleaned_feats[ID] = {}

		for f in feats[ID]:
			if f in skip_features: continue
			cleaned_feats[ID][f] = feats[ID][f]
	return cleaned_feats
def SplitFeatures(feats, split_feats):
	if split_feats is None: return None, feats
	feats1, feats2 = {},{}
	for ID in feats:
		feats1[ID], feats2[ID] = {},{}
		for f in feats[ID]:
			if f in split_feats:
				feats1[ID][f] = feats[ID][f]
			else:
				feats2[ID][f] = feats[ID][f]
	return feats1, feats2
def MakeObservations(feats, classes):
	IDs, Keylist, Observations, Labels = [], [], [], []

	ids = [ ID for ID in feats ]
	Keylist = [ f for f in feats[ids[0]] ]
	#print integer_labels, classes
	for ID in feats:
		IDs.append(ID)
		#print classes[ID]
		Labels.append(integer_labels[classes[ID]])
		obs = []
		for k in Keylist:
			obs.append(feats[ID][k])
		Observations.append(obs)
	return np.array(IDs), np.array(Keylist), np.array(Observations), np.array(Labels)



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

def score_features(features, pcov_file, iteration=0, N=5000, kind="other"):

	#mag_scaler, mag_model, other_scaler, other_model, \
	other_rootname, mag_rootname,\
	skip_features, mag_features, vartypes_to_classify, \
	other_keylist, mag_keylist = pickle.load(open(get_classifier_fname(iteration), 'rb'))



	#scalers = { 'mag' : mag_scaler, 'other' : other_scaler}
	models  = { 'mag' : BaggedModel(), 'other' : BaggedModel()}
	
	models['mag'].load(mag_rootname)
	models['other'].load(other_rootname)

	observs = { 'mag' : None, 'other' : None}
	#scaler  =   scalers[kind]
	model   =   models[kind]

	feats = get_mc_fit_features(features,pcov_file,N=N)
	Feats = { i : f for i, f in enumerate(feats)  }

	F = CleanFeatures(Feats)
	F = AddCustomFeatures(F)
	
	observs['mag'], observs['other'] = SplitFeatures(F, mag_features)

	observations = make_obs(observs[kind], keylist=other_keylist)

	#if not scaler is None:
	#	observations = scaler.transform(observations)
	scores = model.predict_proba(observations)

	return np.array([ s[1] for s in scores ])

