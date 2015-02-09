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
from scipy.interpolate import interp1d
import feature_selection as fs 
from utils import *
import cPickle as pickle
from sklearn.metrics import confusion_matrix

match_cat_dt = np.dtype([
		('id', 'S15'),
		('vartype', 'S10'),
		('dist', np.float_)
	])
num = None
min_ndets = 20
nfolds = 20
#vartypes_to_classify = [ 'M', 'LB', 'EA', 'SR', 'EW', 
#			'RRAB', 'BY', 'UV', 'EA/SD', 'SRB' , 'EB', 
#			'RS', 'L', 'SRA', 'DSCTC', 'E', 'RRC', 
#			'ACV', 'DSCT', 'SRD' , 'INT', 'RR', 'INS',
#			'CST', 'DCEP', 'BE', 'IB', 'EA/DM', 'none']
#vartypes_to_classify = [ 'RRAB', 'none', 'RRC', 'EA', 'EW', 'EB' ]
#vartypes_to_classify = [ 'Pulsating', 'Eclipsing binary', 'none' ]
vartypes_to_classify = [ 'RR Lyrae','none' ]
skip_vartypes = [ '*' ]# , 'R' , 'RR', 'R(B)']
classify_categories = True
#show_confusion_matrix= False
#default_clfr = lambda : SVC(probability=True, kernel='rbf')
default_clfr = QDA
#vartypes_to_classify = [ 'LB', 'EA' ]
#skip_features = [ 'vmag' ]
# 'V', 'chi2_pf', p1_chi2_pf', 'p2_chi2_pf', 'std', 'median_absolute_deviation','p1_frac_of_rms_explained', 'p2_frac_of_rms_explained',  
#'resid_chi2_pf','raw_lsp_peak1_period', 'raw_lsp_peak2_period', 'raw_lsp_peak3_period', 'raw_lsp_peak4_period', 'raw_lsp_peak5_period',
#'resid_lsp_peak1_period', 'resid_lsp_peak2_period', 'resid_lsp_peak3_period', 'resid_lsp_peak4_period', 'resid_lsp_peak5_period',
#
# 'chi2_pf', 'p1_chi2_pf', 'p2_chi2_pf', 'resid_chi2_pf', 'std', 
#				'median_absolute_deviation','p1_frac_of_rms_explained', 'p2_frac_of_rms_explained',  
#				'raw_lsp_peak1_power', 'raw_lsp_peak2_power', 'raw_lsp_peak3_power', 'raw_lsp_peak4_power', 'raw_lsp_peak5_power',
#				'resid_lsp_peak1_power', 'resid_lsp_peak2_power', 'resid_lsp_peak3_power', 'resid_lsp_peak4_power', 'resid_lsp_peak5_power',
#				'p1_lsp_peak1_power', 'p1_lsp_peak2_power', 'p1_lsp_peak3_power', 'p1_lsp_peak4_power', 'p1_lsp_peak5_power',
#				'p2_lsp_peak1_power', 'p2_lsp_peak2_power', 'p2_lsp_peak3_power', 'p2_lsp_peak4_power', 'p2_lsp_peak5_power',
skip_features = [ 'ra', 'dec', 'V', 'std', 
				'raw_lsp_peak1_stetson_strlen', 'raw_lsp_peak2_stetson_strlen', 'raw_lsp_peak3_stetson_strlen', 'raw_lsp_peak4_stetson_strlen',
				'raw_lsp_peak5_stetson_strlen', 
				'raw_lsp_peak1_dworetsky_strlen', 'raw_lsp_peak2_dworetsky_strlen', 'raw_lsp_peak3_dworetsky_strlen', 'raw_lsp_peak4_dworetsky_strlen',
				'raw_lsp_peak5_dworetsky_strlen',
				'resid_lsp_peak1_stetson_strlen', 'resid_lsp_peak2_stetson_strlen', 'resid_lsp_peak3_stetson_strlen', 'resid_lsp_peak4_stetson_strlen',
				'resid_lsp_peak5_stetson_strlen', 
				'resid_lsp_peak1_dworetsky_strlen', 'resid_lsp_peak2_dworetsky_strlen', 'resid_lsp_peak3_dworetsky_strlen', 'resid_lsp_peak4_dworetsky_strlen',
				'resid_lsp_peak5_dworetsky_strlen',
				'p1_lsp_peak1_dworetsky_strlen', 'p1_lsp_peak2_dworetsky_strlen', 'p1_lsp_peak3_dworetsky_strlen', 'p1_lsp_peak4_dworetsky_strlen', 
				'p1_lsp_peak5_dworetsky_strlen', 
				'p1_lsp_peak1_stetson_strlen', 'p1_lsp_peak2_stetson_strlen', 'p1_lsp_peak3_stetson_strlen', 'p1_lsp_peak4_stetson_strlen', 
				'p1_lsp_peak5_stetson_strlen', 
				'p2_lsp_peak1_dworetsky_strlen', 'p2_lsp_peak2_dworetsky_strlen', 'p2_lsp_peak3_dworetsky_strlen', 'p2_lsp_peak4_dworetsky_strlen', 
				'p2_lsp_peak5_dworetsky_strlen', 
				'p2_lsp_peak1_stetson_strlen', 'p2_lsp_peak2_stetson_strlen', 'p2_lsp_peak3_stetson_strlen', 'p2_lsp_peak4_stetson_strlen', 
				'p2_lsp_peak5_stetson_strlen', 
				   ]  
mag_features = [ 'R-V', 'I-V', 'J-V', 'H-V', 'K-V' ]
#mag_features = [ 'V' ]
#mag_features = [ 'vmag' ]
#mag_features = [  'R-V', 'I-V', 'J-V', 'H-V', 'K-V', 'p1_total_amplitude', 'p1h1_amplitude', 'p1h1_phase', 'p1h2_amplitude', 'p1h2_phase' ]
overwrite = False
svm_params = dict(
	kernel = 'rbf',
	class_weight = 'auto',
	probability = True
)

rfc_params = dict(
	n_estimators=50, 
	criterion='gini', 
	max_depth=None, 
	min_samples_split=2, 
	min_samples_leaf=1, 
	max_features='auto', 
	bootstrap=True, 
	oob_score=False, 
	n_jobs= -1, 
	random_state=None, 
	verbose=0, 
	min_density=None, 
	compute_importances=None
)
gridsearch_params = {
	"estimator__gamma": [0.0001, 0.001, 0.01, 0.1, 1.0],
	"estimator__C": [1,2,4,8],
	"estimator__kernel": ["poly","rbf"],
	"estimator__degree":[1, 2, 3, 4],
}
full_gcvs_match_cat_name = "Full_GCVS_Cross-matched_HATIDS_maxdist1.00.catalog"
integer_labels = {'none' : 0}
label_names = { 0 : 'none' }

for t in vartypes_to_classify:
	if t == 'none' : continue
	I = max([ integer_labels[l] for l in integer_labels ])
	integer_labels[t] = I+1
	label_names[I+1] = t

def AddCustomFeature(features):
	new_features = {}
	for f in features:
		new_features[f] = features[f]
	for p in range(1,3):
		for h in range(1, 7):
			new_features['p%d R%d/R1'%(p, h)] = features['p%dh%d_amplitude'%(p, h)]/features['p%dh1_amplitude'%(p)]
			#new_features['normed_p%dh%d_phase'%(p, h)] = features['p%dh%d_phase'%(p, h)] = features['p1h1_phase']
	for p in range(1, 6):
		#new_features['raw_lsp_peak%d_power'%(p)]/=features['raw_lsp_peak1_power']
		new_features['resid_lsp_peak%d_power'%(p)]/=features['resid_lsp_peak1_power']
		new_features['p1_lsp_peak%d_power'%(p)]/=features['p1_lsp_peak1_power']
		new_features['p2_lsp_peak%d_power'%(p)]/=features['p2_lsp_peak1_power']
	return new_features
def AddCustomFeatures(features):
	new_features = {}
	for ID in features:
		new_features[ID] = AddCustomFeature(features[ID])
	return new_features
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
		if not os.path.exists(hat_fname(ID)): ids_to_fetch.append(ID)
	return ids_to_fetch
def GetCategoriesForEachHATID(cat):
	categories = {}
	for source in cat:
		categories[source['id']] = GCVS_GetVartypeNameToUse(source['vartype'])
	return categories
def GetGoodHATIDs(IDs, categs ):
	good_ids = []
	for ID in IDs:
		if not ID in phs13_list: continue
		if categs[ID] not in vartypes_to_classify: continue
		if os.path.exists(hat_fname(ID)) \
			and phs13_list[ID]['ndet'] > min_ndets \
			and not LoadFeatures(ID) is None\
			and not ( ID in bad_ids or ID in look_at_again): 
			#and feats[ID]['V'] < 11.5
				good_ids.append(ID)
	return good_ids
hat_features_fname = lambda hatid : "/Users/jah5/Documents/Fall2014_Gaspar/work/saved_features/%s-features-v2.pkl"%(hatid)

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

def LoadFeatures(ID):
	feat_fname = hat_features_fname(ID)
	if os.path.exists(feat_fname) and not overwrite:
		try:
			return pickle.load(open(feat_fname, 'rb'))
		except: 
			raise Exception("Cannot load features for HAT-ID %s (filename: '%s')"%(ID, feat_fname))
	else:
		LC = rhlc.read_hatlc(hat_fname(ID))
		features = fs.get_features(LC)
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
		cleaned_feats[ID] = {}
		for f in feats[ID]:
			if f in skip_features: continue
			cleaned_feats[ID][f] = feats[ID][f]
	return cleaned_feats
def SplitFeatures(feats, split_feats):
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
	for ID in feats:
		IDs.append(ID)
		Labels.append(integer_labels[classes[ID]])
		obs = []
		for f in feats[ID]:
			obs.append(feats[ID][f])
		Observations.append(obs)
	return np.array(IDs), np.array(Keylist), np.array(Observations), np.array(Labels)
def PlotMagnitudeModelResults(probs, labels):
	#print probs, labels
	fpr, tpr, _ = roc_curve( labels, probs )

	probs_per_class = {}
	for c in label_names: 
		probs_per_class[label_names[c]] = [] 

	for p,l in zip(probs, labels):
		probs_per_class[label_names[l]].append(p)

	f = plt.figure()
	ax = f.add_subplot(111)
	ax.set_title("classifier using color values")
	colors = [ 'r', 'k', 'b', 'g', 'c' ]
	for i, vt in enumerate(vartypes_to_classify):
		ax.hist(probs_per_class[vt], color=colors[i], alpha=0.3, normed=True, label=vt, range=(0,1), bins=20)
	
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


#for f in Features[hatids[0]]: print f
#def CheckFeatures(ID, feats):
#	LC = rhlc.read_hatlc(hat_fname(ID))
#	detrend(LC)
#	t, x = get_t_x(LC)
#	tpf, xpf = ge

# Print out all of the variable types; so you can figure out which ones to skip!
#gcvs_match_catalog = np.loadtxt(full_gcvs_match_cat_name, dtype=match_cat_dt)
#vts = [ g['vartype'] for g in gcvs_match_catalog ]
#print np.unique(vts)

gcvs_crossmatches 					= np.loadtxt(full_gcvs_match_cat_name, dtype=match_cat_dt)
ids_to_fetch 						= GetHATIDsToFetch(gcvs_crossmatches['id'])

categories 							= GetCategoriesForEachHATID(gcvs_crossmatches)

print "Pruning out bad ID's"
hatids 								= GetGoodHATIDs(gcvs_crossmatches['id'], categories)
'''
mean_mags = {}
for ID in hatids:
	lc = rhlc.read_hatlc(hat_fname(ID))
	t, x = get_t_x(lc)
	mean_mags[ID] = np.mean(x)

pickle.dump(mean_mags, open("mean_mags.dict", 'wb'))
'''
mean_mags = pickle.load(open("mean_mags.dict", 'rb'))


inds = [ i for i in range(len(gcvs_crossmatches)) if gcvs_crossmatches[i]['id'] in hatids ]
gcvs_crossmatches 					= GCVS_GetRandomSample(num, gcvs_crossmatches[inds])

print "TOTAL HATID's SELECTED: %d"%(len(hatids))

for vclass in vartypes_to_classify: 
	print vclass, sum([ 1 for ID in hatids if categories[ID] == vclass ])

print "Getting features..."
Full_Features 						= LoadAllFeatures(hatids)
#PlotMagnitudeDist(Full_Features)
obs_mag = {}
for ID in mean_mags:
	obs_mag[ID] = { 'V' : mean_mags[ID] }
PlotMagnitudeDist(obs_mag)
PlotMagnitudeDist(Full_Features)
plt.show()
'''
Vmins, Vmaxes, Vavgs = VmagSplit(Full_Features,nbins=2)

I = 0
print Vmins[I], Vmaxes[I], Vavgs[I]
Full_Features = SelectVmagRange(Full_Features, Vmin=Vmins[I], Vmax=Vmaxes[I])
nrr = 0
nnone = 0
for ID in Full_Features:
	if categories[ID] == 'RR Lyrae': nrr +=1
	else: nnone+=1
print nnone, nrr
'''
print "Cleaning features..."
Features 							= CleanFeatures(Full_Features)
Features                            = AddCustomFeatures(Features)
MagFeatures, OtherFeatures 			= SplitFeatures(Features, mag_features)

#flist = [ f for f in OtherFeatures[hatids[0]] ]
#print flist, len(flist)
# Convert to observations.

if not MagFeatures is None: 
	MagIDs, MagKeylist, MagObservations, MagLabels = MakeObservations(MagFeatures, categories)
OtherIDs, OtherKeylist, OtherObservations, OtherLabels = MakeObservations(OtherFeatures, categories)


cols = GetBestColor(MagObservations, MagLabels)
LDA_ColorValues = np.array([ np.dot(cols, C) for C in MagObservations ])

fprs_mag, fprs_other, fprs_comp, tprs_mag, tprs_other, tprs_comp  = [], [], [], [], [], []
aucs_mag, aucs_other, aucs_comp = [], [], []


def MakeCompositeModelFromScratch(xtrain, ytrain, CLFR=lambda : RandomForestClassifier(**rfc_params)):
	

	clfr = CLFR()
	clfr.fit(xtrain, ytrain)
	return None, clfr

def MakeCompositeModelFromProbs(xtrain1, xtrain2, ytrain, CLFR=LDA):
	X = np.array(zip(xtrain1, xtrain2))
	clfr = CLFR()
	scaler = StandardScaler().fit(X)
	X = scaler.transform(X)
	clfr.fit(X, ytrain)
	return scaler, clfr

def EvalModel(Xtrain, Xtest, scaler, clfr):
	if scaler is None:
		PPtrain, PPtest = clfr.predict_proba(Xtrain), clfr.predict_proba(Xtest)
	else:
		PPtrain, PPtest = clfr.predict_proba(scaler.transform(Xtrain)), clfr.predict_proba(scaler.transform(Xtest))
	
	TestProbs = np.array([ p[1] for p in PPtest ])
	TrainProbs = np.array([ p[1] for p in PPtrain ])
	
	return TrainProbs, TestProbs
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
fold_number = 0
for train_index, test_index in StratifiedKFold(MagLabels, n_folds=nfolds,shuffle=True):
	#print len(train_index), len([ l for l in MagLabels[train_index] if l == 1 ])
	fold_number += 1
	print "Training fold %d; %d training RR Lyrae"%(fold_number, len([ l for l in MagLabels[train_index] if l == 1 ]))
	X_mag_train,   X_mag_test 		= MagObservations[train_index], MagObservations[test_index]
	X_other_train, X_other_test 	= OtherObservations[train_index], OtherObservations[test_index]

	Y_train, Y_test 				= MagLabels[train_index], MagLabels[test_index]
	# We assume all labels are the same!! 

	MagScaler, MagModel = MakeMagnitudeModel(X_mag_train, Y_train)
	OtherScaler, OtherModel = MakeOtherModel(X_other_train, Y_train)
	
	OtherTrainProbs, OtherTestProbs = EvalModel(X_other_train, X_other_test, OtherScaler, OtherModel)
	MagTrainProbs, MagTestProbs = EvalModel(X_mag_train, X_mag_test, MagScaler, MagModel)

	#CompositeScaler, CompositeModel = MakeCompositeModelFromProbs(MagTrainProbs, OtherTrainProbs, Y_train)
	
	#X_comp_train = ZipAndMerge(MagTrainProbs, X_other_train)
	#X_comp_test = ZipAndMerge(MagTestProbs, X_other_test)
	X_comp_train = ZipAndMerge(X_mag_train, X_other_train)
	X_comp_test = ZipAndMerge(X_mag_test, X_other_test)

	CompositeScaler, CompositeModel = MakeCompositeModelFromScratch(X_comp_train, Y_train)
	CompositeTrainProbs, CompositeTestProbs = EvalModel(X_comp_train,X_comp_test,CompositeScaler, CompositeModel)

	plot_both=False
	if plot_both:
		f = plt.figure()

		MP_rr  = [ MP for l, MP in zip(Y_test, MagTestProbs)   if l == 1 ]
		MP_nrr = [ MP for l, MP in zip(Y_test, MagTestProbs)   if l == 0 ]

		OP_rr  = [ OP for l, OP in zip(Y_test, OtherTestProbs) if l == 1 ]
		OP_nrr = [ OP for l, OP in zip(Y_test, OtherTestProbs) if l == 0 ]

		ax = f.add_subplot(111)
		ax.scatter(MP_rr, OP_rr, color='r', alpha=0.5)
		ax.scatter(MP_nrr, OP_nrr, color='k', alpha=0.1)
		ax.set_ylabel("LC Prob.")
		ax.set_xlabel("Mag Prob.")
		ax.set_xlim(0,1)
		ax.set_ylim(0,1)
		plt.show()
	if fold_number == 1:
		#PlotMagnitudeModelResults( MagProbs, Y_mag_test)
		PlotRandomForestImportances( OtherModel, OtherKeylist )

	tpr, fpr, _ = roc_curve(Y_test, 1. - MagTestProbs)
	tprs_mag.append(tpr)
	fprs_mag.append(fpr)
	aucs_mag.append(auc( fpr, tpr ))

	tpr, fpr, _ = roc_curve(Y_test,1. - OtherTestProbs)
	tprs_other.append(tpr)
	fprs_other.append(fpr)
	aucs_other.append(auc( fpr, tpr ))

	tpr, fpr, _ = roc_curve(Y_test,1. - CompositeTestProbs)
	tprs_comp.append(tpr)
	fprs_comp.append(fpr)
	aucs_comp.append(auc( fpr, tpr ))

print "AUC (Mag)   = %.3f +/- %.3f"%(np.mean(aucs_mag), np.std(aucs_mag))
print "AUC (Other) = %.3f +/- %.3f"%(np.mean(aucs_other), np.std(aucs_other))
print "AUC (Both)  = %.3f +/- %.3f"%(np.mean(aucs_comp), np.std(aucs_comp))
f = plt.figure(figsize=(12,4),tight_layout=True )
ax_m = f.add_subplot(131)
ax_o = f.add_subplot(132)
ax_b = f.add_subplot(133)


AddROCtoAxis(ax_m, tprs_mag, fprs_mag)
AddROCtoAxis(ax_o, tprs_other, fprs_other, tpr_range=[0.9, 1.])
AddROCtoAxis(ax_b, tprs_comp, fprs_comp, color='b', tpr_range=[0.9, 1])
ax_m.set_title("Magnitude data")
ax_o.set_title("Other data")

#f.subplots_adjust(bottom=0.18)
plt.show()








