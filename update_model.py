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

import numpy as np
import os, sys
from scipy.interpolate import interp1d

import utils.feature_selection as fs
from utils.miscutils import *
from utils.featureutils import *
import utils.readhatlc as rhlc
from settings import *
import cPickle as pickle
from sklearn.metrics import confusion_matrix

# What iteration are we on?
iteration = get_iteration_number() 
logprint("Updating model; iteration %d"%(iteration))

# Get the labeled HATIDs.
logprint("Loading labeled HATIDs")
hatids, categories = LoadLabeledHatIDs()

logprint("Getting features...")
Full_Features 						= LoadAllFeatures(hatids)
#PlotMagnitudeDist(Full_Features)
#obs_mag = {}
#for ID in mean_mags:
#	obs_mag[ID] = { 'V' : mean_mags[ID] }
#PlotMagnitudeDist(obs_mag)
#PlotMagnitudeDist(Full_Features)
#plt.show()
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
logprint("Cleaning features...")
Features 							= CleanFeatures(Full_Features)
logprint("Adding custom features...")
Features                            = AddCustomFeatures(Features)
logprint("Splitting features into magnitudes and (others)...")
MagFeatures, OtherFeatures 			= SplitFeatures(Features, mag_features)

#flist = [ f for f in OtherFeatures[hatids[0]] ]
#print flist, len(flist)
# Convert to observations.

if not MagFeatures is None:
	logprint("Making training samples for the magnitude data")
	MagIDs, MagKeylist, MagObservations, MagLabels = MakeObservations(MagFeatures, categories)
logprint("Making training samples for the other data")
OtherIDs, OtherKeylist, OtherObservations, OtherLabels = MakeObservations(OtherFeatures, categories)
#for o in OtherKeylist: print o
#sys.exit()

cols = GetBestColor(MagObservations, MagLabels)
LDA_ColorValues = np.array([ np.dot(cols, C) for C in MagObservations ])

fprs_mag, fprs_other, fprs_comp, tprs_mag, tprs_other, tprs_comp  = [], [], [], [], [], []
aucs_mag, aucs_other, aucs_comp = [], [], []
all_pos, all_neg = [ ], []
false_negatives, false_positives = [],[]

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

fold_number = 0
logprint("TRAINING:")
for train_index, test_index in StratifiedKFold(MagLabels, n_folds=nfolds,shuffle=True):

	#print len(train_index), len([ l for l in MagLabels[train_index] if l == 1 ])
	fold_number += 1
	logprint("  training fold %d; %d training RR Lyrae"%(fold_number, len([ l for l in MagLabels[train_index] if l == 1 ])))

	X_mag_train,   X_mag_test 		= MagObservations[train_index], MagObservations[test_index]
	X_other_train, X_other_test 	= OtherObservations[train_index], OtherObservations[test_index]
	IDs_train, IDs_test 			= OtherIDs[train_index], OtherIDs[test_index]
	Y_train, Y_test 				= MagLabels[train_index], MagLabels[test_index]
	assert(all(Y_train == OtherLabels[train_index]))
	# We assume all labels are the same!! 

	MagScaler, MagModel = MakeMagnitudeModel(X_mag_train, Y_train)
	MagTrainProbs, MagTestProbs = EvalModel(X_mag_train, X_mag_test, MagScaler, MagModel)
	
	OtherScaler, OtherModel = MakeOtherModel(X_other_train, Y_train)
	OtherTrainProbs, OtherTestProbs = EvalModel(X_other_train, X_other_test, OtherScaler, OtherModel)
	#CompositeScaler, CompositeModel = MakeCompositeModelFromProbs(MagTrainProbs, OtherTrainProbs, Y_train)
	
	#X_comp_train = ZipAndMerge(MagTrainProbs, X_other_train)
	#X_comp_test = ZipAndMerge(MagTestProbs, X_other_test)
	X_comp_train = ZipAndMerge(X_mag_train, X_other_train)
	X_comp_test = ZipAndMerge(X_mag_test, X_other_test)

	CompositeScaler, CompositeModel = MakeCompositeModelFromScratch(X_comp_train, Y_train)
	CompositeTrainProbs, CompositeTestProbs = EvalModel(X_comp_train,X_comp_test,CompositeScaler, CompositeModel)

	for i, ID, Prob in zip(np.arange(len(IDs_test)),IDs_test, CompositeTestProbs):
		assert(categories[ID] == label_names[Y_test[i]])
		if categories[ID] != 'none':
			all_pos.append(Prob)
		else:
			all_neg.append(Prob)

		#if categories['ID']
		if abs(Y_test[i] - Prob) > 1 - cutoff:
			if categories[ID] is 'none':
				false_positives.append((ID, Prob))
			else:
				false_negatives.append((ID, Prob))
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
	#if fold_number == 1:
		
		#PlotRandomForestImportances( OtherModel, OtherKeylist )
		#PlotProbHistogram(CompositeTestProbs, Y_test)
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
AddROCtoAxis(ax_b, tprs_comp, fprs_comp, color='b', tpr_range=[0.99, 1], fpr_range=[0.0, 0.4])
ax_m.set_title("Magnitude data")
ax_o.set_title("Other data")

#for failure in failed_ids:
print len(false_negatives), "False negatives (RR Lyrae with P(RRLyr) < %.2f)"%(cutoff)
print len(false_positives), "False positives (non RR with P(RR) > %.2f)"%(cutoff)

gcvs_crossmatches 					= np.loadtxt(full_gcvs_match_cat_fname , dtype=match_cat_dt)

for ID, prob in false_negatives:
	j = gcvs_crossmatches['id'].tolist().index(ID)
	print ID, prob, gcvs_crossmatches['vartype'][j], categories[ID]
for ID, prob in false_negatives:
	print ID,
print ""

fhist = plt.figure()
axhist = fhist.add_subplot(111)
axhist.hist(all_pos, color='r', alpha=0.3, range=(0,1), bins=50, normed=True)
axhist.hist(all_neg, color='k', alpha=0.3, range=(0,1), bins=50, normed=True)

pickle.dump(MakeOtherModel(MagObservations, MagLabels) + MakeMagnitudeModel(OtherObservations, OtherLabels) + \
	(skip_features, mag_features,vartypes_to_classify, OtherKeylist, MagKeylist), open(get_classifier_fname(iteration + 1), 'wb'))

#f.subplots_adjust(bottom=0.18)
plt.show()
