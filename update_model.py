from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.lda import LDA
from sklearn.qda import QDA
from sklearn.metrics import roc_curve, auc, accuracy_score

use_matplotlib = True
if use_matplotlib: import matplotlib.pyplot as plt
import masterslave as msl
import numpy as np
import os, sys
from scipy.interpolate import interp1d
import utils.feature_selection as fs
from utils.miscutils import *
from utils.featureutils import *
import utils.readhatlc as rhlc
from settings import *
import cPickle as pickle
from mpi4py import MPI

fit_model_weights = True

# Get rank and size of mpi process
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
ROOT = (rank == 0)

# What iteration are we on?
iteration = get_iteration_number() 
logprint("Updating model; iteration %d"%(iteration))

# Get the labeled HATIDs.
logprint("Loading labeled HATIDs")
hatids, categories = LoadLabeledHatIDs()

# Get the necessary lightcurves
if iteration == 0:
	if not os.path.isdir(data_dir): os.makedirs(data_dir)

	# We're using the online server for this one...
	# slow, but we don't yet have the infrastructure to get an arbitrary hatid...
	get_lc_fname = get_gzipped_csv_lc_fname
	load_lightcurve = lambda hatid : rhlc.read_hatlc(get_lc_fname(hatid))

	# Find missing ID's
	missing_hatids = []
	for ID in hatids:
		if not os.path.exists(get_lc_fname(ID)): missing_hatids.append(ID)

	# Download them.
	if len(missing_hatids) > 0: fetch_lcs(missing_hatids)

# Now distribute the workload to load/make the features...
batch_size = 10
work_packets = [ hatids[i:min([ len(hatids), i+batch_size])].tolist() for i in range(len(hatids)) ]
if len(work_packets[-1]) == 0: work_packets.pop()

# Master/slave workload distribution
logprint("Getting features!")
num_bags = size
num_iterations = 5

btprs_other, btprs_mag, btprs_comp = [], [], []
bfprs_other, bfprs_mag, bfprs_comp = [], [], []


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
if ROOT:
	if size > 1:
		
		all_ft = msl.master(work_packets)
		
		# Combine the feature vectors!
		All_Features = {}
		for Ft in all_ft:
			for ID in Ft:
				All_Features[ID] = Ft[ID]
		
	else:
		All_Features = LoadAllFeatures(hatids)
else:
	msl.slave(LoadAllFeatures)
	
def process(feats):
	logprint("  Cleaning features...")
	feats = CleanFeatures(feats)

	logprint("  Adding custom features...")
	feats = AddCustomFeatures(feats)

	logprint("  Splitting features...")
	magfeats, otherfeats = SplitFeatures(feats, mag_features)

	logprint("  Making observations...")
	if not magfeats is None:
		magids, magkl, magobs, maglabs = MakeObservations(magfeats, categories)
	otherids, otherkl, otherobs, otherlabs = MakeObservations(otherfeats, categories)

	return feats, magfeats, otherfeats, magids, magkl, magobs, maglabs, otherkl, otherobs


for Iter in range(num_iterations):
	if rank > num_bags: break # right now we're operating on a 1 core/bag system.
	if ROOT:

		bags = get_bagged_samples(categories, num_bags + 1)
		for i,b in enumerate(bags):
			ncats = len(np.unique([ categories[ID] for ID in b ]))
			if ncats < 2: raise Exception("bag number %d only has one class!"%(i))
			if i == 0 or i == len(bags) - 1: continue
			FFs = { ID : All_Features[ID] for ID in b }
			comm.send(obj=FFs, dest=i, tag=0)

		Testing_Features = { ID : All_Features[ID] for ID in bags[-1] }
		Full_Features = { ID : All_Features[ID] for ID in bags[0] }
		SVM_Features = { ID : All_Features[ID] for ID in bags[1] }

		# Get cross-validation data!
		logprint("Generating cross validation data", all_nodes=True)
		Testing_Features, Testing_MagFeatures, Testing_OtherFeatures, MagIDs, MagKeylist, \
		Testing_MagObservations, Testing_MagLabels, OtherKeylist, Testing_OtherObservations = process(Testing_Features)

		logprint("Generating SVM training data")
		SVM_Features, SVM_MagFeatures, SVM_OtherFeatures, MagIDs, MagKeylist, \
		SVM_MagObservations, SVM_MagLabels, OtherKeylist, SVM_OtherObservations = process(SVM_Features)


		
	else:
		status = MPI.Status()
		Full_Features = comm.recv(source=ROOT, tag=MPI.ANY_TAG, status=status)

	# Process features
	logprint("Processing features...")
	Features, MagFeatures, OtherFeatures, MagIDs, MagKeylist, \
	MagObservations, MagLabels, OtherKeylist, OtherObservations = process(Full_Features)

	OtherIDs = MagIDs
	OtherLabels = MagLabels

	logprint("  Using LDA to get the best color index")
	# Use an LDA to get the best color index (linear combination of X_i - V)
	cols = GetBestColor(MagObservations, MagLabels)
	LDA_ColorValues = np.array([ np.dot(cols, C) for C in MagObservations ])

	fprs_mag, fprs_other, fprs_comp, tprs_mag, tprs_other, tprs_comp  = [], [], [], [], [], []
	aucs_mag, aucs_other, aucs_comp = [], [], []
	all_pos, all_neg = [], []
	false_negatives, false_positives = [], []



	fold_number = 0
	logprint("TRAINING:")
	for train_index, test_index in StratifiedKFold(MagLabels, n_folds=nfolds,shuffle=True):

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

			if abs(Y_test[i] - Prob) > 1 - cutoff:
				if categories[ID] is 'none':
					false_positives.append((ID, Prob))
				else:
					false_negatives.append((ID, Prob))

		plot_both=False
		if plot_both and use_matplotlib:
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

	logprint( "Model %d of %d \n\t\
AUC (Mag)   = %.3f +/- %.3f\n\t\
AUC (Other) = %.3f +/- %.3f\n\t\
AUC (Both)  = %.3f +/- %.3f\n\t\
%d False negatives (RR Lyrae with P(RRLyr) < %.2f)\n\t\
%d False positives (non RR with P(RR) > %.2f)\
	"%(rank+1,size, 
		np.mean(aucs_mag), np.std(aucs_mag), 
		np.mean(aucs_other), np.std(aucs_other), 
		np.mean(aucs_comp), np.std(aucs_comp),
		len(false_negatives), cutoff, 
		len(false_positives), cutoff
		), 
	all_nodes=True )



	if ROOT:

		nmodels = 1

		X_svm_train_other = SVM_OtherObservations
		X_svm_train_mag = SVM_MagObservations
		X_svm_train_comp = ZipAndMerge(SVM_MagObservations, SVM_OtherObservations)
		Y_svm_train = SVM_MagLabels

		X_comp_train = ZipAndMerge(MagObservations, OtherObservations)
		X_comp_test = ZipAndMerge(Testing_MagObservations, Testing_OtherObservations)

		root_other_scaler, root_other_model = MakeOtherModel(OtherObservations, OtherLabels)
		root_mag_scaler, root_mag_model = MakeMagnitudeModel(MagObservations, MagLabels)
		root_comp_scaler, root_comp_model = MakeCompositeModelFromScratch(X_comp_train, MagLabels)

		other_models = [ root_other_model ]
		other_scalers = [ root_other_scaler ]

		mag_models = [ root_mag_model ]
		mag_scalers = [ root_mag_scaler ]

		comp_models = [ root_comp_model ]
		comp_scalers = [ root_comp_scaler ]

		status = MPI.Status()
		while nmodels < num_bags:

			mags, magm, others, otherm, comps, compm = comm.recv(obj=None, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
			nmodels += 1

			other_models.append(otherm)
			other_scalers.append(others)

			mag_models.append(magm)
			mag_scalers.append(mags)

			comp_models.append(compm)
			comp_scalers.append(comps)

		BaggedMagModel = BaggedModel(mag_scalers, mag_models)
		print np.array(Y_svm_train).shape, np.array(X_svm_train_mag).shape, np.array(X_svm_train_other).shape, np.array(X_svm_train_comp).shape
		if fit_model_weights: BaggedMagModel.fit(X_svm_train_mag, Y_svm_train)

		BaggedOtherModel = BaggedModel(other_scalers, other_models)
		if fit_model_weights: BaggedOtherModel.fit(X_svm_train_other, Y_svm_train)

		BaggedCompModel = BaggedModel(comp_scalers, comp_models)
		if fit_model_weights: BaggedCompModel.fit(X_svm_train_comp, Y_svm_train)

		bagged_mag_model_rootname = "%s/bagged_mag_clfr_iter%04d"
		bagged_other_model_rootname = "%s/bagged_other_clfr_iter%04d"
		bagged_comp_model_rootname = "%s/bagged_comp_clfr_iter%04d"

		JunkScaler = lambda x : x

		logprint("Cross validating!!", all_nodes=True)

		logprint(" (magnitude model)")
		MagTestProbs = BaggedMagModel.predict_proba(Testing_MagObservations)
		MagTestProbs = np.array([ p[1] for p in MagTestProbs ])

		logprint(" (other model)")
		OtherTestProbs = BaggedOtherModel.predict_proba(Testing_OtherObservations)
		OtherTestProbs = np.array([ p[1] for p in OtherTestProbs ])

		logprint(" (composite model)")
		CompTestProbs = BaggedCompModel.predict_proba(X_comp_test)
		CompTestProbs = np.array([ p[1] for p in CompTestProbs ])

		logprint("evaluating cross validation results!", all_nodes=True)
		Y_test = Testing_MagLabels

		tpr_mag, fpr_mag, _ = roc_curve(Y_test, 1. - MagTestProbs)
		tpr_other, fpr_other, _ = roc_curve(Y_test, 1. - OtherTestProbs)
		tpr_comp, fpr_comp, _ = roc_curve(Y_test, 1. - CompTestProbs)

		btprs_comp.append(tpr_comp)
		btprs_other.append(tpr_other)
		btprs_mag.append(tpr_mag)

		bfprs_comp.append(fpr_comp)
		bfprs_other.append(fpr_other)
		bfprs_mag.append(fpr_mag)	

	else:
		X_comp_train = ZipAndMerge(MagObservations, OtherObservations)

		comm.send(obj=  MakeMagnitudeModel(MagObservations, MagLabels) +\
						 MakeOtherModel(OtherObservations, OtherLabels) +\
						 MakeCompositeModelFromScratch(X_comp_train, MagLabels),

						 dest=ROOT)



if use_matplotlib and ROOT:
	f = plt.figure(figsize=(12,4),tight_layout=True )
	ax_m = f.add_subplot(131)
	ax_o = f.add_subplot(132)
	ax_b = f.add_subplot(133)

	AddROCtoAxis(ax_m, btprs_mag, bfprs_mag)
	AddROCtoAxis(ax_o, btprs_other, bfprs_other, tpr_range=[0.9, 1.])
	AddROCtoAxis(ax_b, btprs_comp, bfprs_comp, color='b', tpr_range=[0.99, 1], fpr_range=[0.0, 0.4])
	ax_m.set_title("Magnitude data")
	ax_o.set_title("Other data")
	ax_b.set_title("ALL data")
	f.suptitle("Bagged model (%d bags)"%(num_bags))
	plt.show()


if ROOT:
	aucs_mag = [ auc(fpr, tpr) for fpr, tpr in zip(bfprs_mag, btprs_mag) ]
	aucs_other = [ auc(fpr, tpr) for fpr, tpr in zip(bfprs_other, btprs_other) ]
	aucs_comp = [ auc(fpr, tpr) for fpr, tpr in zip(bfprs_comp, btprs_comp) ]

	
	print "BAGGED MODEL AUC's:\n\t\
Mag:   %.3f +/- %.5f\n\t\
Other: %.3f +/- %.5f\n\t\
Comp:  %.3f +/- %.5f\n\t\
"%(np.mean(aucs_mag), np.std(aucs_mag),
	np.mean(aucs_other), np.std(aucs_other),
	np.mean(aucs_comp), np.std(aucs_comp))
	
	pickle.dump(("%s.pkl"%(bagged_other_model_rootname),) + ("%s.pkl"%(bagged_mag_model_rootname),) + \
			(skip_features, mag_features,vartypes_to_classify, OtherKeylist, MagKeylist), 
			open(get_classifier_fname(iteration + 1), 'wb'))

