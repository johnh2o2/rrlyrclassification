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

alpha = 1.0
beta = 0.001

initial_model_rootname = ""

# Get rank and size of mpi process
comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
rank = comm.Get_rank()
ROOT = (rank == 0)

# TODO: get ALL hatids
hatids = [ hatids_in_field ]

def in_equilibrium(scores):
	kl = [ k for k in scores ]
	if len(scores[kl[0]]) <= 1: return False

	dscores = [ scores[k][-1] - scores[k][-2] for k in scores ]
	if np.mean(dscores) < beta: return True
	else: return False

btprs_comp, bfprs_comp = [], []

def make_model(xtrain, ytrain, CLFR=lambda : RandomForestClassifier(**rfc_params)):
	clfr = CLFR()
	clfr.fit(xtrain, ytrain)
	return None, clfr

def process(feats):
	logprint("  Cleaning features...")
	feats = CleanFeatures(feats)

	logprint("  Adding custom features...")
	feats = AddCustomFeatures(feats)

	hatids, keylist, observations, labels = MakeObservations(feats, categories)

	return feats, hatids, keylist, observations, labels


# Load all features
if ROOT:
	if nprocs > 1:
		
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




for Iter in range(num_iterations):
	if rank > num_bags: break # right now we're operating on a 1 bag/core system.
	if ROOT:

		# An extra bag for testing (and possibly the SVM decision maker)
		if fit_model_weights:
			bags = get_bagged_samples(categories, num_bags + 2)
		else:
			bags = get_bagged_samples(categories, num_bags + 1)

		for i,b in enumerate(bags):
			ncats = len(np.unique([ categories[ID] for ID in b ]))
			if ncats < 2: raise Exception("bag number %d only has one class!"%(i))

			# First and last bag are reserved for:
			#    first bag :  The feature bag that Root uses
			#    last bag  :  The feature bag that's used for testing/cross-validation
			if i == 0 or i == len(bags) - 1: continue
			FFs = { ID : All_Features[ID] for ID in b }
			comm.send(obj=FFs, dest=i, tag=0)

		# Features to test this particular bagged model.
		Testing_Features = { ID : All_Features[ID] for ID in bags[-1] }

		# (This is the list of full features for the ROOT node)
		Full_Features = { ID : All_Features[ID] for ID in bags[0] }

		# Get cross-validation features!
		logprint("Generating cross validation data", all_nodes=True)
		test = {}
		test['features'], test['hatids'], test['keylist'], test['observations'], test['labels'] = process(Testing_Features)
		
		if fit_model_weights:
			# Features to use for testing the SVM decision-maker
			SVM_Features = { ID : All_Features[ID] for ID in bags[1] }
			# Get the SVM decision-maker training features
			logprint("Generating SVM training data")
			svm = {}
			svm['features'], svm['hatids'], svm['keylist'], svm['observations'], svm['labels'] = process(SVM_Features)
			
		
	else:
		status = MPI.Status()
		Full_Features = comm.recv(source=ROOT, tag=MPI.ANY_TAG, status=status)

	# Process features
	logprint("Processing features...")
	rfc = {}
	rfc['features'], rfc['hatids'], rfc['keylist'], rfc['observations'], rfc['labels'] = process(Full_Features)
	
	fprs, tprs, aucs = [], [], []
	all_pos, all_neg = [], []
	false_negatives, false_positives = [], []
	
	fold_number = 0
	logprint("TRAINING:")
	for train_index, test_index in StratifiedKFold(rfc['labels'], n_folds=nfolds,shuffle=True):

		fold_number += 1
		logprint("  training fold %d; %d training RR Lyrae"%(fold_number, len([ l for l in MagLabels[train_index] if l == 1 ])))

		X_train,   X_test 				= rfc['observations'][train_index], rfc['observations'][test_index]
		IDs_train, IDs_test 			= rfc['hatids'][train_index], rfc['hatids'][test_index]
		Y_train, Y_test 				= rfc['labels'][train_index], rfc['labels'][test_index]
		
		Scaler, Model = MakeCompositeModelFromScratch(X_train, Y_train)
		TrainProbs, TestProbs = EvalModel(X_train,X_test,Scaler, Model)

		for i, ID, Prob in zip(np.arange(len(IDs_test)),IDs_test, TestProbs):
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

		
		tpr, fpr, _ = roc_curve(Y_test, 1. - TestProbs)
		tprs.append(tpr)
		fprs.append(fpr)
		aucs.append(auc( fpr, tpr ))

		
	logprint( "Model %d of %d \n\t\
AUC (Both)  = %.3f +/- %.3f\n\t\
%d False negatives (RR Lyrae with P(RRLyr) < %.2f)\n\t\
%d False positives (non RR with P(RR) > %.2f)\
	"%(rank+1,nprocs,  
		np.mean(aucs), np.std(aucs),
		len(false_negatives), cutoff, 
		len(false_positives), cutoff
		), 
	all_nodes=True )



	if ROOT:

		nmodels = 1

		# Now train on everything (except for the testing data)
		X_train = rfc['observations']
		X_test = test['observations']

		root_scaler, root_model = MakeCompositeModelFromScratch(X_train, rfc['labels'])

		models = [ root_comp_model ]
		scalers = [ root_comp_scaler ]

		status = MPI.Status()

		# Get other models from the slaves
		while nmodels < num_bags:

			s, m = comm.recv(obj=None, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
			nmodels += 1
			models.append(m)
			scalers.append(s)
			

		# Put all of these models into bags
		BaggedModel = BaggedModel(scalers, models)

		# You can fit an SVM decision maker if you want; otherwise it'll select the maximum score.
		if fit_model_weights: 
			X_svm_train = svm['observations']
			Y_svm_train = svm['labels'] 

			BaggedModel.fit(X_svm_train,Y_svm_train)

		# Root file names of each model
		bagged_model_rootname = "%s/bagged_comp_clfr_iter%04d"%(model_output_dir, iteration)

		# Save these models!!
		# TODO: re-train on ALL labeled data (these are only trained on MOST of the labeled data for
		#       cross-validation purposes.)
		logprint("Saving bagged model.")
		BaggedModel.save(bagged_model_rootname)


		# Now cross validate!
		logprint("Cross validating!!", all_nodes=True)

		TestProbs = BaggedModel.predict_proba(X_test)
		TestProbs = np.array([ p[1] for p in TestProbs ])

		logprint("evaluating cross validation results!", all_nodes=True)
		Y_test = test['labels']

		tpr, fpr, _ = roc_curve(Y_test, 1. - CompTestProbs)

		btprs.append(tpr)
		bfprs.append(fpr)
		
	else:
		X_train = rfc['observations']

		comm.send(obj= MakeCompositeModelFromScratch(rfc['observations'], rfc['labels']), dest=ROOT)



if use_matplotlib and ROOT:
	f = plt.figure(figsize=(4,4),tight_layout=True )
	ax = f.add_subplot(111)

	AddROCtoAxis(ax, btprs, bfprs, tpr_range=[0.99, 1], fpr_range=[0.0, 0.4])
	ax.set_title("Model iteration %d"%(iteration))
	f.suptitle("Bagged model (%d bags)"%(num_bags))
	plt.show(block=True)

	PlotRandomForestImportances( BaggedModel.models,[  rfc['keylist'] for i in range(len(BaggedModel.models)) ] )


if ROOT:
	aucs = [ auc(fpr, tpr) for fpr, tpr in zip(bfprs_comp, btprs_comp) ]

	print "BAGGED MODEL AUC:  %.3f +/- %.5f"%(np.mean(aucs), np.std(aucs))
	
	pickle.dump(("%s"%(bagged_model_rootname), skip_features, vartypes_to_classify, rfc['keylist']), 
			open(get_classifier_fname(iteration), 'wb'))

