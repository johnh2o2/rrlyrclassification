from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.lda import LDA
from sklearn.qda import QDA
from sklearn.metrics import roc_curve, auc, accuracy_score, roc_auc_score, log_loss
from sklearn.learning_curve import learning_curve, validation_curve

from settings import *
use_matplotlib = True
if use_matplotlib: 
	if RUNNING_ON_DELLA:
		print "update_model2: setting matplotlib backend to agg"
		import matplotlib as mpl
		mpl.use('Agg')
	import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy.interpolate import interp1d
import utils.feature_selection as fs
from utils.miscutils import *
from utils.featureutils import *
import utils.readhatlc as rhlc
import cPickle as pickle

iteration = get_iteration_number() + 1 # We're making the next model!

# Get the labeled HATIDs.
logprint("Loading labeled HATIDs")
hatids, categories = LoadLabeledHatIDs()

# Process the features
def process(feats):
	logprint("  Cleaning features...")
	feats = CleanFeatures(feats)

	logprint("  Adding custom features...")
	feats = AddCustomFeatures(feats)

	hatids, keylist, observations, labels = MakeObservations(feats, categories)

	return feats, hatids, keylist, observations, labels

# Load all the features
All_Features = LoadAllFeatures(hatids)

# Process the features
Features, Hatids, Keylist, Observations, Labels = process(All_Features)

# Train the model
logprint("TRAINING:")
fold_number = 0
tprs, fprs, aucs, oob_scores = [], [], [], []

for train_index, test_index in StratifiedKFold(Labels, n_folds=nfolds, shuffle=True):

	fold_number 					+= 1
	logprint("  training fold %d; %d training RR Lyrae"%(fold_number, len([ l for l in Labels[train_index] if l == 1 ])))

	X_train, X_test 				= Observations[train_index], Observations[test_index]
	Y_train, Y_test 				= Labels[train_index], Labels[test_index]
	
	Model 							= RandomForestClassifier(**rfc_params).fit(X_train, Y_train)
	TrainProbs, TestProbs 			= EvalModel(X_train, X_test, None, Model)
	
	oob_scores.append(Model.oob_score_)

	tpr, fpr, thresholds = roc_curve(Y_test, 1. - TestProbs)
	tprs.append(tpr)
	fprs.append(fpr)
	aucs.append(auc( fpr, tpr ))

	# Find the right threshold!
	logprint("   %-10s %-20s %-20s"%("threshold", "tpr", "fpr"))
	for i,t in enumerate(thresholds):
		logprint("  %-10.4f %-6.4f      %-6.4f"%(t, tpr[i], fpr[i]))

	logprint("   OOB_SCORE = %.6f"%(Model.oob_score_))

# Now train full model.
Full_Model = RandomForestClassifier(**rfc_params).fit(Observations, Labels)
final_oob_score = Full_Model.oob_score_

# Save model + relevant information
pickle.dump(Full_Model, open(get_classifier_fname(iteration), 'wb'))
pickle.dump((skip_features, vartypes_to_classify, Keylist ), open(get_classifier_information_fname(iteration), 'wb'))

# Print OOB scores
print "OOB  Scores"
print "Avg     = %.5f +/- %.5f"%(np.mean(oob_scores), np.std(oob_scores))
print "Final   = %.5f"%(final_oob_score)

# Do a more careful assessment of the learning curve
scoring_func = log_loss
scoring = lambda model, x, y : scoring_func(y, model.predict_proba(x))
estimator = RandomForestClassifier(**rfc_params)
train_sizes = np.linspace(0.1, 1.0, 10)
train_sizes_abs, train_scores, test_scores = learning_curve(estimator, Observations, Labels, scoring=scoring, train_sizes=train_sizes)

# TODO: add learning_validation to see how this depends upon number of estimators...

# Plot results (if need be.)
if use_matplotlib:

	# Plot learning curve information
	f_lc = plt.figure(figsize=(4,4), tight_layout=True)
	ax_lc = f_lc.add_subplot(111)
	train_scores_avg =[  np.mean(ts) for ts in train_scores ]
	test_scores_avg = [ np.mean(ts) for ts in test_scores] 
	ax_lc.plot(train_sizes_abs, train_scores_avg, label='Training')
	ax_lc.plot(train_sizes_abs, test_scores_avg, label='Testing')
	ax_lc.legend(loc='best')
	ax_lc.set_xlabel("Number of training samples")
	ax_lc.set_ylabel("Log-loss score")
	f_lc.suptitle("RandomForestClassifier (%d estimators)"%(rfc_params['n_estimators']))
	f_lc.savefig("%s/model_iter%d_learning_curve"%(parent_dir, iteration))

	# Plot ROC results.
	f = plt.figure(figsize=(4,4),tight_layout=True )
	ax = f.add_subplot(111)
	AddROCtoAxis(ax, tprs, fprs, tpr_range=[0.99, 1], fpr_range=[0.0, 0.3])
	ax.set_title("Model iteration %d"%(iteration))
	f.suptitle("RandomForestClassifier (%d estimators)"%(rfc_params['n_estimators']))
	f.savefig('%s/model_iter%d_cross_validation.png'%(parent_dir, iteration))

	# Plot which variables are the most important
	PlotRandomForestImportances( Full_Model, Keylist, savefig='%s/model_iter%d_feature_imps.png'%(parent_dir, iteration ))


