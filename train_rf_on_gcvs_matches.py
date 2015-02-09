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
nfolds = 6
#vartypes_to_classify = [ 'M', 'LB', 'EA', 'SR', 'EW', 
#			'RRAB', 'BY', 'UV', 'EA/SD', 'SRB' , 'EB', 
#			'RS', 'L', 'SRA', 'DSCTC', 'E', 'RRC', 
#			'ACV', 'DSCT', 'SRD' , 'INT', 'RR', 'INS',
#			'CST', 'DCEP', 'BE', 'IB', 'EA/DM', 'none']
#vartypes_to_classify = [ 'RRAB', 'none', 'RRC', 'EA', 'EW', 'EB' ]
#vartypes_to_classify = [ 'Pulsating', 'Eclipsing binary', 'none' ]
vartypes_to_classify = [ 'RR Lyrae','none' ]
classify_categories = True
show_confusion_matrix= False
#vartypes_to_classify = [ 'LB', 'EA' ]
#skip_features = [ 'vmag' ]
skip_features = [ 'ra', 'dec', 'V', 'chi2_pf', 'p1_chi2_pf', 'p2_chi2_pf', 'resid_chi2_pf', 'std', 
				'median_absolute_deviation','p1_frac_of_rms_explained', 'p2_frac_of_rms_explained',  
				'raw_lsp_peak1_power', 'raw_lsp_peak2_power', 'raw_lsp_peak3_power', 'raw_lsp_peak4_power', 'raw_lsp_peak5_power',
				'resid_lsp_peak1_power', 'resid_lsp_peak2_power', 'resid_lsp_peak3_power', 'resid_lsp_peak4_power', 'resid_lsp_peak5_power',
				'raw_lsp_peak1_stetson_strlen', 'raw_lsp_peak2_stetson_strlen', 'raw_lsp_peak3_stetson_strlen', 'raw_lsp_peak4_stetson_strlen',
				'raw_lsp_peak5_stetson_strlen', 
				'raw_lsp_peak1_dworetsky_strlen', 'raw_lsp_peak2_dworetsky_strlen', 'raw_lsp_peak3_dworetsky_strlen', 'raw_lsp_peak4_dworetsky_strlen',
				'raw_lsp_peak5_dworetsky_strlen',
				'resid_lsp_peak1_stetson_strlen', 'resid_lsp_peak2_stetson_strlen', 'resid_lsp_peak3_stetson_strlen', 'resid_lsp_peak4_stetson_strlen',
				'resid_lsp_peak5_stetson_strlen', 
				'resid_lsp_peak1_dworetsky_strlen', 'resid_lsp_peak2_dworetsky_strlen', 'resid_lsp_peak3_dworetsky_strlen', 'resid_lsp_peak4_dworetsky_strlen',
				'resid_lsp_peak5_dworetsky_strlen',
				'p1_lsp_peak1_power', 'p1_lsp_peak2_power', 'p1_lsp_peak3_power', 'p1_lsp_peak4_power', 'p1_lsp_peak5_power',
				'p1_lsp_peak1_dworetsky_strlen', 'p1_lsp_peak2_dworetsky_strlen', 'p1_lsp_peak3_dworetsky_strlen', 'p1_lsp_peak4_dworetsky_strlen', 
				'p1_lsp_peak5_dworetsky_strlen', 
				'p1_lsp_peak1_stetson_strlen', 'p1_lsp_peak2_stetson_strlen', 'p1_lsp_peak3_stetson_strlen', 'p1_lsp_peak4_stetson_strlen', 
				'p1_lsp_peak5_stetson_strlen', 
				'p2_lsp_peak1_power', 'p2_lsp_peak2_power', 'p2_lsp_peak3_power', 'p2_lsp_peak4_power', 'p2_lsp_peak5_power',
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



# Load catalog


full_match_cat_name = "Full_GCVS_Cross-matched_HATIDS_maxdist1.00.catalog"
all_match_cat = np.loadtxt(full_match_cat_name, dtype=match_cat_dt)
instances = {}
#for 

# Get a random subset if num is not None
if num is not None:
	all_match_cat = np.loadtxt(full_match_cat_name, dtype=match_cat_dt)
	inds = np.arange(0, len(all_match_cat))
	np.random.shuffle(inds)
	num_per_cat = {}
	for cat in vartypes_to_classify:
		num_per_cat[cat] = 0
	full = lambda : all([ num_per_cat[cat] >= num for cat in vartypes_to_classify ])
	selected_inds = []
	i=0
	for i in inds:
		if full(): break

		m = all_match_cat[i]
		if m['vartype'] in vartypes_to_classify \
						and num_per_cat[m['vartype']] < num\
						and os.path.exists(hat_fname(m['id'])): 

			num_per_cat[m['vartype']] += 1
		else:
			num_per_cat['none'] += 1
			selected_inds.append(i)

	full_match_cat = all_match_cat[selected_inds]
else:
	full_match_cat = np.loadtxt(full_match_cat_name, dtype=match_cat_dt)


categories = {}
Features = {}

hatids_to_fetch = []
for f in full_match_cat:
	hatid = f['id']
	categories[hatid] = 'none'
	if not os.path.exists(hat_fname(hatid)): hatids_to_fetch.append(hatid)
	if classify_categories:
		for v in vartypes_to_classify:
			if v == 'none' : continue
			if f['vartype'] in variable_star_classes[v]:
				categories[hatid] = v
	else:	
		if f['vartype'] in vartypes_to_classify:
			categories[hatid] = f['vartype']

hatids = []
for hatid in full_match_cat['id']:
	if not hatid in phs13_list: continue
	#print categories[hatid]
	if os.path.exists(hat_fname(hatid)) \
			and phs13_list[hatid]['ndet'] > min_ndets \
			and categories[hatid] in vartypes_to_classify\
			and not ( hatid in bad_ids or hatid in look_at_again): \
				hatids.append(hatid)

# Print number of each category
numbers = {}
for var in vartypes_to_classify: numbers[var] = 0
for hid in categories: numbers[categories[hid]]+=1
for var in numbers: print var, numbers[var]


print "TOTAL HATID's SELECTED: %d"%(len(hatids))

def hat_features_fname(hatid):
	return "saved_features/%s-features-v2.pkl"%(hatid)

print "Getting features..."
for nhid, hatid in enumerate(hatids):
	assert(os.path.exists(hat_fname(hatid)))
	feat_fname = hat_features_fname(hatid)
	#print "%s (%d dets) (%d of %d)"%(hatid, phs13_list[hatid]['ndet'], nhid + 1, len(hatids))
	if os.path.exists(feat_fname) and not overwrite:
		#print "  loading saved features.."
		try:
			Features[hatid] = pickle.load(open(feat_fname, 'rb'))
			continue
		except: 
			print "   error loading .. "
	
	#print "  Reading Lightcurve .. "
	LC = rhlc.read_hatlc(hat_fname(hatid))
	#print "  Getting features .. "
	Features[hatid] = fs.get_features(LC)
	#print "  Saving features .. "
	pickle.dump(Features[hatid], open(feat_fname, 'w'))
	#if Features[hatid] is not None:
	#	print Features[hatid]['std']

print "Cleaning features.."
Features_cleaned = {}
for i,hatid in enumerate(Features):
	if Features[hatid] is None: continue
	Features_cleaned[hatid] ={}
	for F in Features[hatid]:
		if F not in skip_features:
			Features_cleaned[hatid][F] = Features[hatid][F]

	if i==0: fdict = Features_cleaned[hatid]
# Convert to observations.

mag_observations = []
other_observations = []
labels = []

mag_key_list = [ key for key in fdict if key in mag_features ]
other_key_list = [ key for key in fdict if key not in mag_features ]
integer_labels = {}
print np.unique([ categories[hid] for hid in Features_cleaned ])
for i, cat in enumerate(np.unique( [ categories[hid] for hid in Features_cleaned ] )):
	integer_labels[cat] = i
label_names = {}
for l in integer_labels:
	label_names[integer_labels[l]] = l

for hid in Features_cleaned:
	fdict = Features_cleaned[hid]
	mag_fvec = []
	other_fvec = []
	#print hid
	for key in mag_features:
		mag_fvec.append(fdict[key])

	for key in fdict:
		if key not in mag_features:
			other_fvec.append(fdict[key])
	#print mag_fvec
	mag_observations.append(mag_fvec)
	other_observations.append(other_fvec)
	labels.append(categories[hid])
#print mag_observations[0], " before arrayed"
mag_observations = np.array(mag_observations)
other_observations = np.array(other_observations)
labels = np.array(labels)
#labels = np.array([ integer_labels[l] for l in labels ])

print mag_observations[0]
AUCs = {}
fprs = {}
tprs = {}
ACCs = []
for l in integer_labels:
	AUCs[l] = { 'mag' : [], 'lc' : [], 'product' : [], 'final' : []}
	
	fprs[l] = { 'mag' : [], 'lc' : [], 'product' : [], 'final' : []}
	tprs[l] = { 'mag' : [], 'lc' : [], 'product' : [], 'final' : []}
	
fno = 1

def handle_magnitudes(xtrain, xtest, ytrain, ytest, CLFR=QDA ):
	#Scale the data to z-scores.
	scaler = StandardScaler().fit(xtrain)
	xtrain_scaled = scaler.transform(xtrain)
	xtest_scaled = scaler.transform(xtest)
	clfr = CLFR()
	clfr.fit(xtrain_scaled, ytrain)

	probs =[ p[0] for p in clfr.predict_proba(xtest_scaled) ]

	probs_per_class = {}
	for c in label_names: 
		probs_per_class[label_names[c]] = [] 

	for p,y in zip(probs, ytest):
		probs_per_class[label_names[y]].append(p)


	fpr,tpr,_ = roc_curve( ytest,  1 - np.array(probs) )
	roc_auc = auc(fpr, tpr)
	print "Magnitude classifier, roc_auc = %f"%(roc_auc)
	f = plt.figure()
	ax = f.add_subplot(111)
	ax.set_title("classifier using color values")
	colors = [ 'r', 'k', 'b', 'g', 'c' ]
	for i, vt in enumerate(vartypes_to_classify):
		ax.hist(probs_per_class[vt], color=colors[i], alpha=0.3, normed=True, label=vt, range=(0,1),bins=20)	
	
	ax.set_ylabel("pdf")
	ax.set_xlabel("Prob score")
	
	if CLFR == LDA:
		all_x, all_y = [], []
		all_x.extend(xtrain)
		all_x.extend(xtest)

		all_y.extend(ytrain)
		all_y.extend(ytest)

		all_xy = zip(all_x, all_y)
		np.random.shuffle(all_xy)
		all_x = [ X for X, Y in all_xy ]
		all_y = [ Y for X, Y in all_xy ]

		scaler = StandardScaler().fit(all_x)
		all_x = scaler.transform(all_x)
		stds = [  np.std([ a[i] for a in all_x ]) for i in range(len(all_x[0]))  ]

		lda = CLFR()
		lda.fit(all_x, all_y)

		W = [ l[0] for l in  lda.scalings_]
		#xbar = lda.xbar_

		# Q = ((x - <x>)/<x^2> - xb)*w
		#   = sum_i x_i*w_i/<x_i^2> - xb_i*w_i - w_i<x_i>/<x_i^2>
		#   = sum_i x_i*w_i/<x_i^2> + Const.
		const = np.divide(W, stds)
		C = np.dot(const, const)
		const/=sqrt(C)
		for c, C in zip(const, mag_key_list):
			print "%.3f(%s)"%(c, C)
	return probs



def plot_importances( RFC ):
	imps = RFC.feature_importances_
	inds = np.argsort(imps)
	imps = np.sort(imps)
	labels = [ other_key_list[i] for i in inds]

	x = np.arange(len(imps)) + 1
	f = plt.figure(figsize=(10,9))
	ax = f.add_subplot(111)

	ax.barh(x - 0.5,imps, height=1.0, alpha=0.5 )
	ax.set_yticks(x)
	ax.set_ylim(min(x) - 0.5, max(x) + 0.5)
	ax.set_yticklabels(labels, fontsize=10)
	f.subplots_adjust(left=0.3, top=0.95, bottom=0.05)
	plt.show()

print len([ l for l in labels if l == 'RR Lyrae' ]), len([ l for l in labels if l == 'none' ])
for train_index, test_index in StratifiedKFold(labels, n_folds=nfolds,shuffle=True):
	print len([ l for l in labels[train_index] if l == 'RR Lyrae' ])
	print "Fold %d..."%(fno)
	#mag_classifier = OneVsRestClassifier(SVC(**svm_params))
	lc_classifier = RandomForestClassifier(**rfc_params)
	#lc_classifier = SVC(**svm_params)

	X_mag_train, X_mag_test = mag_observations[train_index], mag_observations[test_index]
	X_other_train, X_other_test = other_observations[train_index], other_observations[test_index]
	train_labels, test_labels = labels[train_index], labels[test_index]

	Y_train = [ integer_labels[l] for l in train_labels ]
	Y_test = [ integer_labels[l] for l in test_labels]
	

	# Scale the magnitude features
	print "  scaling"
	#mag_scaler = StandardScaler().fit(X_mag_train)
	#X_mag_train = mag_scaler.transform(X_mag_train)
	#X_mag_test = mag_scaler.transform(X_mag_test)
	# Train the magnitude classifier
	#mag_classifier.fit(X_mag_train, Y_train)
	#print X_mag_train[0], X_mag_train[1]
	
	#print LDA_components

	magprobs = handle_magnitudes(X_mag_train, X_mag_test, Y_train, Y_test)
	

	# Train the lc classifier
	print "  training LC classifier"
	lc_classifier.fit(X_other_train, Y_train)

	plot_importances(lc_classifier)
	# Get probability predictions
	#Y_prob_mag = mag_classifier.predict_proba(X_mag_test)
	print "  predicting"
	Y_prob_lc = lc_classifier.predict_proba(X_other_test)
	Y_pred = lc_classifier.predict(X_other_test)
	#print Y_prob_lc

	#print Y_prob_mag
	#print Y_prob_lc
	#print mag_classifier.estimators_, mag_classifier.multilabel_
	#print mag_classifier.classes_
	#sys.exit()
	Y_prob_final = Y_prob_lc
	#Y_prob_final = np.multiply(Y_prob_mag, Y_prob_lc)

	Y_prob_final = np.array(Y_prob_final)
	Y_test  = np.array(Y_test)

	probs_rr = [ p[0] for p, l in zip(Y_prob_final, Y_test) if l == 0 ]
	probs_nonrr = [ p[0] for p, l in zip(Y_prob_final, Y_test) if l == 1 ]

	magprobs_rr = [ p for p, l in zip(magprobs, Y_test) if l == 0]
	magprobs_nonrr = [ p for p, l in zip(magprobs, Y_test) if l == 1]
	F = plt.figure()
	AX = F.add_subplot(111)
	#AX.hist(probs_rr, color='r', alpha=0.3, normed=True, range=(0,1), bins=20)
	#AX.hist(probs_nonrr, color='k', alpha=0.3, normed=True, range=(0,1),bins=20)
	AX.scatter(magprobs_nonrr, probs_nonrr, color='k', alpha=0.1)
	AX.scatter(magprobs_rr, probs_rr, color='r', alpha=0.5)
	AX.set_xlim(0,1)
	AX.set_ylim(0,1)
	plt.show()

	for l in integer_labels:
		#print len(Y_test[:,integer_labels[l]]), len(Y_prob_final[:,integer_labels[l]])
		test = np.zeros(len(Y_prob_final))
		for n, yt in enumerate(Y_test):
			if yt == integer_labels[l]:
				test[n] = 1.
		pred = [ Y_prob_final[n][integer_labels[l]] for n in range(len(Y_prob_final))]
		
			
		
		#print len(test), len(pred)
		fpr,tpr,_ = roc_curve( test,  pred )
		
		#print len(fpr), len(tpr)
		roc_auc = auc(fpr, tpr)

		
		#acc = 0
		#print roc_auc
		#sys.exit()
		fprs[l]['final'].append(fpr)
		tprs[l]['final'].append(tpr)
		AUCs[l]['final'].append(roc_auc)
		
	fno += 1
	def conv(prob_list):
		best_inds = [ i for i in range(len(prob_list)) if prob_list[i] == max(prob_list)]
		#print best_inds
		return min(best_inds)
	pred_float = [ conv(prd) for prd in Y_prob_final ]
	acc = accuracy_score(Y_test, pred_float )
	ACCs.append(acc)
	#print Y_prob_final, Y_pred
	if show_confusion_matrix:
		cm0 = confusion_matrix(Y_test, Y_pred)
		cm = np.zeros(shape=cm0.shape)
		for i in range(len(cm0)):
			cm[i,:] = np.array([ float(cm0[i][j]) for j in range(len(cm0)) ])/float(sum(cm0[i,:]))


		#print(cm)

		# Show confusion matrix in a separate window
		f = plt.figure()
		ax = f.add_subplot(111)
		ms = ax.matshow(cm)
		plt.colorbar(ms)
		ax.set_xticks(.5 + np.arange(len(label_names)))
		#print [ label_names[LN] for LN in range( len(label_names) )  ]
		ax.set_xticklabels( [ label_names[LN] for LN in range( len(label_names) )  ], rotation=45)

		ax.set_yticks(.5 + np.arange(len(label_names)))
		ax.set_yticklabels( [ label_names[LN] for LN in range( len(label_names))  ], rotation=45)
		#ax.set_xlim(1,len(label_names) + 1)
		#ax.set_ylim(1,len(label_names) + 1)
		plt.show()

class_err = 1 - np.mean(ACCs)
class_err_unc = np.std([ 1- a for a in ACCs])
print "Misclassification rate : %.4f +/- %.5f"%(class_err, class_err_unc)

print "%-20s %-6s     %-5s "%("Label", "AUC", "unc.")
print "%-20s %-6s     %-5s "%("-------", "------","-----")

aucs_mean = {}


for l in integer_labels:
	
	#n_thresh = len(fprs[l]['final'][0])
	#n_cv = len(fprs[l]['final'])
	aucs_mean[l] = np.mean(AUCs[l]['final'])
	
	aucs_std = np.std(AUCs[l]['final'])
	
	#fprs_mean[l] = [ np.mean([ fprs[l]['final'][cvno][thresh] for cvno in range(n_cv) ]) for thresh in range(n_thresh) ]
	#tprs_mean[l] = [ np.mean([ tprs[l]['final'][cvno][thresh] for cvno in range(n_cv) ])  for thresh in range(n_thresh) ]
	print "%-20s %-6.3f +/- %-5.3f" %(l, aucs_mean[l], aucs_std)


roc_fig = plt.figure()
ax_roc = roc_fig.add_subplot(111)
ax_roc.set_ylabel('True positive rate')
ax_roc.set_xlabel('False positive rate')
colors = [ 'r', 'g', 'b', 'c', 'y', 'm' ]
ls = [ '-', '--' , ':']
for nl,l in enumerate(integer_labels):
	Color = colors[nl%len(colors)]
	Linestyle = ls[(nl/len(colors)) % len(ls) ]
	for i in range(len(fprs[l]['final'])):
		#print fprs[l]['final'][i], tprs[l]['final'][i]
		#sys.exit()
		if i ==0:

			ax_roc.plot(fprs[l]['final'][i], tprs[l]['final'][i], label="%s (area = %.3f)"%(l, aucs_mean[l]), alpha=0.5, color=Color, ls=Linestyle)
		#else:
			#ax_roc.plot(fprs[l]['final'][i], tprs[l]['final'][i], alpha=0.5,color=Color, ls=Linestyle)


ax_roc.legend(loc='lower right', fontsize=10)
plt.show()









