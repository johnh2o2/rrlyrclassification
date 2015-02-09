from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import StratifiedKFold
from sklearn.decomposition import RandomizedPCA
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
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

# This is just a copy of train_rf... for now.

match_cat_dt = np.dtype([
		('id', 'S15'),
		('vartype', 'S10'),
		('dist', np.float_)
	])
num = None
min_ndets = 20
nfolds = 3
#vartypes_to_classify = [ 'M', 'LB', 'EA', 'SR', 'EW', 
#			'RRAB', 'BY', 'UV', 'EA/SD', 'SRB' , 'EB', 
#			'RS', 'L', 'SRA', 'DSCTC', 'E', 'RRC', 
#			'ACV', 'DSCT', 'SRD' , 'INT', 'RR', 'INS',
#			'CST', 'DCEP', 'BE', 'IB', 'EA/DM', 'none']
#vartypes_to_classify = [ 'RRAB', 'none', 'RRC', 'EA', 'EW', 'EB' ]
#vartypes_to_classify = [ 'Pulsating', 'Eclipsing binary', 'none' ]
vartypes_to_classify = [ 'RRAB', 'none' ]
classify_categories = False
show_confusion_matrix= False
#vartypes_to_classify = [ 'LB', 'EA' ]
#skip_features = [ 'vmag' ]
skip_features = [ ]
#mag_features = [ 'R-V', 'I-V', 'J-V', 'H-V', 'K-V' ]
#mag_features = [ 'vmag' ]
mag_features = [ 'vmag' ]
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

variable_star_classes = {
	'Eruptive' : ['FU', 'GCAS', 'I', 'IA', 'IB', 
			'IN', 'INA', 'INB', 'INT', 'IT', 'IN(YY)', 
			'IS', 'ISA', 'ISB', 'RCB', 'RS', 'SDOR', 'UV', 'UVN', 'WR'],

	'Pulsating' : ['ACYG', 'BCEP', 'BCEPS', 'CEP', 'CEP(B)', 'CW', 'CWA', 'CWB', 'DCEP', 'DCEPS',
           'DSCT', 'DSCTC', 'GDOR', 'L', 'LB', 'LC', 'M', 'PVTEL', 'RPHS', 'RR', 'RR(B)', 'RRAB',
           'RRC', 'RV', 'RVA', 'RVB', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'SXPHE', 'ZZ', 'ZZA', 'ZZB'],

    'Rotating' : [ 'ACV', 'ACVO', 'BY', 'ELL', 'FKCOM', 'PSR', 'SXARI' ],

    'Cataclysmic' : ['N', 'NA', 'NB', 'NC', 'NL', 'NR',
           'SN', 'SNI', 'SNII', 'UG', 'UGSS', 'UGSU', 'UGZ', 'ZAND'],

    'Eclipsing binary' : ['E', 'EA', 'EB', 'EW', 'GS', 'PN', 'RS', 'WD', 'WR', 'AR', 'D', 'DM',
           'DS', 'DW', 'K', 'KE', 'KW', 'SD']
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
			
bad_ids = [ 'HAT-079-0000101', 'HAT-128-0000156', 'HAT-141-0001285', 'HAT-141-0004548', 'HAT-142-0004019'
'HAT-150-0012878', 'HAT-168-0002894', 'HAT-189-0002780', 'HAT-196-0018339', 'HAT-207-0011053', 
'HAT-248-0000036', 'HAT-277-0004093', 'HAT-287-0017860', 'HAT-292-0028865', 'HAT-339-0136924',
'HAT-362-0002588', 'HAT-388-0000557' ]
look_at_again = [ 'HAT-223-0003186', 'HAT-237-0002943', 'HAT-242-0026174', 'HAT-242-0034689','HAT-256-0005695',
'HAT-292-0100671', 'HAT-332-0001158', 'HAT-339-0101490', 'HAT-363-0012214','HAT-431-0000070', 'HAT-437-0000456' ]


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
	return "saved_features/%s-features.pkl"%(hatid)

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

	for key in mag_features:
		mag_fvec.append(fdict[key])
	for key in fdict:
		if key in mag_features:
			mag_fvec.append(fdict[key])
		else:
			other_fvec.append(fdict[key])

	mag_observations.append(mag_fvec)
	other_observations.append(other_fvec)
	labels.append(categories[hid])
mag_observations = np.array(mag_observations)
other_observations = np.array(other_observations)
labels = np.array(labels)
#labels = np.array([ integer_labels[l] for l in labels ])


AUCs = {}
fprs = {}
tprs = {}
ACCs = []
for l in integer_labels:
	AUCs[l] = { 'mag' : [], 'lc' : [], 'product' : [], 'final' : []}
	
	fprs[l] = { 'mag' : [], 'lc' : [], 'product' : [], 'final' : []}
	tprs[l] = { 'mag' : [], 'lc' : [], 'product' : [], 'final' : []}
	
fno = 1
for train_index, test_index in StratifiedKFold(labels, n_folds=nfolds,shuffle=True):
	print "Fold %d..."%(fno)
	mag_classifier = OneVsRestClassifier(SVC(**svm_params))
	lc_classifier = RandomForestClassifier(**rfc_params)

	X_mag_train, X_mag_test = mag_observations[train_index], mag_observations[test_index]
	X_other_train, X_other_test = other_observations[train_index], other_observations[test_index]
	train_labels, test_labels = labels[train_index], labels[test_index]

	Y_train = [ integer_labels[l] for l in train_labels ]
	Y_test = [ integer_labels[l] for l in test_labels]
	#for i,l in enumerate(test_labels):
	#	lab_vec = []
	#	for j in range(len(integer_labels)):
	#		lab_vec.append(0)
	#	lab_vec[integer_labels[l]] = 1
	#	Y_test.append(lab_vec)

	#for i,l in enumerate(train_labels):
	#	lab_vec = []
	#	for j in range(len(integer_labels)):
	#		lab_vec.append(0)
	#	lab_vec[integer_labels[l]] = 1
	#	Y_train.append(lab_vec)
		
	#print Y_train
	#Y_train = [ integer_labels[l] for l in train_labels ]
	#Y_test = [ integer_labels[l] for l in test_labels ]


	#print Y_train

	# Scale the magnitude features
	print "  scaling"
	mag_scaler = StandardScaler().fit(X_mag_train)
	X_mag_train = mag_scaler.transform(X_mag_train)
	X_mag_test = mag_scaler.transform(X_mag_test)
	# Train the magnitude classifier
	#mag_classifier.fit(X_mag_train, Y_train)

	# Train the lc classifier
	print "  training LC classifier"
	lc_classifier.fit(X_other_train, Y_train)


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









