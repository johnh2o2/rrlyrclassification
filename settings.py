import re, sys, os
from lcutils.lcutils_config import *

RUNNING_ON_DELLA = True
model_prefix = "rrab_v0"
field_to_analyze = '219'

if RUNNING_ON_DELLA:
	parent_dir = '/home/jah5/rrlyr_search/rrlyrclassification'
	SCRATCH = '/tigress/jah5/rrlyr_scratch'
else:
	parent_dir = '/Users/jah5/Documents/Fall2014_Gaspar'
	SCRATCH = '%s/rrlyr_classification'%(parent_dir)
	DYLD_LIBRARY_PATH = '/opt/local/lib'
data_dir = '%s/data'%(SCRATCH)
model_dir = "%s/models"%(SCRATCH)
model_output_dir = "%s/%s"%(model_dir, model_prefix)
feat_dir = "%s/features"%(SCRATCH)


labeled_hatids_fname = "%s/labeled_hatids.txt"%(model_output_dir)
dt_labeled_hatids = np.dtype([
	('ID', 'S15'), ('label', np.int_), ('iter_detected', np.int_)
	])
keylist_dt = np.dtype([
		('TSTF', 'S8'),
		('EXP', float),
		('BJD', float),
		('FLT', str),
		('unknown1', float),
		('unknown2', float),
		('unknown3', float),
		('unknown4', float),
		('unknown5', float),
		('unknown6', float),
		('unknown7', float),
	])
twomass_dt = np.dtype([ 	('hatid', 'S15'),
						('ra', float),
						('dec', float),
						('jmag', float),
						('jerr', float),
						('hmag', float),
						('herr', float),
						('kmag', float),
						('kerr', float),
						('flags', 'S3'),
						('Bmag', float),
						('Vmag', float),
						('Rmag', float),
						('Imag', float),
						('umag', float),
						('gmag', float),
						('rmag', float),
						('imag', float),
						('zmag', float),
						('hatfield', int),
						('hatfield objid', int)
				])

num = None
min_ndets = 20
nfolds = 10
cutoff = 0.05
overwrite = False

COL_TYPE = 'TF'
COL_SELECTION = 'locally-brightest'

nharmonics = 8
npers      = 1

NPEAKS_TO_SAVE = 5

n_peaks    = 1
DELTA_I    = 10
delta_P    = .01

DPHASE     = 1.0
NSEARCH    = 1

ofac       = 3
hifac      = 2
MACC       = 3
NBINS      = 50
eps        = 10E-7
max_per    = 10.

use_bootstrap = True
boot_size  = 100
n_boots = 1
fraction_of_smallest_stds_to_use = 0.5


time_col = 'BJD'
grpsize = 10
bad_ids = [ 'HAT-079-0000101', 'HAT-128-0000156', 'HAT-141-0001285', 'HAT-141-0004548', 'HAT-142-0004019'
'HAT-150-0012878', 'HAT-168-0002894', 'HAT-189-0002780', 'HAT-196-0018339', 'HAT-207-0011053', 
'HAT-248-0000036', 'HAT-277-0004093', 'HAT-287-0017860', 'HAT-292-0028865', 'HAT-339-0136924',
'HAT-362-0002588', 'HAT-388-0000557', ' HAT-135-0007139 ', 'HAT-189-0006202', 'HAT-239-0006835','HAT-241-0014081',
'HAT-241-0018480']
look_at_again = [ 'HAT-223-0003186', 'HAT-237-0002943', 'HAT-242-0026174', 'HAT-242-0034689','HAT-256-0005695',
'HAT-292-0100671', 'HAT-332-0001158', 'HAT-339-0101490', 'HAT-363-0012214','HAT-431-0000070', 'HAT-437-0000456', 
'HAT-142-0004019', 'HAT-384-0061789', 'HAT-199-1738692', 'HAT-341-0078624', 'HAT-230-0003941', 'HAT-190-0005199',
'HAT-167-0025662']
#HAT-384-0061789 -- wrong period (v3); possibly improved by doing min(stetson). Also a dip prior to the peak?
#HAT-199-1738692 -- RR Lyr but noisy
#HAT-341-0078624 -- ^ same
#HAT-230-0003941 -- ^ same
#HAT-190-0005199 -- not sure what to make of this...
#HAT-167-0025662 -- noisy
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
           'DS', 'DW', 'K', 'KE', 'KW', 'SD'],
    'RR Lyrae' : [ 'RRAB', 'RRC', 'R', 'RR', 'RR(B)' ]
}

hat_features_fname = lambda hatid, model_name : "%s/%s-features-%s.pkl"%(feat_dir, hatid, model_name)
get_pcov_file    = lambda HATID  : "%s/%s.pcov"%(LCCACHE, HATID )
get_model_name = lambda iteration : "%s_iter%04d"%(model_prefix, iteration)
get_scores_fname = lambda HATID, iteration  : "%s/%s-%s.scores"%(LCCACHE, HATID, get_model_name(iteration))

get_labeled_hatids_fname = lambda : "%s/labeled_hatids.dat"%(model_output_dir)
get_mystery_hatids_fname = lambda iteration : "%s/uncertain_labels_iter%d.dat"%(model_output_dir, iteration)
get_classifier_fname = lambda iteration : "%s/classifier_iter%04d.pkl"%(model_output_dir, iteration)
get_keylist_dir = lambda field : "~/2007_hatnet_phot/G%s/BASE"%(field)
get_local_keylist_dir = lambda : model_dir
get_local_keylist_fname = lambda field : "%s/keylist_field%s.txt"%(get_local_keylist_dir(), field)


get_remote_2mass_dir = lambda : "~"
get_remote_2mass_fname = lambda field : "%s/colors_field%s.dat"%(get_remote_2mass_dir(), field)
get_local_2mass_dir = lambda : model_output_dir
get_local_2mass_fname = lambda field : "%s/colors_field%s.dat"%(get_local_2mass_dir(), field)


API_KEY = {
    'lcdirect' : 'ZjNmZjQ0NzY4MTQxNzQ0Zjk1OTdlNzY1MTAxOTY1YTQyNDNlMzZlZmE2MWE3M2E3YTY0OWE1MDM5ZDU5NmRjYQ'
}

datafiles = os.listdir(data_dir)
regexp = re.compile("(HAT-...-.......)-hatlc\.csv\.gz")
available_hatids = []
for df in datafiles:
	result = regexp.match(df)
	if result is None: continue
	elif result.groups(1) is None: continue
	else: hatid = result.groups(1)[0]
	available_hatids.append(hatid)
available_hatids = np.array(available_hatids)
#print available_hatids[0]
#print len(available_hatids)

#def get_n_random_available_hatids(n):


phs13_lcdir = "/nfs/phs3/ar1/lcserver/lightcurves"
phs13_list_dt = np.dtype([
		('hatid', 'S15'),
		('ndet', np.int_),
		('relative_path','S200')
	])

phs13_list_file = 'lcs_on_phs13_gcvs.txt'
phs13_list = {}
with open(phs13_list_file,'r') as f:
	head = 0
	for line in f:
		if head < 2:
			head += 1
			continue
		cols = [ l.strip() for l in line.split('|') ]
		phs13_list[cols[0]] = {
			'ndet' : int(cols[1]),
			'relpath' : cols[2]
			}

#print "PHS13 files = %s"%(len(phs13_list))
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
ids_to_use = [ hid for hid in gcvs_m if gcvs_m_types[hid] in types_to_use and os.path.exists(hat_fname(hid)) and not hid in bad_ids and not hid in look_at_again ]



match_cat_dt = np.dtype([
		('id', 'S15'),
		('vartype', 'S10'),
		('dist', np.float_)
	])
vartypes_to_classify = [ 'RRAB','none' ]
skip_vartypes = [ '*' , 'R' , 'RR', 'R(B)']
classify_categories = False
default_clfr = QDA
skip_features = [ 'ra', 'dec', 'V', #'std', 'V'
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
	#max_features=None,#'auto', 
	max_features='auto',
	bootstrap=True, 
	oob_score=False, 
	n_jobs= -1, 
	random_state=None, 
	verbose=0
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

FeatureVectorFileName = lambda hatid, model_name : hat_features_fname(hatid, model_name)
FitCovarianceMatrixFileName = lambda hatid, model_name : "%s/%s-covmat-%s.pkl"%(feat_dir, hatid, model_name)
