import os, sys
from utils import *
from featureutils import *
from settings import *
import cPickle as pickle

# Make the directory for this model if it doesn't already exist
if not os.path.exists(model_output_dir): os.makedirs(model_output_dir)

# Get the filename for the labeled HAT ids.
fname = get_labeled_hatids_fname()

# Check if the labeled HAT ids file already exists. If so, grumble and exit
if os.path.exists(fname):
	print "%s already exists. You should delete %s if you're sure you want to overwrite it."%(fname, fname)
	sys.exit()

gcvs_crossmatches 					= np.loadtxt(full_gcvs_match_cat_name, dtype=match_cat_dt)
ids_to_fetch 						= GetHATIDsToFetch(gcvs_crossmatches['id'])

categories 							= GetCategoriesForEachHATID(gcvs_crossmatches)
good_gcvs_hatids_fname = "%s/rrlyr_classification/good_gcvs_hatids.list"%(parent_dir)
print "Pruning out bad ID's"
if overwrite or not os.path.exists(good_gcvs_hatids_fname):
	hatids 								= GetGoodHATIDs(gcvs_crossmatches['id'], categories)
	pickle.dump(hatids, open(good_gcvs_hatids_fname, 'wb'))
else:
	hatids = pickle.load(open(good_gcvs_hatids_fname, 'rb'))

inds = [ i for i in range(len(gcvs_crossmatches)) if gcvs_crossmatches[i]['id'] in hatids ]
gcvs_crossmatches 					= GCVS_GetRandomSample(num, gcvs_crossmatches[inds])

print "TOTAL HATID's SELECTED: %d"%(len(hatids))

for vclass in vartypes_to_classify: 
	print vclass, sum([ 1 for ID in hatids if categories[ID] == vclass ])

# Write initial labeled ids!
f = open(fname, 'w')
for ID in hatids:
	f.write("%-20s%-10i%-20s\n"%(ID, 0, categories[ID]))
f.close()

