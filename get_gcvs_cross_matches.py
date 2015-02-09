from math import *
import numpy as np
import gcvs
from astroML.crossmatch import crossmatch_angular
#from utils import *

max_match_dist = 1. # arcsecond

def ang_diff(a1, a2):
    diff_a = a1 - a2 
    while diff_a < 0: diff_a += 360.0
    while diff_a > 360: diff_a -= 360.0
    return diff_a

def get_dist(ra1,dec1,ra2,dec2):
    # In arcseconds
    # this is only approximate for small angular distances, which is what we care about

    return sqrt(pow(ang_diff(ra1,ra2),2) + pow(ang_diff(dec1,dec2),2))*3600.0

# Load GCVS catalog and HAT sources
print "Loading GCVS catalog..."
gcvs_cat = gcvs.get_conv_cat()
var_types = np.unique( gcvs_cat['VarType'] )
matched_var_types = {}
for v in var_types:
    matched_var_types[v.strip().replace(':','')] = 0



print "Loading list of HAT sources..."
sources = np.loadtxt("/Users/jah5/Documents/Fall2014_Gaspar/complete_catalog.txt", dtype = np.dtype([('id', 'S15'), ('ra', np.float_), ('dec', np.float_)]))


#print sources.iloc[50:60]['#hat_id']
# Cross-match the two.
print "Cross-matching GCVS and HAT sources..."
hat_radec = np.zeros((len(sources),2), dtype=np.float64)
gcvs_radec = np.zeros((len(gcvs_cat),2), dtype=np.float64)

hat_radec[:,0] = sources['ra'][:]
hat_radec[:,1] = sources['dec'][:]

gcvs_radec[:,0] = gcvs_cat['RA']
gcvs_radec[:,1] = gcvs_cat['DEC']

gcvs_match_distances, gcvs_match_indices = crossmatch_angular(hat_radec, gcvs_radec,max_distance=1.0)
huhs = 0
matched_inds = []

for i,d in enumerate(gcvs_match_distances*3600.): 
    if d < max_match_dist: 
        if gcvs_match_indices[i] in np.arange(0,len(gcvs_cat)): 
            matched_inds.append(i)
        else:
            print "Huh..."
            huhs += 1
print "Huston, we have %d problems."%(huhs)
nmatches = len(matched_inds)
matched_inds = np.array(matched_inds)
print "%d matches found during cross-matching"%(nmatches)

for i in matched_inds:
    matched_var_types[gcvs_cat[gcvs_match_indices[i]]['VarType'].strip().replace(':','')] += 1

matched_var_types_arr = np.zeros(len(matched_var_types), dtype=np.dtype([('type', 'S20'), ('num', np.int_)]))
for i,v in enumerate(matched_var_types):
    matched_var_types_arr[i]['type'] = v
    matched_var_types_arr[i]['num'] = matched_var_types[v]

matched_var_types_arr = np.sort(matched_var_types_arr,order='num')[::-1]

nodets = []
dets = []
for v in matched_var_types:
    if matched_var_types[v] > 0: dets.append(v)
    else: nodets.append(v)

#for i,v in enumerate(matched_var_types_arr['type']):
#    if v in dets:
#        print "%-20s %-10d"%( v, matched_var_types[v] )
#print nodets

full_cat = open("Full_GCVS_Cross-matched_HATIDS_maxdist%.2f.catalog"%(max_match_dist),'w')
full_cat.write("# HAT_ID    VAR_TYPE    MATCH_DIST_ARCSEC\n")
for i in matched_inds:
    #full_cat.write("%-20s %-20s %.5e\n"%(sources[i]['id'], gcvs_cat[gcvs_match_indices[i]]['VarType'].strip().replace(':',''), gcvs_match_distances[i]*3600.))
    vt = gcvs_cat[gcvs_match_indices[i]]['VarType'].strip().replace(':','')
    if vt == 'RCB' : print "RCB: ", sources[i]['id']
    if vt == '' or vt == ' ': vt = "NaN"
    full_cat.write("%s %s %.5e\n"%(sources[i]['id'], vt, gcvs_match_distances[i]*3600.))
full_cat.close()

#matched_sources = sources[matched_inds]

