import pandas as pd
from math import *
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import readhatlc as rhlc
import gcvs
from astroML.crossmatch import crossmatch_angular
from utils import *
mag_cols = [ 'TF1' ]#, 'TF2', 'TF3' ]
mag_type = 'vmag'
lightcurves = {}
stats = {}
stats['rms'] = {}
#stats['f1'] = {}
#stats['a1'] = {}
#stats['f2'] = {}
#stats['a2'] = {}
#stats['f3'] = {}
#stats['a3'] = {}
for key in stats:
    for mtype in mag_cols:
        stats[key][mtype] = {}

nsrc = 100
step = 1000
max_match_dist = 5. # arcsecond



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
print "Loading list of HAT sources..."
sources = pd.read_csv("../data/hat-field-219-object-list.csv")
sources = sources.iloc[:-1]
sources = sources.sort('vmag')

#print sources.iloc[50:60]['#hat_id']
# Cross-match the two.
print "Cross-matching GCVS and HAT sources..."
hat_radec = np.zeros((len(sources),2), dtype=np.float64)
gcvs_radec = np.zeros((len(gcvs_cat),2), dtype=np.float64)

hat_radec[:,0] = np.array([ ra for ra in sources.loc[:,'ra'] ])
hat_radec[:,1] = np.array([ dec for dec in sources.loc[:,'decl'] ])

gcvs_radec[:,0] = gcvs_cat['RA']
gcvs_radec[:,1] = gcvs_cat['DEC']

gcvs_match_distances, gcvs_match_indices = crossmatch_angular(hat_radec, gcvs_radec,max_distance=1.0)
f_gcvs_match = open("gcvs_matches.catalog","w")
matched_inds = []
for i,d in enumerate(gcvs_match_distances*3600.): 
    if d < max_match_dist: 
        matched_inds.append(i)
        if gcvs_match_indices[i] in np.arange(0,len(gcvs_cat)): 
		f_gcvs_match.write("%s %s\n"%(sources.iloc[i]['#hat_id'],gcvs_cat[gcvs_match_indices[i]]['VarType']))
f_gcvs_match.close()
sys.exit()
nmatches = len(matched_inds)
matched_inds = np.array(matched_inds)
print "%d matches found during cross-matching"%(nmatches)

###
print "Now selecting %d random HAT sources"%(nsrc)
# Pick nsrc random HAT sources from the full list

indices = np.arange(0,len(sources),step)
#np.random.shuffle(indices)
rand_inds = indices[0:nsrc]

indices = np.arange(0,len(matched_inds))
np.random.shuffle(indices)
rminds = indices[0:min([len(indices),nsrc])]

rand_inds = [ ri for ri in rand_inds ]
for mi in matched_inds[rminds]: 
    if not mi in rand_inds: rand_inds.append(mi)
rand_inds = np.array(rand_inds)

hatids = [ hid for hid in sources.iloc[rand_inds]['#hat_id'] ]
#print hatids

# Now get the distance to and type of the closest variable star for each of these HAT sources
match_dists = {}
match_vartypes = {}
magnitudes = {}
for i in rand_inds: 
    HAT_ID = sources.iloc[i]['#hat_id']
    magnitudes[HAT_ID] = sources.iloc[i][mag_type]
    if HAT_ID not in hatids: 
        print "Asking for unknown hatid: %s"%(HAT_ID)
        sys.exit()
    j = gcvs_match_indices[i]
    if j == len(gcvs_cat):
        d = 100*3600.# in arcseconds
        match_dists[HAT_ID] =  d 
        match_vartypes[HAT_ID] = "None"
    else:
        d = gcvs_match_distances[i]*3600.# in arcseconds
        match_dists[HAT_ID] =  d 
        match_vartypes[HAT_ID] = gcvs_cat[j]['VarType']


def get_ref_val(hatid,col):
    ind = -1
    for i in range(0,len(sources)): 
        if sources['#hat_id'][i] == hatid:
            ind = i
            break
    if ind == -1:
        print "ERROR, can't find index for %s"%(hatid)
    return sources[col][ind]



# FETCH lightcurves from the server
print "Fetching lightcurves ... "
fetch_lcs(hatids)

# Calculate stats (rms, frequency fits, etc)
print "Calculating stats..."
skipped_matches = 0
dt_stats = np.dtype([('id','S15'),('rms',np.float_)])
rmsfile = "stats_rms.dat"
if os.path.exists(rmsfile):
    f = open(rmsfile,'r')
    sdata = np.loadtxt(rmsfile,dtype=dt_stats)
    for i,hid in enumerate(sdata['id']):
        if hid not in hatids: continue
        for mcol in mag_cols:
            stats['rms'][mcol][hid] = sdata[i]['rms']
    f.close()

else: sdata = None
for ID in hatids:

    if not sdata is None:
        if ID in sdata['id']: continue
    print "Analyzing %s"%(ID)
    matched_id = False
    if match_dists[ID] < max_match_dist: 
        #print "   [matched]"
        matched_id = True
    fname = "%s/%s-hatlc.csv.gz"%(data_dir,ID)
    lightcurves[ID] = rhlc.read_hatlc(fname)
    #print lightcurves[ID]
    for mcol in mag_cols:
        lc = lightcurves[ID][mcol]
        lc_nonans = [ l for l in lc if not isnan(l) ]
        print "  %.2f%% of %s values are NaN's"%(100*float(len(lc) - len(lc_nonans))/float(len(lc)), mcol)
        if len(lc_nonans) < 2:
            if matched_id: skipped_matches += 1
            print "   (skipping)"
            stats['rms'][mcol][ID] = -1.0
            continue 
        else:
            mu = np.mean(lc_nonans)
            stats['rms'][mcol][ID] = sqrt(np.mean(np.power(lc_nonans - mu,2)))
print "%d out of %d matched sources are all NaN's :("%(skipped_matches,nmatches)


f = open(rmsfile,'w')
mcol = mag_cols[0]
for ID in hatids:
    f.write("%s %e\n"%(ID,stats['rms'][mcol][ID]))
f.close()


rmsarr = {}
mags = {}
colors = {}
def color(hatid):
    if match_dists[hatid] < max_match_dist: return 'r'
    else: return 'k'
okids = [ ID for ID in hatids if stats['rms'][mcol][ID] > 0 ]

print "Making plotting arrays"
for mcol in mag_cols: 
    mags[mcol] = [ magnitudes[ID] for ID in okids ]
    rmsarr[mcol] = [ stats['rms'][mcol][ID] for ID in okids ]
    colors[mcol] = [ color(ID) for ID in okids ]
print "(done)"

f = plt.figure()
ax = f.add_subplot(111)
mcol = mag_cols[0]
ax.set_yscale('log')
ax.scatter(mags[mcol], rmsarr[mcol], c=colors[mcol],alpha=0.5)
ax.set_ylabel("RMS")
ax.set_xlabel("vmag")
f.savefig("hat219_with_gcvs_rms.png")

plt.show()



