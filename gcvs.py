from math import *
import numpy as np
import pandas as pd
import re, sys, os

gcvs_dir = "/Users/jah5/Documents/Fall2014_Gaspar/gcvs_cat"

gcvs_file = gcvs_dir + "/gcvs_cat.dat"
gcvs_conv_file = gcvs_dir + "/gcvs_conv_cat.dat"

REGEXP = "([0-9]{6})"			# VarNum
REGEXP = "%s(.)\|"%(REGEXP)		# m_VarNum
REGEXP = "%s(.{10})"%(REGEXP)	# GCVS
REGEXP = "%s(.)\|"%(REGEXP)		# n_GCVS
REGEXP = "%s([0-9]{2})"%(REGEXP)		# RAh
REGEXP = "%s([0-9]{2})"%(REGEXP)		# RAm
REGEXP = "%s([0-9\.]{4})"%(REGEXP)	# RAs
REGEXP = "%s([\+0-9]{3})"%(REGEXP)	# DEd
REGEXP = "%s([0-9]{2})"%(REGEXP)	# DEm
REGEXP = "%s([0-9]{2})"%(REGEXP) 	# DEs
REGEXP = "%s(.)\|"%(REGEXP)		    # u_DEs
REGEXP = "%s(.{10})\|"%(REGEXP)		# VarType
REGEXP = "%s.*"%(REGEXP)
vsre = re.compile(REGEXP)
nentries = 0
lines = []
with open(gcvs_file,'r') as f:
	for line in f: 
		if vsre.search(line):
			lines.append(line)
		
nentries = len(lines)
dt_gcvs_red = np.dtype([
	('VarNum',np.int_),
	('m_VarNum','S1'),
	('GCVS','S10'),
	('n_GCVS','S1'),
	('RAh',np.float_),
	('RAm',np.float_),
	('RAs',np.float_),
	('DEd',np.float_),
	('DEm',np.float_),
	('DEs',np.float_),
	('u_DEs','S1'),
	('VarType','S10')
])

dt_gcvs_conv = np.dtype([
	('VarNum',np.int_),
	('GCVS','S10'),
	('RA',np.float_),
	('DEC',np.float_),
	('VarType','S10')
])
def get_ra_dec(entry):
		RA = entry['RAh']*15.0 + entry['RAm']*0.25 + entry['RAs']*(1.0/240.0)
		DEC = entry['DEd'] + entry['DEm']*(1./60.) + entry['DEs']*(1./60.)*(1./60.)
		return RA, DEC

def get_conv_cat():
	reduced_catalog = np.zeros(nentries, dtype=dt_gcvs_red)
	converted_catalog = np.zeros(nentries, dtype=dt_gcvs_conv)
	

	for i,line in enumerate(lines):
		mat = vsre.search(line)
		

		reduced_catalog[i]['VarNum'] = int(mat.group(1))
		reduced_catalog[i]['m_VarNum'] = mat.group(2)
		reduced_catalog[i]['GCVS'] = mat.group(3)
		reduced_catalog[i]['n_GCVS'] = mat.group(4)
		reduced_catalog[i]['RAh'] = float(mat.group(5))
		reduced_catalog[i]['RAm'] = float(mat.group(6))
		reduced_catalog[i]['RAs'] = float(mat.group(7))
		reduced_catalog[i]['DEd'] = float(mat.group(8))
		reduced_catalog[i]['DEm'] = float(mat.group(9))
		reduced_catalog[i]['DEs'] = float(mat.group(10))
		reduced_catalog[i]['u_DEs'] = mat.group(11)
		reduced_catalog[i]['VarType'] = mat.group(12)
		#print get_ra_dec(reduced_catalog[i])
	for i,entry in enumerate(reduced_catalog):
		converted_catalog[i]['VarNum'] = reduced_catalog[i]['VarNum']
		converted_catalog[i]['GCVS'] = reduced_catalog[i]['GCVS']
		converted_catalog[i]['RA'],converted_catalog[i]['DEC'] =  get_ra_dec(reduced_catalog[i])
		converted_catalog[i]['VarType'] = reduced_catalog[i]['VarType']
	return converted_catalog

if not os.path.exists(gcvs_conv_file):
	converted_catalog = get_conv_cat()
	f = open(gcvs_conv_file,'w')
	for s in converted_catalog:
		f.write('%d,%s,%.10f,%.10f,%s\n'%(s['VarNum'], s['GCVS'], s['RA'], s['DEC'], s['VarType']))
	f.close()
else:
	converted_catalog = np.loadtxt(gcvs_conv_file,dtype=dt_gcvs_conv,delimiter=',')
    