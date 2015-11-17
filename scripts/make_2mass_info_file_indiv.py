# Running this file will
#	Build twomass_info dictionary for all HAT-IDs (instead of going field by field)
# NOT DONE.
print "importing miscutils.."
from utils.miscutils import *
print "importing settings.."
from settings import *
print "importing everything else.."
import cPickle as pickle
from time import time
import subprocess, os, shlex, sys, gzip, sys


nthreads = 4 # Number of threads to use on phn1.
fields = 'all' 
fname = "twomass_info.pklz" 
FORCE_REDO = True

print "loading field list.."
field_list = pickle.load(open("field_info2.pkl", 'rb'))

print "opening ssh connection.."
# Open ssh connection to phn1
client, sftp = open_ssh_connection()



if __name__ == '__main__':

	twomass_info = {}
	hatid_field_list = load_hatid_field_list(all_fields)

	for field in all_fields:
		twomass_info = {}


	get_2mass_data_for_hatid_over_ssh(hatid, sftp)
	fields = [ field for field in field_list ]
	#fields.append("gcvs")

	generate_color_files( fields )
	tmd = download_and_add_colorfiles( fields )





