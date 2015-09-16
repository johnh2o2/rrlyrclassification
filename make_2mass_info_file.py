# Running this file will
#    * Load in the twomass_info file, which contains twomass info for all hatids 
#	 * Iterates through a list of fields
#			* checks if a color_file for the field is available on phn1
#			* if not,
#				a) Makes a txt file containing list of hatids in that field
#				b) Transfers this to phn1
#				c) Runs a bash script on phn1 which generates a twomass color_fileFIELD.dat file.
#	* Iterates through a list of fields (again)
#			* Downloads the colorfile for that field from phn1
#			* Loads in the colorfile; adds it to the twomass_info dict
#			* Saves the twomass_info dict.
#
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
FORCE_REDO = False

print "loading field list.."
field_list = pickle.load(open("field_info2.pkl", 'rb'))

print "opening ssh connection.."
# Open ssh connection to phn1
client, sftp = open_ssh_connection()

def load_twomass_safely(filename):
	ncols = len(twomass_dt_arr)
	f = open(filename, 'r')

	lines = []
	bad_lines = []
	i=0
	for line in f:
		i+=1
		if '#' in line: continue
		cols = line.replace('\n','').split()
		#print len(cols)
		# Is the number of columns ok?
		if len(cols) != ncols: 
			bad_lines.append((i,line))
			continue

		# Is the value of each column ok?	
		col_is_bad = False
		for j,col in enumerate(cols):
			name, dt = twomass_dt_arr[j]
			try:
				conv(col, dt)
			except:
				col_is_bad = True
				break
		if col_is_bad:
			bad_lines.append((i,line))
			continue

		# OK!
		lines.append(cols)

	f.close()
	arr = np.empty(len(lines), dtype=twomass_dt)
	for i,line in enumerate(lines):
		#print line, len(bad_lines)
		for j, col in enumerate(line):
			colname, dt = twomass_dt_arr[j]
			arr[i][colname] = conv(col, dt)

	return arr, bad_lines

def download_and_add_colorfiles(fields):

	# Now just download the colorfiles from phn1 and load them into the tmd dict.
	for fno, field in enumerate(fields):
		print "   - reading twomass_info file for field %s.."%(field)

		# Read in the twomass_info file.
		fname = "twomass_info_field%s.pklz"%(field)
		tmd = {}
		if os.path.exists(fname):
			tmd_l = pickle.load(gzip.open(fname, 'rb'))
			for hatid in tmd_l:
				tmd[hatid] = tmd_l[hatid]

		print "   - read twomass info for %d hatids"%(len([hatid for hatid in tmd]))

		color_file = "colors_field%s.dat"%(field)
		color_file_l = "%s/%s"%(keylist_dir, color_file)
		color_file_r = "%s/%s"%("/home/jhoffman",color_file)


		# Download
		if not os.path.exists(color_file_l) or FORCE_REDO:
			print "Downloading colorfile for field %s"%(field)
			sftp.get(color_file_r, color_file_l)

		print "   - loading colorfile for field %s (%d/%d)"%(field, fno, len(fields) )

		# Convert
		dat, bad_lines = load_twomass_safely(color_file_l)
		if len(bad_lines) > 0: 
			print "   - %d/%d lines are bad (field %s)"%(len(bad_lines), len(dat) + len(bad_lines), field)

		for d in dat:
			if d['hatid'] in tmd: continue
			tmd[d['hatid']] = {}
			for column, dt in twomass_dt_arr:
				tmd[d['hatid']][column] = conv(d[column], dt)

		# save.
		print "   - saving %d hatids"%(len([hatid for hatid in tmd]))
		pickle.dump(tmd,gzip.open(fname, 'wb'))

	return tmd


def new_proc(field):
	color_file = "colors_field%s.dat"%(field)
	color_file_l = "%s/%s"%(keylist_dir, color_file)
	color_file_r = "%s/%s"%("/home/jhoffman",color_file)
	is_done_file =  "/home/jhoffman/is_done%s.dat"%(field)

	# If it's already done, skip it.
	if rexists(sftp,is_done_file) and not FORCE_REDO:
		if sftp.lstat(is_done_file).st_size > 0:
			return None, None, field
	
	# make file of hatids in the field
	hatids_fname = "%s/hatids_field%s.dat"%(SCRATCH, field)
	if not os.path.exists(hatids_fname) or FORCE_REDO:
		
		all_hatids = [ hatid for hatid in hatid_field_list if hatid_field_list[hatid] == field and not hatid in tmd ]
		print "making hatids file for field %s..."%(field)
		f = open(hatids_fname, 'w')
		for hatid in all_hatids:
			f.write("%s\n"%(hatid))
		f.close()

	# Move that file to phn1
	sftp.put(hatids_fname, "/home/jhoffman/hatids_field%s.txt"%(field))
	
	# Make colorfile
	print "running batch script for field %s..."%(field)
	t0 = time()
	stdin, stdout, stderr = client.exec_command('bash /home/jhoffman/color-file.sh %s %s > is_done%s.dat'%(field, field_list[field], field))
	channel = stdout.channel

	return channel, t0, field

	

def generate_color_files(Fields, nthreads=nthreads):
	fields = [ field for field in Fields ]
	
	procs = []
	while True:

		# Are we done yet?
		if len(fields) == 0 and len(procs) == 0: 
			print "Done!"
			return True

		# Iterate through running processes
		if len(procs) > 0:
			inds_to_pop = []
			i=0
			while i < len(procs):
				channel, t0, field = procs[i]

				# If this field already has a colorfile, act accordingly
				if channel is None: 
					print " field %s already done."%(field)
					procs.pop(i)
					continue

				# If this is finished, note it
				elif channel.exit_status_ready():
					print " field %s finished in %.3f minutes"%(field, (time() - t0)/60.)
					procs.pop(i)
					continue

				else:
					i+=1
					continue
		
		# Add processes
		while len(procs) < nthreads and len(fields) > 0:
			procs.append(new_proc(fields.pop()))

if __name__ == '__main__':

	
	fields = [ field for field in field_list ]
	#fields.append("gcvs")

	generate_color_files( fields )
	tmd = download_and_add_colorfiles( fields )





