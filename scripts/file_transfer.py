
import subprocess, os, shlex, sys
from time import time
from settings import *

nthreads = 4
RSYNC_RSH="ssh -c arcfour -o Compression=no"
#remote_server="jah5@phn1.astro.princeton.edu"
remote_server="phn1"

def new_proc(rdir, ldir):
	print "Starting process to move %s to %s"%(rdir, ldir)
	if ldir[-1] == '/':
		name = ldir[:-1].split('/')[-1]
	else:
		name = ldir.split('/')[-1]
	
	full_ex = "rsync -aurvhW  --include='./' --include='*.tfalc' --exclude='*' --stats -e '%s' %s:%s %s"%(RSYNC_RSH, remote_server, rdir, ldir)
	logfile = open('rsync-%s.log'%(name), 'w')
	
	
	ex = shlex.split(full_ex)

	proc = subprocess.Popen(ex, stdout=logfile, stderr=logfile)
	print "  --> started %s to %s"%(rdir, ldir)

	return proc, ldir

def transfer_gcvs():

	print "GCVS TRANSFER"
	print "importing..."
	from utils.miscutils import open_ssh_connection, close_ssh_connection, get_field_of, rexists
	print "loading hatids.."
	hatids 				= pickle.load(open("%s/good_gcvs_hatids.list"%(information_dir), 'rb'))
	print "opening ssh connection.."
	client, sftp = open_ssh_connection()
	for i,hatid in enumerate(hatids):
		print hatid, " (GCVS) %d/%d"%(i+1, len(hatids))
		field = get_field_of(hatid)
		if field is None: 
			print " ++++++ ",hatid, " (GCVS) %d/%d"%(i+1, len(hatids)), " field is None!"
			continue
		directory = "%s/%s"%(LCCACHE, field)

		fname = "%s.tfalc"%(hatid)
		fname_r = "%s/%s"%(field_list[field], fname)
		fname_l = "%s/%s"%(directory, fname)
		if os.path.exists(fname_l): continue
		if not rexists(sftp, fname_r ): 
			print "-----  ",hatid, " (GCVS) %d/%d"%(i+1, len(hatids)), "is not available on phn1"
			continue
		if not os.path.isdir(directory): os.makedirs(directory)
		print "Getting ",hatid
		sftp.get(fname_r, fname_l)

	close_ssh_connection(client, sftp)


def transfer(remote_dirs, local_dirs, nthreads=nthreads):
	assert(len(remote_dirs) == len(local_dirs))

	procs = []
	while True:

		# Are we done yet?
		if len(remote_dirs) == 0 and len(procs) == 0: 
			print "Done!"
			return True

		# Iterate through running processes
		if len(procs) > 0:
			inds_to_pop = []
			for i in range(len(procs)):
				proc, ldir = procs[i]
				rc = proc.poll()
				# If this is finished, note it
				if not rc is None: 
					print "==> %s has been copied! :)"%(ldir)
					inds_to_pop.append(i)
			# pop finished processes from the procs array
			for i in inds_to_pop:
				procs.pop(i)
		
		# Add processes
		while len(procs) < nthreads and len(remote_dirs) > 0:
			procs.append(new_proc(remote_dirs.pop(), local_dirs.pop()))

if __name__ == '__main__':
	import cPickle as pickle
	field_list = pickle.load(open(field_info_fname, 'rb'))
	#fields = [ field for field in field_list ]
	#fields =fields_to_analyze
	fields = [ 'gcvs' , '145']
	if 'gcvs' in fields: transfer_gcvs()
	ldir_pfix = "/tigress/jah5/rrlyr_scratch/LCCACHE"
	#ldir_pfix = "testing_file_transfer"
	#ldir_pfix = "/Users/jah5/Documents/Fall2014_Gaspar/rrlyr_classification/SCRATCH/LCCACHE"

	remote_dirs = [ "%s/"%(field_list[field]) for field in fields if not field == 'gcvs' ]
	local_dirs = [ "%s/%s"%(ldir_pfix, field) for field in fields if not field == 'gcvs']

	transfer(remote_dirs, local_dirs)

