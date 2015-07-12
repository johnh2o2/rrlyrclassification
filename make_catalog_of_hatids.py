import os, sys
import numpy as np
import cPickle as pickle
from contextlib import closing
from paramiko import SSHConfig, SSHClient
from mpi4py import MPI
import masterslave as msl
from settings import *
from utils.miscutils import logprint

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
ROOT = (rank == 0)

field_info = pickle.load(open('field_info.dict', 'rb'))
all_fields = [ F for F in field_info ]

def rexists(sftp, path):
    """os.path.exists for paramiko's SCP object
    """
    try:
        sftp.stat(path)
    except IOError, e:
        if e[0] == 2:
            return False
        raise
    else:
        return True

def get_field_ids(field, sftp):
	logprint(" getting hat ids in field %s"%(field), all_nodes=True)
	remote_path = field_info[field]['path']
	if not rexists(sftp, remote_path):
		logprint(" field %s is not where we think it is (%s)"%(field, remote_path), all_nodes=True)
		return { field : None}

	sftp.chdir(remote_path)
	fdict = {field : []}
	files = sftp.listdir()

	for f in files:
		# Skip non-tfalc lightcurves
		if not 'tfalc' in f: continue

		# if more than 1 period in filename, skip.
		if len(f.split('.')) > 2: continue
		
		fdict[field].append(f.replace('.tfalc', ''))
	return fdict

def print_list(l):
	pstr = ""
	for i,item in enumerate(l):
		if i == len(l) - 1: pstr = "%s%s"%(pstr, item)
		else: pstr = "%s%s,"%(pstr, item)

	return pstr
	
def make_hatid_lists(fields = None):
	if fields is None: 
		fields = []
		for field in all_fields:
			fname = get_hatids_in_field_fname(field)
			if not overwrite and os.path.exists(fname): continue
			fields.append(field)
		if len(fields) == 0: 
			logprint("No fields are missing!")
			return None

	if ROOT:
		all_field_infos = msl.master(fields)

		# save information.
		bad_fields = []
		for fi in all_field_infos:
			for field in fi:
				if fi[field] is None: 
					bad_fields.append(field)
					continue
				fname = get_hatids_in_field_fname(field)
				pickle.dump(fi[field], open(fname, 'wb'))
		logprint("%d out of %d fields (%d total fields) have invalid paths!"%(len(bad_fields), len(fields), len(all_fields)))
		logprint(print_list(bad_fields))

	else:
		# SSH connection
		config = SSHConfig()
		with open(os.path.expanduser('~/.ssh/config')) as config_file:
		    config.parse(config_file)
		d = config.lookup(ssh_host_name)

		# Get list of hatids for each 
		with closing(SSHClient()) as ssh:
			ssh.load_system_host_keys() #NOTE: no AutoAddPolicy() 
			ssh.connect(d['hostname'], username=d.get('user'))
			with closing(ssh.open_sftp()) as sftp:
				msl.slave(lambda field : get_field_ids(field, sftp))

def get_field_of(hatid):
	# exhaustive search. slower than it needs to be.
	for field in all_fields:
		fname = get_hatids_in_field_fname(field)
		hatids = pickle.load(open(fname, 'rb'))
		if hatid in hatids: return field
	return None

def get_remote_lc_filename(hatid):
	field = get_field_of(hatid)
	return "%s/%s.tfalc"%(field_info[field]['path'], hatid)



if __name__ == "__main__": make_hatid_lists()
	