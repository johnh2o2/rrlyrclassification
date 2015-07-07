import os, sys, glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import utils.feature_selection as fs
import utils.readhatlc as rhlc
from utils.miscutils import *
from utils.featureutils import *
from settings import *
import cPickle as pickle
from contextlib import closing
from paramiko import SSHConfig, SSHClient
from mpi4py import MPI

config = SSHConfig()
with open(os.path.expanduser('~/.ssh/config')) as config_file:
    config.parse(config_file)
d = config.lookup(hostname)

field_info = pickle.load(field_info_fname)



def Generate2mass(field_name):


def ListLCDir(remote_path, stfp):
    sftp.chdir(remote_path)
    return stfp.listdir()

def remote_path_exists(path, stfp):
	try:
		stfp.stat(path)
		return True
	except IOError, e:
		return False

def LoadRemoteKeylist(field_name, stfp):
	fname = remote_keylist_fname(field_name)
	
	if not remote_path_exists(fname, stfp): return None
	
	stfp.get(fname, '%s/keylist_field%s.txt'%(keylist_dir, field_name)

	return np.loadtxt(fname, dtype=keylist_dt)

def FindLightcurvePathOnDella(ID,MPI_COMM):
	with closing(SSHClient()) as ssh:
	    ssh.load_system_host_keys() #NOTE: no AutoAddPolicy() 
	    ssh.connect(d['hostname'], username=d.get('user'))
	    with closing(ssh.open_sftp()) as sftp:
			for field in field_info:
				lcpath = "%s/%s.tfalc"%(field_info[field]['path'], ID)
				if remote_path_exists(lcpath, stfp): return lcpath
	return None
