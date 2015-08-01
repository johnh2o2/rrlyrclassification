from utils.miscutils import *
from settings import *
import os, sys
from contextlib import closing
from paramiko import SSHConfig, SSHClient

ids = [ "HAT-173-0000034", "HAT-173-0000016", "HAT-173-0000037" ]


config = SSHConfig()
with open(os.path.expanduser('~/.ssh/config')) as config_file:
    config.parse(config_file)
d = config.lookup(ssh_host_name)


with closing(SSHClient()) as ssh:
	ssh.load_system_host_keys() #NOTE: no AutoAddPolicy() 
	ssh.connect(d['hostname'], username=d.get('user'))
	with closing(ssh.open_sftp()) as sftp:
		results = make_2mass_data_for_field("173", ssh, sftp)
		print results
