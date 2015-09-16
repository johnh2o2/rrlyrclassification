# Downloads all of the keylist files from phn1

from utils.miscutils import *
from settings import *
import gzip
import cPickle as pickle

ssh, sftp = open_ssh_connection()

#fields = np.unique([ hatid_field_list[hatid] for hatid in hatid_field_list ])
fields = [ '145' ]
bad_fields = []
for field in fields:
	print "%s"%(field)
	klfname_r = "%s/keylist.txt"%(get_keylist_dir(field))
	klfname_l = get_local_keylist_fname(field)
	klfname_d = get_local_keylist_dict_fname(field)

	if os.path.exists(klfname_d): continue

	# Download if the keylist doesn't exist.
	if not os.path.exists(klfname_l):
		if not rexists(sftp, klfname_r): 
			logprint("Cannot find keylist %s for field %s on remote system."%(klfname_r, field))
			bad_fields.append(field)
			continue


		klfilesize = sftp.stat(klfname_r).st_size
		logprint("               : remote keylist.txt filesize "+`klfilesize`)
		
		sftp.get(klfname_r, klfname_l)

	# If that failed, there's no keylist...
	if not os.path.exists(klfname_l):
		logprint("Warning: keylist data for field %s is not available!"%(field))
		bad_fields.append(field)
		continue

	# Now try to carefully open keylist.
	try:
		keylist_data = safe_open_keylist(klfname_l)
	except:
		logprint(" %s: Can't even open keylist data with safe_open_keylist..."%(field))
		bad_fields.append(field)
		continue

	# Now try to make a dictionary out of that keylist
	try:
		keylist = make_keylist_dict(keylist_data)
		logprint(" %s: Made keylist dict."%(field))
	except:
		logprint(" %s: Can't make keylist dict :(")
		bad_fields.append(field)
		continue

	# Now try to save that dictionary!
	try:
		f = gzip.open(klfname_d, 'wb')
		pickle.dump(keylist, f)
		logprint(" %s: Saved gzipped keylist dict."%(field))
		f.close()
	except:
		logprint(" %s: Can't save keylist dict :(")
		bad_fields.append(field)
		continue

close_ssh_connection(ssh, sftp)

print "done."


