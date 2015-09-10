from settings import *
import gzip, sys
import cPickle as pickle
from utils.miscutils import open_ssh_connection, close_ssh_connection, get_2mass_data_for_hatid_over_ssh

hatids = pickle.load(open("%s/good_gcvs_hatids.list"%(parent_dir),'rb'))
hatid_field_list = pickle.load(open("%s/hatid_field_list.pkl"%(parent_dir), 'rb'))
norig = len(hatids)
pruned_hatids = [ hatid for hatid in hatids if hatid in hatid_field_list ]
print norig - len(pruned_hatids), " hatids pruned because they have no known field."
hatids = pruned_hatids

twomass_info_gcvs = {}

fields = np.unique([ hatid_field_list[hatid] for hatid in hatids ])

ssh, sftp = open_ssh_connection()
for i,field in enumerate(fields):
	print " Obtaining info for field %s (%d/%d)"%(field, i+1, len(fields))
	twomass_info = twomass_info_file(field)
	
	field_hatids = [ hatid for hatid in hatids if hatid_field_list[hatid] == field ]
	twomass_hatids = [ hatid for hatid in twomass_info ]
	all_field_ids = [ hatid for hatid in  hatid_field_list if hatid_field_list[hatid] == field ]
	missing_hatids = [ hatid for hatid in field_hatids if not hatid in twomass_hatids ]
	
	print "%d/%d hatids in list of twomass ids for field %s"%(len(twomass_hatids), len(all_field_ids), field)
	print "%d gcvs hatids supposedly in field %s"%(len(field_hatids), field)
	print "%d gcvs hatids have twomass data for field %s"%(len(field_hatids) - len(missing_hatids), field)
	resave= False
	nfound = 0
	for hatid in missing_hatids:

		twomass_info_gcvs[hatid] = get_2mass_data_for_hatid_over_ssh(hatid, ssh)
		if not twomass_info_gcvs[hatid] is None:
			nfound += 1
			resave= True
			twomass_info[hatid] = twomass_info_gcvs[hatid]

	if resave: 
		pickle.dump(twomass_info, gzip.open("%s/twomass_info_field%s.pklz"%(twomass_dir, field), 'wb'))
		print " -- %d twomassinfo files found!"%(nfound)
	#sys.exit()
	for hatid in field_hatids:
		if not hatid in twomass_info:
			continue
		twomass_info_gcvs[hatid] = twomass_info[hatid]

close_ssh_connection(ssh, sftp)
print " Done! we have twomass info for %d out of %d hatids"%(len([ hatid for hatid in twomass_info_gcvs ]), norig)
pickle.dump(twomass_info_gcvs,gzip.open("%s/twomass_info_fieldgcvs.pklz"%(twomass_dir), 'wb'))