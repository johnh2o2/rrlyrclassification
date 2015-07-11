
# June 17, 2015

* Revisiting project; trying to remember where I was at!
* Thoughts : Need a coherent logger structure
* Also need to make things as automated and as intuitive as I can.
* Iterative improvements of the model should be documented, etc.
* Build modular code so that it can be extended to all fields as easily as possible.
* Integrate with the existing python system as much as possible

# June 29-31, 2015

* Working on **label_candidates.py** -- should just integrate with the vislc program you wrote...
	* **done**(?) not tested...
* Link of possible interest: http://www.robots.ox.ac.uk/~szheng/CRFasRNN.html
* **snag** the Visualizer can only work with the Gzipped csv files and not the raw lightcurves.
	* This means that I have to do an awkward `fetch_lcs()` call, ostensibly while on Della, when the data might be right there or accessible in other, faster ways. The `fetch_lcs` call gets the files via an 
	http GET request. This will be MUCH slower than `scp`/`rsync` from the LC server...
		* Talk to Waqas!!!
* Should have a seamless way to just fetch a lightcurve -- the module/function should be able to deal with `.tfalc`, `csv.gz`, etc. Maybe it could just return a standardized dictionary of information
* Debugging the `create_initial_labeled_hatids.py` script; what a nightmare. 

#July 8, 2015

* Porting everything to della
* Problems arising -- need to install fastlombscargle, etc. Maybe there's a global setup.py file that would work?
* Paramiko wasn't installed...I thought it was? I ran `python setup.py install --user` in the `paramiko` directory. That seems to have solved the problem...
* `update_model.py` is NOT parallelized (yet; working on it).
	* Needs to be master-slave for getting features.
	* Needs to NOT just load the csv file from the server (which defeats the entire purpose of what we're doing here...)
	* Memory concerns:
		* 89 features (as of now); assume double precision (8 bytes). There are <~ 6 * 10^6 total hatids. Loading all of the features into a single object will take: (89 features) * (8 bytes/feature) * (6 * 10^6 objects) / 10^9 Gb = 4.272 Gb. So, for individual fields (for the time being) we should be fine on memory (4-8 Gb per core on Della).
* To make things work more parallel, I decided to go with bagging models (#bags = #cores at the moment)
	* Debugging status: things are hanging up in the `masterslave.py` file...still not sure what the cause of this is...
		* This was caused by an empty sublist in the `work_packets` list.
	* Need to extend pickle for I/O
	* 