
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
	* Did this; 

#July 27, 2015

* Making `create_initial_labeled_hatids.py` download the lightcurves. In order to do this, we need a function that *generates* twomass files.
	* Done; just uses binary twomass file to get twomass data for individual hatids.
* t-SNE **LOOK AT THIS!!!** https://www.youtube.com/watch?v=RJVL80Gg3lA&list=UUtXKDgv1AVoG88PLl8nGXmw

#August 4, 2015

* Problem with ssh'ing into phn1 on della; it's because there was some change on della and I'm now jah5@della4 not jah5@della3; made a new dsa key, added to `authorized_keys` on `phn1`; can log in without password now
* new error:

```bash
ssh.connect(d['hostname'], username=d.get('user'))
  File "build/bdist.linux-x86_64/egg/paramiko/client.py", line 237, in connect
socket.gaierror: [Errno -3] Temporary failure in name resolution
```

* **SNAG:** they do NOT allow outside connections on the compute nodes (makes sense, security reasons)
	* Will need to write an rsync script to get things moved from phn1 to della;
	* Space: assume 2Mb per tfalc lc; 5*10^6 lc's = 10 Tb of space (!!!)
	* I have a 10 Tb allocation on /tigress...275Gb are taken up by my N-body sims...might want to delete them?
	* A message from bill:
```
Della and Tiger have their own scratch spaces - /scratch/gpfs which
is parallel, not backed up and seen by all nodes including tigressdata;
/scratch/network which is serial, not backed up and seen by all nodes;
and each compute node's individual /scratch.

The default quota on /scratch/gpfs is also 500G. How much data do
you need to pull in? Which cluster are you using?
```

#August 12, 2015
* Need to write script that transfers data via `rsync` to **tigress**;
	* Run in parallel (4 threads at once?)

* E-mail from Waqas:

```txt
Hi John,

Direct NFS to della/tigress from our machines is probably out of the
question (but I would ask if they'd allow us to mount their scratch
space). rsync over ssh if that's allowed might be your best option. I
wouldn't use a big single rsync though; you'll need to parallelize it by
breaking up by directories (perhaps one rsync process per LC directory),
with a maximum of 4 or so running at the same time. We only have 1-Gbit
uplinks from our racks to the rest of campus (we should ask to upgrade
this too), so this will take some time.

- see this for ideas:
https://wiki.ncsa.illinois.edu/display/~wglick/Parallel+Rsync
- pull from their side rather than push from our side
- use weak SSH ciphers to speed things up
- use the --whole-file flag (-W) to skip rsync delta checking (this
slows things down a lot)
- note that rsync is terrible at lots of tiny files
- do an ssh-copy-id command from their end to ours (probably on phn1) to
copy over your della account SSH key to our server so no password
prompts will be required

Here's an rsync over ssh command I've used to transfer things around at
70 MB/sec or so:

on della/tigress:

$ nohup rsync -aurvhW --stats --rsh="ssh -c arcfour -o Compression=no"
phnX:/source/dir/path/ /path/to/destination/ > rsync-process-XX.log &

rsync with these options over NFS can do a bit better (around 100
MB/sec) with our 1-Gbit links.
```

* Note to self: **delete `.tfalc` files once a `-full.tfalc` file is made**
* Rsync pull: `rsync [OPTION...] [USER@]HOST:SRC... [DEST]` (from http://linux.die.net/man/1/rsync)
* Thoughts
	* Avoid mpi, no scheduler so this isn't necessary
	* Use python (i.e. the `subprocess` module) to spawn rsync commands. (use `subprocess.check_call`, which waits for the command to complete or throws error!)
* Script to transfer data written (and a script to get directory sizes written for use on `phn1`)
* **THERE ARE ~8.2 TB OF DATA**
	* So, enough that I can transfer EVERYTHING to tigress.
* Currently pulling from tigressdata:
	* Transfer rate ~ 110 MB/s (running on 4 threads...sooo I assume this is about 30 MB/s per thread)
* It looks like `*-full.tfalc` files are ~3 times the size of the original `tfalc` files...this is a problem. 
	* Solved if we use gzip to compress/decompress files. *smaller* than original LC.

# August 17

* Wrote script to download twomass information for all hatids; runs `nthreads` processes on `phn1`. This finished last night!
* Strange problem when trying to create twomass info file for gcvs sources -- almost every gcvs crossmatched source was missing from the twomass info files..
	* I was worried that this meant I had somehow deleted all of the gcvs tfalc files on phn1...
	* I checked to see if this was true by testing one of the missing hatids:
```
jhoffman@phn1:/nfs/phn15/ar1/H/4KRED/4KAP_LC/062_PR5_D20$ ls *HAT-094-0001548*
HAT-094-0001548.epdlc  HAT-094-0001548.epdlog  HAT-094-0001548.rlc  HAT-094-0001548.tfalc
```
	* So, it seems like this is *not* the case; somehow the gcvs sources were removed from the twomass info files. How??!?
	* The twomass info for these sources IS there if you run the `2massread` binary on them.
	* Wrote script to take results from `2massread` binary and put them back in the twomass info files... hopefully this should solve the problem.

# August 23
* Came across lots of empty `hatids_in_field[].list` files; e.g. 258.
	* Empty directory in corresponding tigress directory
	* Empty phn1 directory (`/nfs/phn1/ar2/EDRIVE/EDRIVE3/LC/258`).
		* **ALL GOOD**.
* Should I make a "pseudo"-field called gcvs and move all gcvs sources to that directory? That way there's no ambiguity and I think everything would be a lot easier...then again maybe it's best to preserve the filestructure?

# August 25
* OK, the following seems to work fine on my mac:
	`mpirun -n 4 python create_initial_labeled_hatids.py`
	`mpirun -n 4 python update_model.py`
* The PROBLEM is that `get_candidates.py` can't load the model. I need to make sure
	1. That the `update_model.py` script ends by training the model on ALL hatids
	2. That the `update_model.py` script saves this bagged model in a consistent way that can be loaded later on.
* I also need to document `update_model.py` via comments...that thing is a bit cryptic.

#September 10
* Saved bagged models, but these are not trained on the COMPLETE set of labeled data (basically I save them in the middle of cross validation)
	* This will need to be changed, but I'm just trying to get the pipeline to WORK for now.
* The 145 base in the /jhoffman/G145 directory (or whatever) has a lot of bad symlinks. I had picked this field arbitrarily, and the keylist couldn't be found since the BASE directory was not there!! The correct location is "/nfs/phn1/ar2/EDRIVE/EDRIVE51/ANL/145/BASE", and I fixed this sloppily by adding an ``if`` condition in the ``get_remote_keylist_dir`` function in ``settings.py``.
	* **This needs to be looked into more** : how many keylists are missing because of this?
* ``score_hatids`` function was outdated and only used other/mag models and not composite models. I 
* Ran the ``get_candidates.py`` script on the first 10 hatids from 145; ALL OF THEM HAVE SCORES OF **EXACTLY** 0.5!!!!
	* This is obviously a problem.
	* Tested with GCVS sources, STILL A PROBLEM
	* **realized** generate_features doesn't **return** feature vector. Made ``feature_vector`` function that does.
	* **NOT THE PROBLEM** generate_features doesn't have to return feature vector. 
	* It's likely because you aren't giving the model the appropriate translation of the feature dict.
	* **FOUND THE PROBLEM**
		* this is a model-persistence issue.
		* Doing a save/load in ``update_model.py`` generated the exact same problem.
* Trying to use ``dill`` python library to see if that fixes the issue (this will likely need to be installed on Della)
* Problem is that the composite_prediction function was not doing what it was supposed to.
	* Loading of models, scalers, SVM's are fine (except for maybe the SVM's)
	* Doing away with the SVM's in favor of just a simple mean score don't appreciably change things.

#September 11
* Everything runs now, but the model doesn't seem to perform as well as advertised by cross-validation...
	* Even doing away with MC and using the GCVS sources, the model only picks up something like 3/4 of the GCVS RRlyr

#September 17
* Tried running on Della: 
	* `hatids_in_field_().list` files don't exist. Writing code in miscutils to generate them.
		* Wrote something in miscutils, but realized the issue is as simple as scp'ing the folder to della...
	* Also had to copy "other_data" file over to della -- need a simpler system for this!!
	* Needed to do "pip install dill --user"

#September 18 -- meeting notes
* Told them about doing 145 and showed them the sample RR-Lyrae that I found
	* It *could* be an RRab, but might be RRc or contact binary
* Gaspar wants to know timing
* Gaspar pointed out that it might be useful to use as much existing data for a given source as is available
	* ASAS, Hipparcos, etc.
* **Proper motions** -- Joel mentioned that there's some sort of "corrected"(?) proper motion value, and that this gives a reasonable distance estimate (!).
* **Multiple sources**
	* Overlapping fields
	* HAT-144-... might be in field 145 -- this means that the *primary* location of the source is in field 144, but the observed field is 145
	* SHOULD ALL BE NAMED THE SAME THING
	* *Chelsea figured out how to stitch these together* -- you have to detrend them. Basically add a free parameter to your fit!
	* Ignoring this for now; smaller amplitude sources might be missed, but that's ok (for now)
* Worry about how robust your cross-validation is
	* You're not using a representative sample
* Gaspar once-again cautioned that the project can very quickly become too technical
* Gaspar talked about how it might be very useful to involve some sort of citizen science

# September 20
* Another idea that might be more efficient and maybe will achieve the same results (or similar results):
	* Instead of *manually* labeling a flagged subset of the unlabeled data at each iteration, what if we *automatically* label sources that are > some threshold as being RR lyrae? Then retrain, etc.
	* We could re-evaluate scores and then *unlabel* some data that falls below the threshold...
* Yet another idea:

``` python

		# First scores are from initial model
		scores = { source : [ initial_model.score(source) ] }

		while not in_equilibrium(scores):
			# update the scores
			for source in sources:

				# train set of models with stochastic labels
				# on the other sources
				Models = []

				for n in range(N_monte_carlo):
					labels, features = {}, {}

					# Add stochastic labels
					for source2 in sources:
						if source2 == source: continue

						a = random(0,1)
						if a > scores[source][-1]: labels[source] = 'RRlyr'
						else: labels[source] = 'none'

					# train this model
					model = train( features, labels )
					Models.append(model) 

				# new score is the mean of the model scores
				scores[source].append( mean([ model.score(source[i]) for model in Models]) )

```


# TODO
* Collect all relevant non-lc files into a single tarball
* Write code to get `hatids_in_field_gcvs.list`
* Test all aspects of the pipeline.
* Write script to setup directory structure once scratch is erased. (This should also work on remote systems too)