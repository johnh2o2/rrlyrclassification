
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

* Attempted that last idea; not working *at all* (using an SVM). It's working *so* poorly that I suspect there's a bug in the code...
* Yet *another* idea
	* This is effectively a (4d-ish) problem: 2d kohonen for color + 2d kohonen for LC shape
	* You could probably train
		1. Kmap_color   
		2. Kmap_lcshape 
		3. TSVM(1, 2) -> score

# September 23
* Goal for today:
	* Streamline labeling process
		* Implementing simple process to transfer -> label -> transfer results back. 
		* Done!

# Oct 1
* Modified vislc/`label_candidates.py` to work more smoothly; corrected some annoying cosmetic bugs in `label_candidates.py`
* Label candidates now works on both local and remote machines:
	1. `scp della:/tigress/jah5/rrlyr_scratch/candidates_N.tar .`
	2. `python label_candidates.py --tarfile candidates_N.tar --dont-save --visualize`
	3. `scp candiate_results.dat della:/home/jah5/rrlyr_search/rrlyrclassification/candidate_results_N.dat`
	4. `ssh della`
	5. `cd rrlyr_scratch/rrlyrclassification`
	6. `module load anaconda; module load mpi4py; module load openmpi`
	7. `python label_candidates --save --results candidate_results_N.dat`
* Then to update the model:
	1. `salloc --ntasks=12 --ntasks-per-socket=12 -t 2:00:00` (to run an interactive job for this, it shouldn't take long)
	2. `module load anaconda; module load mpi4py; module load openmpi` (load necessary modules)
	3. `time python update_model2.py` (or n = 1, whatever; number of RF classifiers to train)
* Then to get candidates:
	1. 
* Now update_model.py isn't working...
	* **REASON** -- BECAUSE THIS IS THE OLD VERSION. new version is `update_model2.py`.

# Oct 3
* Now I need to simplify the update_model.
* Bagged random forest classifiers is a completely redundant idea.
* I should experiment with AdaBoost/other classifiers

# Oct 9
* Tried to run get_candidates on a number of fields (~10); allocated 12 (12 procs/node) nodes for 12 hours; waited ~30 hours in the queue, only for it to break in < 10 minutes!!!!
	* "Bus error" was the reason.
* Reran on 145 in interactive mode; min score = 0.4, with p(score > min score) > 0.6. 5 candidates found; 
	* HAT-145-0000507 (RR c); other four are probably sidereal noise.
	* Am I sure that I'm not training the algorithm on RRc's??
		* Need to look through all GCVS RRAB sources again...
* Where are these fields? Show map!
	* 145 is far outside the galactic plane

# Oct 11
* Going to get get_candidates to work today if I have to stay up all night.
* Running on interactive node for a bunch of 2__ fields. So far no problems....hmm....
* Running with 12 nodes/12 proc per node on a 1 hr allocation...

# Oct 13
* Results (fields: 145, 219, 216, 214, 215, 212, 213) 
```
410197 total hatids
98645 bad ids
33 candidates
               : done.

real    215m27.890s
```

* Large number of 'bad ids'. Are we sure that we're using all of the twomass data?
* 12 RRab's found; 
* DT Gem is in this list of candidates: DT Gem is in GCVS as RRc!! Joel confirms that this is RRab (per > 0.6 d, not sinusoidal).
* Updating model
* **Meeting with Joel**:
	* Testing our algorithm:
		1. Inject/recover fake RR-lyrae:
			a. Use large sample (i.e. from OGLE or ASAS (or both))
			b. Use smaller sample from 1 or 2 GC's
				* These should be relatively bias-free, since more attention is paid to each star...
		2. Use regions that have Kepler observations; how many Kepler RR Lyrae do we recover?
	* Adding features?
		* Dust maps to estimate extinction and color corrections (Schlegel and Finkbeiner, e.g.)
		* Reduced Proper Motion
	* Use a simple, traditional "cut" model to see how that does, and if you're missing any.
* Submitted job on ALL fields (12 * 12 for 6 hrs) (I don't want to wait too long for it to run, so I'm going to do things in shorter runs).
	* Broke; some field_ids are None, and miscutils tries to iterate through them (line 92); hopefully fixed this now.
	* resubmitted.

# Oct 17
* 22 candidates from the last run
	* couple of problems during labeling:
		* about 5-10 lightcurves were problematic: looking at HAT-129-0016713-full.tfalc.gz, the BJD column is *all strings*!! added line to pyvislc that casts all of the elements to floats.
		* BUT, in that lightcurve, we have that the len(TF1) > len(BJD)!! Trying other columns for flux gives us a different number.
		* No idea what the hell is going on here.
	* 6 found!!
* Submitted 12 hr job on line 3 of fields in settings.py

# Oct 21
* Didn't update in a little while; on iteration 6 right now.
	* Trying to do more fields at a time.
	* Did a 12hr/12node*12core run yesterday, >3hrs in the queue, only to have python throw an exception (lc was None in load_tfalc i think?).
		* SRUN DIDN'T BREAK.
		* E-mailed CSES with the relevant log information.
			* They are insisting that it didn't crash because I caught the exception, but I *didnt* catch the exception...

* Working on getting simulated RR lyrae;
	* Read in and fit fourier components to M5 LC's
	* Need to fetch a set of random HAT lightcurves now.
		* Evaluate them to see if they're
* It looks like features, when compressed, should take up about 2E3B * 5E6 lcs = 10 GB.
* Should get a current list of all RR lyrae
* FOUND DRAKE ET AL 2013: 12303 HALO RRAB....
	* [Drake 2013](http://cdsads.u-strasbg.fr/abs/2013ApJ...763...32D)
	* Downloaded catalog.
	* How did I crossmatch against GCVS again?
	* (2181 cross-matches!!!!) <-- nvm, code bug; only 341.

* **0 crossmatches with OGLE III** (to be expected).
* What about ASAS?
* How about a catalog of boring stars?
	* Cross-match with Catalina survey!!
		* [Catalina Survey DR2](http://nesssi.cacr.caltech.edu/DataRelease/)

# Oct 22
* Talked with Waqas at end of day yesterday; 
	* Stripe 82 (SDSS) [here](http://www.astro.washington.edu/users/ivezic/sdss/catalogs/S82variables.html)
		* [Standard stars!!](http://www.astro.washington.edu/users/ivezic/sdss/catalogs/stripe82.html)
		* Also [here](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/MNRAS/424/2528)
	* NSVS [here](http://skydot.lanl.gov/)
	* Loneos [here](http://adsabs.harvard.edu/abs/2008ApJ...678..865M)
	* ASAS (sent csv of catalog)
	* PanSTARRS ("though they've been bad at publicly releasing data")
	* Palomar Transient Factory [here](http://www.ptf.caltech.edu/page/DR2), though I don't know how to get to the data
	* 
* Also found [LINEAR](http://www.astro.washington.edu/users/ivezic/linear/PaperIII/PLV.html) catalog(s)
* Checked:
	* OGLE III
	* GCVS
	* Catalina
		* Drake 14 (non-RRab variable stars)
		* Drake 13 (Halo RRLyr)
		* Drake 13 II
	* 
* Look into using the [Supersmoother](http://www.slac.stanford.edu/cgi-wrap/getdoc/slac-pub-3477.pdf) algorithm; [github repo](https://github.com/jakevdp/supersmoother)
* NEO projects from Wikipedia:
	* ATLAS 
	* . Catalina Sky Survey 
	* . LINEAR 
	* . LONEOS 
	* NEAT 
	* NEOSSat 
	* NEOCam 
	* NEODyS 
	* OGS Telescope 
	* Orbit@home 
	* Pan-STARRS 
	* SCAP 
	* Sentinel Space Telescope Sentry 
	* Spacewatch 
	* WISE
* Also:
	* SEKBO [here](http://cdsbib.u-strasbg.fr/cgi-bin/cdsbib?2012ApJ...756...23A) and [here](http://cdsbib.u-strasbg.fr/cgi-bin/cdsbib?2009ApJ...691..306P)
	* Tsinghua University–NAOC Transient Survey

# Nov 13
* Time to get cracking.
* Goal: paper submitted in january 
* Start working on it! MNRAS sounds like it makes the most sense

# Nov 17
* Goals today:
	* E-mail Gaspar to meet about (1) thesis meeting (2) thesis ideas (3) going forward: goals + timeline
	* **[DONE]** Simplify get_candidates script; maybe incorporate processing?
		* Add command line argument capabilities
		* Edit miscutils...what do you need and what do you not need?
		* Edit settings?
	* **[DONE]** Get large della runs going
	* Continue to write code to evaluate completeness
* **lots of itfalc** lightcurves being used....this seems like it's probably a problem.
* Last della run (6356384: 4 hours 12 * 12 cores)
	* Exceptions thrown -- in load_tfalc. We "except ValueError" when reading the lightcurve file but then raise an exception. 
		* This has been fixed on the della repo by instead printing the error and then returning None.
* Submitted another della run (6365411)
	* Aborted (was testing hatids)
* Miscutils is still messy; made a couple of functions to set 
	* `fields_to_analyze` (so we don't have to change the settings file every single time we want to run things)
	* `hatid_field_list` (which relies upon `fields_to_analyze`)
* Added `how_is` function to `~/.bash_profile` on della.
	* Do `how_is $SLURM_FILE` and it will check for Exceptions.
* How do we port everything to Della without substantial delays + breaking everything???
	* Wait until current run is over.
	* Git pull. Resolve (hopefully) minor issues.
	* Interactive run.
	* Run on pickled list of crossmatched sources!!

# Nov 18
* Della run completed: 29 candidates
* Pulled on della. Let's see if we can't get the new get_candidates to work!

# TODO:
* Look at: http://arxiv.org/pdf/1306.6664v2.pdf <-- conditional entropy; 1.5 orders of magnitude faster than LS and about 1 order of mag more effective.
* Do better cross-validation.
	* You should be using the 2d selection function: FPR(pmin, P(p>pmin)), TPR(pmin, P(p>pmin)); pick a minimum threshold, then choose an optimal selection criteria
* **[DONE]** Run current implementation on 145 + 219 ON DELLA.
* **[DONE]** (~17-18 seconds) Generate timing information. How long / lightcurve?
* **[DONE]** Implement SSRFC method (in its own module). Test on some example datasets.
	* Doesn't work.
* Write a general set of convergence functions (or use anything given by sklearn)
	* Want to know:
		1. How does the error rate depend on number of samples?
		2. How does the ROC curve look as a function of iteration?
* Implement a better candidate selection mechanism: You're not taking advantage of the state of the art.
	* Don't want RRLyr-like objects. You want to choose "candidates" in the most efficient way possible to improve the model!
* Attempt to make Kohonen maps of LC shapes faster.
	* The "naive" way (maybe the only way) to train is len(xi) * len(xi) * N^d * Ntrain * Nsamples ~ (15)^2 * (100)^(2) * (5 * 10^6) * 10000 ~ 1.3E17 FLOPs. 
	* Time = FLOPs/(FLOPs/s) = FLOPs / (1-4 FLOPs/cycle * (clockfreq)) ~ (1.3E17 FLOPs) ( 1 cycle / 1 FLOPs) ( 1 s / 2E9 cycles) ( 1 Hr / 3.6E3 s) ~ (1.3/7.2)E4 ~ 1.8E4 computational hours 
	* Parallelize the training/searching!
