## Automated classification of HAT sources



### Overview

The directory setup is a bit overwhelming at the moment, but the basics of getting things up and running are:

*  Run `setup.py` (build/install) to install the Cython LSP algorithm
*  Install any dependencies, including
	* [PyVislc](https://github.com/johnh2o2/pyvislc)
	* numpy
	* matplotlib
	* scipy
	* scikit-learn
	* mpi4py
	* hashlib
	* passlib
	* psycopg2
	* MySQLdb
	* pyfits
	* shlex
	* tornado
	* paramiko
	* spice/PySPICE: [PySPICE](https://github.com/rca/PySPICE), [CSPICE](http://naif.jpl.nasa.gov/naif/toolkit_C_PC_Linux_GCC_64bit.html)

### Install

* Edit the `settings.py` file to suit your needs _this is also quite messy at the moment and will be cleaned up_
* Run `cd utils; python setup.py install`; if you don't have permissions, you can append `--user` to the end of that.
* If you're running this on a remote system, make sure you have X11 forwarding on. Otherwise matplotlib will crash everything. Sorry; kind of annoying.
* 

### Running
1. Run `create_initial_labeled_hatids.py` to label and organize GCVS-crossmatched sources.
2. Run `update_model.py` to generate a classifier based on the labeled sources
3. Run `get_candidates.py` to search a HAT field for possible RR Lyr (or another user-specified subtype of variable star)
4. Run `label_candidates.py` **to be written** to visualize and manually label the candidates
5. Repeat 2 - 4 as many times as necessary until no more new sources are found.


# Tutorials

###Connecting to `della` via terminal:
	* Use `SonicWALL Mobile Connect` to connect to a VPN
	* From terminal, do: `$ ssh della`

###Setting up environment
	* Clone the `rrlyrclassifiation` project to home directory: `git clone git@github.com:johnh2o2/rrlyrclassification.git`
	* [optional] Clone the `pyvislc` project to home directory: `git@github.com:johnh2o2/pyvislc.git`
	* Install additional modules
	* Run `python setup.py install --user` 
