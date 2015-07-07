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

1. Edit the `settings.py` file to suit your needs _this is also quite messy at the moment and will be cleaned up_
2. Run `create_initial_labeled_hatids.py` to label and organize GCVS-crossmatched sources.
3. Run `update_model.py` to generate a classifier based on the labeled sources
4. Run `get_candidates.py` to search a HAT field for possible RR Lyr (or another user-specified subtype of variable star)
5. Run `label_candidates.py` **to be written** to visualize and manually label the candidates
6. Repeat 3 - 5 as many times as necessary until no more new sources are found.

