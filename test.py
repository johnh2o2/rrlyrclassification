# Tests that mpi_pool.py behaves as expected
# when multiple nodes are involved.
# John Hoffman Feb. 9 2015
from sys import exit
import mpi_pool as mpip
import numpy as np
from math import *
import mpi4py.MPI as mpi

comm = mpip.Comm()
rank = comm.rank
#print rank
#if rank == 0:
NX = 100

wrapper_function = mpip.poolMap


#@wrapper_function
def test_func(cache, x):
	return sum(x)

p = mpip.startPool()

if rank == 0:
	X = [ ]
	for i in range(NX):
		X.append(np.ones(10))
	
	result = p.map('default', test_func, X)
	print result
p.exit()
