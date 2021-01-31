#!/usr/bin/env python
#
#  Taken from:   https://scicomp.aalto.fi/triton/examples/python/python_openmp/python_openmp/
#  E.g.
#  export OMP_PROC_BIND=true
#  export OMP_NUM_THREADS=4

import os
from time import time
import numpy as np

print('SLURM_CPUS_PER_TASK: Using %d processors' % int(os.getenv('SLURM_CPUS_PER_TASK',1)))

nrounds = 5
n = 2000

t_start = time()

for i in range(nrounds):
    a = np.random.random([n,n])
    a = a + a.T
    b = np.linalg.pinv(a)

t_delta = time() - t_start

print('Seconds taken to invert %d symmetric %dx%d matrices: %f' % (nrounds, n, n, t_delta))
