#! /usr/bin/env python
#
#  Comparing python and NEMO/C in reading a (large) table, and performing some
#  basic stats on all columns.
#
#  tabgen tab5 100000 1000
#  /usr/bin/time tabstat tab5 1:1000 median=f
#  /usr/bin/time ./tabstat.py tab5
#
#     size             C     py
#  100000 x  100      1.0   4.5
#   10000 x 1000      1.0   3.3
#  100000 x 1000     10.6  41.3

import sys
import numpy as np
from scipy.stats import skew
from scipy.stats import kurtosis

tab = sys.argv[1]

data = np.loadtxt(tab).T
(nc,nr) = data.shape

for i in range(nc):
    print(i,data[i].min(), data[i].max(), data[i].sum(),
          data[i].mean(), data[i].std(),
          skew(data[i]), kurtosis(data[i]))
