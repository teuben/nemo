#! /usr/bin/env python
#
#   Examples how to read an asci table 
#

import sys
import numpy as np
import pandas as pd
from astropy.io import ascii

if __name__ == "__main__":
    name = sys.argv[1]
    try:
        print("Trying astropy.io.ascii.read()")
        data1 = ascii.read(name)
        print(data1)
    except:
        print("Cannot use astropy")
        
    try:
        print("Trying pandas.read_table()")
        data2 = pd.read_table(name)
        print(data2)
    except:
        print("Cannot use pandas")

    try:
        print("Trying np.loadtxt()")
        data3 = np.loadtxt(name).T
        print(data3)
    except:
        print("Cannot use np.loadtxt")
        
