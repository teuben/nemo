#! /usr/bin/env python
#
#   Example how to read an asci table into astropy or pandas

from astropy.io import ascii
import pandas as pd
import sys

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
