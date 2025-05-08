#! /usr/bin/env python
#
#    test how to parse sections of a CLI
#    Using docopt, but should be no different from e.g. argparse
# For example:
#    ./test_docopt.py -i file1 -x 2 -y 1 -i file2 -x 10 -i file3 

"""Usage: rsr_fit.py [options]

Options:
   -i --input TABLE       Input table, followed by more options. No default.
   -x --xcol XCOL         X Column(s). [Default: 1]
   -y --ycol YCOL         Y Column(s). [Default: 2]
   -h --help              This help
   -d --debug             More debug output
   -v --version           Show version

Parsing ideas for plotting tables
"""

_version = "24-apr-2025"

import os
import sys
import numpy as np
from docopt import docopt
import matplotlib.pyplot as plt

#  only once to process help and version?
#av = docopt(__doc__,options_first=True, version='test_docopt %s' % _version)


# show input
print('argv:',sys.argv)

# find the -i sections in sys.argv

argv = sys.argv[1:]
inputs = [i for i, x in enumerate(argv) if x == '-i' or x == '--input']
print('inputs:',inputs)


# loop over the input files
nf = len(inputs)
if nf == 1:
    inputs2 = [len(argv)]
else:
    inputs2 = inputs[1:] + [len(argv)]

print('inputs2:',inputs2)
        
for i in range(len(inputs)):
    print(i,argv[inputs[i]:inputs2[i]])
    sys.argv[1:] = argv[inputs[i]:inputs2[i]]
    print(sys.argv)
    av = docopt(__doc__,options_first=True, version='test_docopt %s' % _version)
    print(av)