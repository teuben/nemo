#! /usr/bin/env python
#
#    test how to parse sections of a CLI using argparse
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

import sys, argparse

print('argv:', sys.argv)

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

# Defining parser
# add_argument attaches individual assignment specifications to the parser

# Add help
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True) # -i/--input argument
parser.add_argument('-x', '--xcol', default=1, type=int)
parser.add_argument('-y', '--ycol', default=2, type=int)

# Store parsed results
avs = []
for i in range(len(inputs)):
    section_args = argv[inputs[i]:inputs2[i]]
    print(f"\nGroup {i}: {section_args}")
    av = parser.parse_args(section_args)
    print(f"Parsed av[{i}]:", av)
    avs.append(av)


#print(avs)

