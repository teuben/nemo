#! /usr/bin/env python
#
# https://docs.python.org/3/library/argparse.html

import math
import argparse

#             somehow formatting with \n doesn't work
my_help = """
  This script grabs aaa and bbb via argparse. \nThere are several options:\n
  --aaa 10\n
  --aaa=10\n
  -a 10      (note the = option is not allowed in short form)
"""

p = argparse.ArgumentParser(description=my_help)
p.add_argument('--aaa', type = int,   default = 1,   help='The aaa (int) variable')
p.add_argument('--bbb', type = float, default = 2.0, help='The bbb (float) variable')

args = p.parse_args()

if __name__ == '__main__':

    print(args.aaa, args.bbb)

