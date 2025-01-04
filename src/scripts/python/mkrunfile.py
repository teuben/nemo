#! /usr/bin/env python3
#
#

import sys
import numpy as np    

#   simple usage
if len(sys.argv) < 3:
    print("Usage: %s progname p1=v11,v12... p2=..." % sys.argv[0])
    sys.exit(0)

#   parse command line args
progname = sys.argv[1]
pars = sys.argv[2:]

npar = len(pars)
par = list(range(npar))
val = list(range(npar))

# setting
for i in range(npar):
    par[i] = 'par%d' % i
    par[i] = pars[i].split('=')[0]
    val[i] = pars[i].split('=')[1].split(',')

# loop over pars and vals and create the cmds for the runfile
# @todo  generalize this

if npar == 1:
    for p0 in range(len(val[0])):
        cmd = '%s %s=%s' % (progname, par[0], val[0][p0])
        print(cmd)
elif npar == 2:
    for p1 in range(len(val[1])):
        for p0 in range(len(val[0])):
            cmd = '%s %s=%s %s=%s' % (progname, par[0], val[0][p0], par[1], val[1][p1])
            print(cmd)
elif npar == 3:
    for p2 in range(len(val[2])):
        for p1 in range(len(val[1])):
            for p0 in range(len(val[0])):
                cmd = '%s %s=%s %s=%s %s=%s' % (progname, par[0], val[0][p0], par[1], val[1][p1], par[2], val[2][p2])
                print(cmd)
elif npar == 4:
    for p3 in range(len(val[3])):
        for p2 in range(len(val[2])):
            for p1 in range(len(val[1])):
                for p0 in range(len(val[0])):
                    cmd = '%s %s=%s %s=%s %s=%s %s=%s' % (progname, par[0], val[0][p0], par[1], val[1][p1], par[2], val[2][p2], par[3], val[3][p3])
                    print(cmd)
else:
    print("# too many parameters")
                                                          
                                                        
