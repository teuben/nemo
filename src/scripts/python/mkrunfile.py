#! /usr/bin/env python3
#
#

import sys
import numpy as np    

#   simple usage
if len(sys.argv) < 3:
    print("Usage: %s [options[ progname p1=v11,v12... p2=..." % sys.argv[0])
    print("[options]   are not parsed yet")
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


# case 1: -l rundir,run_%04d,1        first case
rundir = 'rundir'
run    = 'run_%04d'
i0     = 1

# case 2: -l rundir,run_              second case
# case 3: -l rundir,run/              this case


# loop over pars and vals and create the cmds for the runfile
# @todo  generalize this

if npar == 1:
    for p0 in range(len(val[0])):
        _rdc = run % i0
        _l = '%s=%s' % (rundir, _rdc)
        i0 = i0 + 1
        cmd = '%s %s %s=%s' % (progname, _l, par[0], val[0][p0])
        print(cmd)
elif npar == 2:
    for p1 in range(len(val[1])):
        for p0 in range(len(val[0])):
            _rdc = run % i0
            _l = '%s=%s' % (rundir, _rdc)
            i0 = i0 + 1 
            cmd = '%s %s %s=%s %s=%s' % (progname, _l, par[0], val[0][p0], par[1], val[1][p1])
            print(cmd)
elif npar == 3:
    for p2 in range(len(val[2])):
        for p1 in range(len(val[1])):
            for p0 in range(len(val[0])):
                _rdc = run % i0
                _l = '%s=%s' % (rundir, _rdc)
                i0 = i0 + 1
                cmd = '%s %s %s=%s %s=%s %s=%s' % (progname, _l, par[0], val[0][p0], par[1], val[1][p1], par[2], val[2][p2])
                print(cmd)
elif npar == 4:
    for p3 in range(len(val[3])):
        for p2 in range(len(val[2])):
            for p1 in range(len(val[1])):
                for p0 in range(len(val[0])):
                    _rdc = run % i0
                    _l = '%s=%s' % (rundir, _rdc)
                    i0 = i0 + 1
                    cmd = '%s %s %s=%s %s=%s %s=%s %s=%s' % (progname, _l, par[0], val[0][p0], par[1], val[1][p1], par[2], val[2][p2], par[3], val[3][p3])
                    print(cmd)
else:
    print("# too many parameters")
                                                          


_future = """    
mkrunfile.py progname a=1,2 b=2,4 c=10
gives:
 progname a=1 b=2 c=10
 progname a=2 b=2 c=10
 progname a=1 b=4 c=10
 progname a=2 b=4 c=10

Of course this lacks an important feature: what if each line needs another parameter (e.g. rundir=) in
which the parameter values, or an ordinal, is encoded. E.g.

 progname a=1 b=2 c=10 rundir=run010
 progname a=2 b=2 c=10 rundir=run011
 progname a=1 b=4 c=10 rundir=run012
 progname a=2 b=4 c=10 rundir=run013

where the base number (run010) can be choosen as well

or
 progname a=1 b=2 c=10 rundir=run_1_2_10
 progname a=2 b=2 c=10 rundir=run_2_2_10
 progname a=1 b=4 c=10 rundir=run_1_4_10
 progname a=2 b=4 c=10 rundir=run_2_4_10

or even
 progname a=1 b=2 c=10 rundir=runs/1/2/10
 progname a=2 b=2 c=10 rundir=runs/2/2/10
 progname a=1 b=4 c=10 rundir=runs/1/4/10
 progname a=2 b=4 c=10 rundir=runs/2/4/10

so perhaps:

 -l rundir,run%03d,10        first case
 -l rundir,run_              second case
 -l rundir,run/              this case
"""
