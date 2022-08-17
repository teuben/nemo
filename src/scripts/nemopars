#! /usr/bin/env python
#
# in bash we can hardcode the variables we need, viz.
#
#   nemopars.sh */nemopars.rc
# with
#   for f in $*; do
#      source $f
#      echo $m $v0 $m16
#   done
#
# but we need a more generic method, e.g.
# 
#   nemopars.py m,v0,m16 */nemopars.rc
#
# hence this more complicated looking script.  Or can we do this in bash too?
#
#    11-aug-2022    toy version in bash - PJT
#    17-aug-2022    more generic python version - PJT


import sys, os
import tempfile

if len(sys.argv) < 3:
    print("Usage: %s par1,par2,...   rc1 rc2 ..." % sys.argv[0])
    print("  loops over all bash rc files that contain par1=val1 par2=val2...")
    print("  and tabulates  par1 par2 ...")
    print("  NEMO convention is that the rc files are called nemopars.rc")
    sys.exit(0)

# grab the parameters to print (no error checking)
pars = sys.argv[1].split(',')

# make a temporary name for the bash script file
s = tempfile.NamedTemporaryFile()
script = s.name
# print('SCRIPT',script)
fp = open(script,'w')
cmd = 'source $1; echo'
for p in pars:
    cmd = cmd + ' $%s' % p
fp.write(cmd + '\n')
fp.close()

# loop over the rc files
for rc in sys.argv[2:]:
    cmd = 'bash %s %s' % (script,rc)
    # print("CMD",cmd)
    os.system(cmd)
