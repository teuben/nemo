#! /usr/bin/env python
#
# To install:    pip install python-unsio -U
# See also:      https://pypi.org/project/python-unsio/
#                https://projets.lam.fr/projects/unsio
#
# Created on Wicked Friday (friday Sep 13, 2019) in a commit war between PJ and JC
#

import unsio.input as uns_in
import os, sys

if len(sys.argv) > 1:
    myfile = sys.argv[1]
else:
    myfile = "p100"

if not os.path.exists(myfile):
    os.system('hackcode1 out=%s' % myfile)


my_in=uns_in.CUNS_IN(myfile,"all")

print(my_in.getFileStructure())

#
# Reading Loop
#

while my_in.nextFrame(): 
  s0,nsel  = my_in.getData("nsel")
  s1,pos   = my_in.getData("all","pos")
  s2,timex = my_in.getData("time")
  if s0 and s1 and s2:
      n1 = len(pos)
      nbody = n1//3          
      pos=pos.reshape(nbody,3)
      print("Time=%g %s nbody=%d" % (timex,str(pos.shape),nsel))
      print("Pos[0]: %s" % str(pos[0]))
  # be aware that snapshots with no PhaseSpace will not be flagged as such, the same pos array is returned
