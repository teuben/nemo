#! /usr/bin/env python
#
# To install:    pip install python-unsio -U
# See also:      https://pypi.org/project/python-unsio/
#                https://projets.lam.fr/projects/unsio
#

import unsio.input as uns_in
import os


if True:
    os.system('hackcode1 out=p100')

myfile="p100" 

my_in=uns_in.CUNS_IN(myfile,"all")

print(my_in.getFileStructure())

#
# Reading Loop
#

while my_in.nextFrame(): 
  s1,pos   = my_in.getData("all","pos")
  s2,timex = my_in.getData("time")
  s0,nsel = my_in.getData("nsel")
  if s1 and s2 and s0:
      n1 = len(pos)
      nbody = n1//3                         #   why not something like this my_in.getData("nbody")
      pos=pos.reshape(nbody,3)
      print("Time=%g %d %d %s nbody=" % (timex,s1,s2,str(pos.shape),nsel))
      print("Pos[0]: %s" % str(pos[0]))
  # be aware that snapshots with no PhaseSpace will not be flagged as such, the same pos array is returned
