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
#my_in=uns_in.CUNS_IN(myfile)

print(my_in.getFileStructure())

#
# Reading Loop
#
while my_in.nextFrame(): # load first snapshot
  s1,pos=my_in.getData("all","pos")
  s2,timex=my_in.getData("time")
  print("Time=%g %d %d %s" % (timex,s1,s2,str(pos.shape)))
