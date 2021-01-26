#! /usr/bin/env python
#
# To install:    pip install python-unsio -U
# See also:      https://pypi.org/project/python-unsio/
#                https://projets.lam.fr/projects/unsio
#
# Created on Wicked Friday (friday Sep 13, 2019) in a commit war between PJ and JC
#

import os, sys

try:
    import unsio.input as uns_in
    import unsiotools.simulations.cfalcon as falcon

    if len(sys.argv) > 1:
        myfile = sys.argv[1]
    else:
        myfile = "p100"

    if not os.path.exists(myfile):
        os.system('hackcode1 out=%s freqout=1' % myfile)
    else:
        print("Re-using existing %s" % myfile)


    my_in=uns_in.CUNS_IN(myfile,"all")

    print(my_in.getFileStructure())

    #
    # Reading Loop
    #

    while my_in.nextFrame(): 
        s0,nsel  = my_in.getData("nsel")
        s1,pos   = my_in.getData("all","pos")
        s2,timex = my_in.getData("time")
        s3,mass  = my_in.getData("all","mass")
  
        if s0 and s1 and s2:
            n1 = len(pos)
            nbody = n1//3          
            pos1=pos.reshape(nbody,3)
            print("Time=%g %s nbody=%d" % (timex,str(pos1.shape),nsel))
            print("Pos[0]: %s" % str(pos1[0]))
            # be aware that snapshots with no PhaseSpace will not be flagged as such, the same pos array is returned

    cf = falcon.CFalcon()
    ok,rho,hsml=cf.getDensity(pos,mass)
    print('pos ',pos[0],pos[1],pos[2])
    print('mass',mass[0])
    print('rho ',rho[0])
    
except:
    print("Failed loading unsio/unsiotools")
    sys.exit(1)
