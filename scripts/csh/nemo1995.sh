#! /bin/bash
#
#   Snippet of code as was presented in 1995ASPC...77..398T
#   "The Stellar Dynamics Toolbox NEMO" , Teuben (1995)

#   Still runs, though it had a bug (the times=80 keyword should not have been used)
#   Even the URL from that 1995 paper   http://astro.umd.edu/nemo   still works.
#
#   Runs in about 1.6sec on a 2021 laptop

# clean up from a previous run
rm -rf mod1 mod2 mod5 mod5.log mod5.out

# the code example - minus the "bug"
mkplummer out=mod1 nbody=256
mkplummer mod2 1024/4
snapstack in1=mod1 in2=mod2 out=mod5 deltar=10,2,0 deltav=-0.6,0,0

hackcode1 in=mod5 out=mod5.out freqout=1 tstop=80 > mod5.log

# "times=80" was added in error to snaptrim
snaptrim in=mod5.out out=- |\
    hackforce - - |\
    snapcenter - - 'weight=-phi*phi*phi' report=f |\
    unbind - - |\
    snapplot  - xrange=0:9 yrange=-3:3 nxticks=5 nyticks=5 \
	      xvar=r 'yvar=(x*vy-y*vx)/sqrt(x*x+y*y)' \
	      times=0,10,20,30,40,50,60,70,80 nxy=3,3
