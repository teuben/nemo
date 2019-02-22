#! /usr/bin/env python
#
#  some fun with spherical nbody models -  30-aug-2018   Peter Teuben
#
#  old default:
#      plummer_stats.csh tmp=exp1 nexp=100 mlo=1 mhi=1
#  takes 90" on NEMO2/3, of which  ~10" in mkplummer loop
#
#  new default:
#      cluster_stats.py out='"Exp1"' nexp=100 mlo=1 mhi=1
#  takes about 3", unequal masses more like 10"
#
#  Once split: nemo loop = 1.0"   python loop = 3.3"   n=100
#                          2.5                 13.3    n=500
#                          4.8                 21.6    n=1000
#                          8.7                 43.8    n=2000
#
# Just running mkplummer nmodel=$nexp, shows the expected behavior O(nbody*nmodel)
# nexp\nbody  100     200    400
#  2000     0.185   0.425  0.737
#  4000     0.429   0.773  1.500
#  8000     0.824   1.480  3.024
# 16000     1.614   3.159  5.576
# -> 6M particles in 5.5" -> PR=1Mpps  for mkhomsph PR=1.7Mpps
# In the init loop-1 this is more like 0.23 Mpps, is that python looping, or small file I/O overhead?

import numpy as np
import numpy.ma as ma
import sys
import os
import scipy
import scipy.stats
import matplotlib.pyplot as plt


out    = 'tmcluster'             # baseline name for files
nbody  = 1000000                 # nbody for large cluster baseline plot
nsmall = 100                     # nbody for small cluster
nexp   = 2                       # number of experiments with small cluster
mlo    = 0.3                     # lower mass cutoff
mhi    = 5                       # lower mass cutoff
pimf   = -2.0                    # power low of IMF (as per massexpr=pow(m,p))
rcut   = 5.0                     # cutoff in virial radii of small cluster
tstop  = 0                       # use 10 or so if evolution is needed
rbin   = np.arange(0,4,0.125)    # radial bins
show   = 0                       # show plots on screen?
mode   = 0                       # 0=all  1=only NEMO init   2=only PYTHON analyis
quick  = 0                       # 1=just simple stats, no r_v and r_c


# poor man's command line parser --- do not change parameters below this ---
for arg in sys.argv[1:]:
    exec(arg)


# administrativia, if needed, clean up old mess 
if mode==0 or mode==1:
    os.system('rm -rf %s.*'  % out)

#@ todo make the big "nbody" sized cluster - but need a model building function for this

nbin = len(rbin)-1

vz_m  = np.zeros(nexp*nbin).reshape(nexp,nbin)
v2_m  = np.zeros(nexp*nbin).reshape(nexp,nbin)
vz_s  = np.zeros(nexp*nbin).reshape(nexp,nbin)
v2_s  = np.zeros(nexp*nbin).reshape(nexp,nbin)

# Loop1: creating data
if mode < 2:
    print("=== Loop over %d experiments of %d bodies each: creating data" % (nexp,nsmall))
    for i in range(nexp):
        file1 = "%s.%d.dat" % (out,i+1)
        if mlo==mhi:
            # faster version for equal mass starts (3" vs. 10")
            cmd = "mkplummer %s %d rfrac=%g seed=%d " \
              %      (file1, nsmall,  rcut,   i+1)
        else:
            # official version
            cmd = "mkplummer %s %d rfrac=%g seed=%d massname='n(m)' masspars=p,%g massrange=%g,%g > /dev/null 2>&1" \
              %      (file1, nsmall,  rcut,   i+1,                       pimf,         mlo, mhi)
        os.system(cmd)

        if tstop > 0:
            # @todo fix this
            cmd = 'hackcode1 $tmp.$i.dat $tmp.$i.edat tstop=$tstop > $tmp.$i.edat.log'
            cmd = 'rm $tmp.$i.dat '
            cmd = 'snaptrim $tmp.$i.edat $tmp.$i.dat  times=$tstop'

if mode==1:
    sys.exit(0)

# Loop2: analysis
print("=== Loop over %d experiments of %d bodies each: analysis" % (nexp,nsmall))
for i in range(nexp):
    file1 = "%s.%d.dat" % (out,i+1)
    file2 = "%s.%d.tab" % (out,i+1)
    cmd = 'snapprint %s r2,m,vx,vy,vz,v2 debug=-1 > %s' % (file1,file2)
    os.system(cmd)

    (r2,m,vx,vy,vz,v2) = np.loadtxt(file2,unpack=True)

    if quick == 0:
        file3 = "%s.stats" % (out)
        cmd = 'hackforce_qp %s - debug=-1 | snapstat - all=t > %s'   % (file1,file3)
        os.system(cmd)

        file4 = "%s.rvtab" % (out)
        cmd = "grep r_v %s | awk '{print $5}' >> %s" % (file3,file4)
        os.system(cmd)

        file5 = "%s.rctab" % (out)
        cmd = "grep r_c %s | awk '{print $3}' >> %s" % (file3,file5)
        os.system(cmd)

    def std(a):         # not used, but never tested if this or the lambda is faster
        return a.std()

    # the eqv. of the htab tables
    bs1 = scipy.stats.binned_statistic(r2,vz, statistic='mean',          bins=rbin)
    bs2 = scipy.stats.binned_statistic(r2,v2, statistic='mean',          bins=rbin)
    bs3 = scipy.stats.binned_statistic(r2,vz, statistic=lambda x:x.std(),bins=rbin)
    bs4 = scipy.stats.binned_statistic(r2,v2, statistic=lambda x:x.std(),bins=rbin)
    vz_m[i,:] = bs1[0]
    v2_m[i,:] = bs2[0]
    vz_s[i,:] = bs3[0]
    v2_s[i,:] = bs4[0]
    #
print("=== done ===")

# since some bins may contains NaN's, mask them so later mean/std won't become NaN
vz_m = ma.masked_invalid(vz_m)    
vz_s = ma.masked_invalid(vz_s)    
v2_m = ma.masked_invalid(v2_m)    
v2_s = ma.masked_invalid(v2_s)    

# bin radii are just average of the edges
r = (rbin[1:] + rbin[:-1])/2.0

# report the core and virial radii averages over the sample
if quick == 0:
    rv = np.loadtxt(file4)
    rc = np.loadtxt(file5)
    print("r_v:",rv.mean(),"+/-",rv.std())
    print("r_c:",rc.mean(),"+/-",rc.std())

# take mean and disp over all the experiments
vz_mm = vz_m.mean(axis=0)
v2_mm = v2_m.mean(axis=0)
vz_ms = vz_m.std(axis=0)
v2_ms = v2_m.std(axis=0)
vz_sm = vz_s.mean(axis=0)
v2_sm = v2_s.mean(axis=0)
vz_ss = vz_s.std(axis=0)
v2_ss = v2_s.std(axis=0)

# plotting, this is awkward, how to better do interactive vs. batch mode

def tabplot(f,out,show,x,y,dy,xrange,yrange,xlab,ylab):
    plt.figure(f)
    plt.errorbar(x,y,yerr=dy,fmt='o')
    plt.xlim(xrange[0],xrange[1])
    plt.ylim(yrange[0],yrange[1])
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    figname = '%s.%s.%s.png' % (out,xlab,ylab)
    plt.savefig(figname)

tabplot(1,out,show,r,vz_mm,vz_ms,[0,4],[-1,1],'r2','vz_mean')
tabplot(2,out,show,r,vz_sm,vz_ss,[0,4],[0,2], 'r2','vz_std')
tabplot(3,out,show,r,v2_mm,v2_ms,[0,4],[0,2], 'r2','v2_mean')
tabplot(4,out,show,r,v2_sm,v2_ss,[0,4],[0,2], 'r2','v2_std')
if show>0:
    plt.show()

if False:
    for i in range(len(r)):
        print(r[i],v2_mm[i],v2_ms[i],vz_mm[i],vz_ms[i])

# end    
