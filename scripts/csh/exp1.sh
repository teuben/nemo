#! /bin/bash
#
#  Some experiments on slightly cooled models
#  1)   Are oscillations in Plummer, OM-models, Polytrope different
#       after the initial contraction
#
#  Examples to run and plot:
#
# ./exp1.sh nbody=10000 tstop=10 run=run4
# ./exp1.sh nbody=10000 tstop=10 run=run5 vscale=1
# ./exp1.sh nbody=10000 tstop=10 run=run6 vscale=0.75
#
# ./exp1.py  run4/p01.cgs/snap.mrad run5/p01.cgs/snap.mrad
# ./exp1.py  run4/p01.cgs/snap.mrad run5/p01.cgs/snap.mrad run6/p01.cgs/snap.mrad
# ./exp1.py  run4/p02.cgs/snap.mrad run5/p02.cgs/snap.mrad run6/p02.cgs/snap.mrad
# ./exp1.py  run4/p04.cgs/snap.mrad run5/p04.cgs/snap.mrad run6/p04.cgs/snap.mrad
# ./exp1.py  run4/p31.cgs/snap.mrad run5/p31.cgs/snap.mrad run6/p31.cgs/snap.mrad
# ./exp1.py  run4/p11.cgs/snap.mrad run5/p11.cgs/snap.mrad run6/p11.cgs/snap.mrad
# ./exp1.py  run4/p13.cgs/snap.mrad run5/p13.cgs/snap.mrad run6/p13.cgs/snap.mrad

#
#  2-apr-2021   quick and dirty, with  n=1e4 and 6 models
#  3-apr-2021   generalized for more models
#

nbody=10000
run=run1
vscale=0.5
g=4.4971       #  check why this isn't 4.30071
seed=123
tstop=10
clean=1

#  simple keyword=value command line parser for bash - don't make any changing below
for arg in $*; do
  export $arg
done


if [ -d $run ]; then
    echo Run directory $run already exists
    exit 1
fi

mkdir -p $run
cd $run

echo nbody=$nbody run=$run vscale=$vscale g=$g seed=$seed tstop=$tstop > exp1.rc

# directories which hold initial, scaled and DF functions
mkdir i j k

# in virial:
# p02 oscillates
# p06 is nuts
# p31,33,35,37 seem to oscillate a bit, especially the earlier ones

anisot k/plum.dat plummer rmax=10
anisot k/k1.dat   king    rmax=10 w0=1 
anisot k/k3.dat   king    rmax=10 w0=3 
anisot k/k5.dat   king    rmax=10 w0=5 
anisot k/k7.dat   king    rmax=10 w0=7

mkplummer                    i/p01 $nbody seed=$seed
mkpolytrope                  i/p02 $nbody seed=$seed
mkplum                       i/p03 $nbody seed=$seed
mkdehnen                     i/p04 $nbody seed=$seed gamma=0
mkdehnen                     i/p05 $nbody seed=$seed gamma=1
mkdehnen                     i/p06 $nbody seed=$seed gamma=2

mkommod k/plum.dat           i/p21 $nbody seed=$seed

#mkommod k/k1.dat            i/p22 $nbody seed=$seed
#mkommod k/k3.dat            i/p23 $nbody seed=$seed
#mkommod k/k5.dat            i/p24 $nbody seed=$seed
#mkommod k/k7.dat            i/p25 $nbody seed=$seed

mkking                       i/p31 $nbody seed=$seed W0=1
mkking                       i/p33 $nbody seed=$seed W0=3
mkking                       i/p35 $nbody seed=$seed W0=5
mkking                       i/p37 $nbody seed=$seed W0=7 

mkommod $NEMODAT/plum.dat    i/p41 $nbody seed=$seed
mkommod $NEMODAT/devauc.dat  i/p42 $nbody seed=$seed

mkommod $NEMODAT/k1isot.dat  i/p11 $nbody seed=$seed
mkommod $NEMODAT/k3isot.dat  i/p13 $nbody seed=$seed
mkommod $NEMODAT/k5isot.dat  i/p15 $nbody seed=$seed
mkommod $NEMODAT/k7isot.dat  i/p17 $nbody seed=$seed

models=$(cd i; ls)

for p in $models; do
    echo ==== $p =====
    snapscale i/$p j/$p vscale="sqrt($g)*$vscale"
    runCGS out=$p.cgs in=j/$p tstop=$tstop maxstep=999999999 dt=0.001 ndtout=5 ndtcmss=400
    snapcopy $p.cgs/snap.out - | snapcenter - - |\
	snapmradii - 0.01,0.1:0.9:0.1,0.99 > $p.cgs/snap.mrad
    snapmradii $p.cgs/snap.out 0.01,0.1:0.9:0.1,0.99 > $p.cgs/snap.mrad0
    # remove big files ?
    if [ $clean = 1 ]; then
	rm -f $p.cgs/snap.out
    fi
done



  
