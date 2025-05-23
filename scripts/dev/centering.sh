#! /usr/bin/env bash
#
#      experiments different ways to center a plummer sphere
#      Mentioned in the man page for hackdens
#
#  ./centering.sh n=10000 > c10000.tab
#  tabstat c10000.tab 1,2,3,7,8,9,13,14,15,19,20,21
#  tabstat c10000.tab 4,5,6,10,11,12,16,17,18,22,23,24

#  @todo   add half mass radius:   mkplummer - nbody=$nbody seed=$seed zerocm=f | snapmradii - 0.5

nbody=100
seed=-3
n=10
shift=0
nfrac=3
model=0     # 0=plummer    via mkommod:  1,3,5 king   -1 plummer -2 devauc     6  dehnen

for _a in $*; do
    export $_a
done

export DEBUG=-1

for i in $(seq $n); do
    if [ $model = 0 ]; then
	mkplummer                          - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r | snapshift - - $shift,0,0 > p.dat
    elif [ $model = -1 ]; then
	mkommod $NEMODAT/plum.dat          - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r | snapshift - - $shift,0,0 > p.dat
    elif [ $model = -2 ]; then
	mkommod $NEMODAT/devauc.dat        - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r | snapshift - - $shift,0,0 > p.dat	
    elif [ $model = 1 ]; then
	mkommod $NEMODAT/k${model}isot.dat - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r | snapshift - - $shift,0,0 > p.dat
    elif [ $model = 3 ]; then
	mkommod $NEMODAT/k${model}isot.dat - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r | snapshift - - $shift,0,0 > p.dat
    elif [ $model = 5 ]; then
	mkommod $NEMODAT/k${model}isot.dat - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r | snapshift - - $shift,0,0 > p.dat
    elif [ $model = 6 ]; then
	mkdehnen - nbody=$nbody seed=$seed gamma=1 | snapsort - - rank=r | snapshift - - $shift,0,0 > p.dat	
    else
	echo "Not a model=$model"
	exit 0
    fi
   c1=$(snapcenter p.dat . report=t)
   c2=$(hackdens p.dat - | snapcenter  - . weight=dens report=t)
   c3=$(hackforce p.dat - fcells=3 | snapcenter - . weight="-phi*phi*phi" report=t)
   #c4=$(hackforce p.dat - fcells=3 | snapcenter - . weight="phi*phi*phi*phi" report=t)
   #c5=$(hackforce p.dat - fcells=3 | snapmnmx - phi min out=- | snapprint -)
   c6=$(hackforce p.dat - fcells=3 | snapsort - - phi | snapcopy - - "i<$nbody/$nfrac" | snapcenter - . report=t)
   #c7=$(hackforce p.dat - fcells=3 | snapsort - - phi | snapcopy - - "i<$nbody/$nfrac" | snapcenter - . weight="-phi*phi*phi" report=t)
   #c8=$(hackforce p.dat - fcells=3 | snapcenter - . weight="phi*phi" report=t)
   #c9=$(hackforce p.dat - fcells=3 | snapcenter - . weight="-phi" report=t)
   echo $c1 $c2 $c3 $c6
   
   #  snapcenterp

   #  c5: center based on the particle with the deepest potential; use snapmnmx | snapprint
   #  c6: center based on some fraction of particles sorted by potential; use snapsort | snapcopy i < nbody/k or so
   #  c7: center based on some fraction of particles sorted by potential and weighted by phi*3
   
   #echo $c1 $c2 $c3 $c4
done
