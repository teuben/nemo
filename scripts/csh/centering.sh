#! /usr/bin/env bash
#
#      experiments different ways to center a plummer sphere
#  ./centering.sh n=10000 > c10000.tab
#  tabstat c10000.tab 1,2,3,7,8,9,13,14,15,19,20,21
#  tabstat c10000.tab 4,5,6,10,11,12,16,17,18,22,23,24


nbody=100
seed=-2
n=10

for _a in $*; do
    export $_a
done

export DEBUG=-1

for i in $(seq $n); do
   mkplummer - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r > p.dat
   c1=$(snapcenter p.dat . report=t)
   c2=$(hackdens p.dat - | snapcenter  - . weight=dens report=t)
   c3=$(hackforce p.dat - fcells=3 | snapcenter - . weight="-phi*phi*phi" report=t)
   c4=$(hackforce p.dat - fcells=3 | snapcenter - . weight="phi*phi*phi*phi" report=t)


   echo $c1 $c2 $c3 $c4
done
