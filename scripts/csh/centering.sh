#! /usr/bin/env bash
#
#      experiments different ways to center a plummer sphere
#      Mentioned in the man page for hackdens
#
#  ./centering.sh nbody=100 n=10000 > c10000.tab
#  tabstat c10000.tab 1:3,7:9,13:15,19:21
#  tabstat c10000.tab 4:6,10:12,16:18,22:24


#  @todo   add half mass radius:   mkplummer - nbody=$nbody seed=$seed zerocm=f | snapmradii - 0.5

nbody=100
seed=-3
n=10
nfrac=3
model=0     # 0=plummer    via mkommod:  1,3,5 king   -1 plummer -2 devauc     6  dehnen
# if gridding, npix > 0 is needed
npix=0
size=2
box=8

for _arg in "$@"; do
    export "$_arg"
done


export DEBUG=-1

for i in $(seq $n); do
    # create a radius sorted model in p.dat
    
    if [ $model = 0 ]; then
	mkplummer                          - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r > p.dat
    elif [ $model = -1 ]; then
	mkommod $NEMODAT/plum.dat          - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r > p.dat
    elif [ $model = -2 ]; then
	mkommod $NEMODAT/devauc.dat        - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r > p.dat	
    elif [ $model = 1 ]; then
	mkommod $NEMODAT/k${model}isot.dat - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r > p.dat
    elif [ $model = 3 ]; then
	mkommod $NEMODAT/k${model}isot.dat - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r > p.dat
    elif [ $model = 5 ]; then
	mkommod $NEMODAT/k${model}isot.dat - nbody=$nbody seed=$seed zerocm=f | snapsort - - rank=r > p.dat
    elif [ $model = 6 ]; then
	mkdehnen                           - nbody=$nbody seed=$seed gamma=1  | snapsort - - rank=r > p.dat	
    else
	echo "Not a model=$model"
	exit 0
    fi

    # output various centroids

    if [ $npix = 0 ]; then
    
	c1=$(snapcenter p.dat . report=t)
	c2=$(hackdens p.dat - | snapcenter  - . weight=dens report=t)
	c3=$(hackforce p.dat - fcells=3 | snapcenter - . weight="-phi*phi*phi" report=t)
	#c4=$(hackforce p.dat - fcells=3 | snapcenter - . weight="phi*phi*phi*phi" report=t)
	#c5=$(hackforce p.dat - fcells=3 | snapmnmx - phi min out=- | snapprint -)
	c6=$(hackforce p.dat - fcells=3 | snapsort - - phi | snapcopy - - "i<$nbody/$nfrac" | snapcenter - . report=t)
	#c7=$(hackforce p.dat - fcells=3 | snapsort - - phi | snapcopy - - "i<$nbody/$nfrac" | snapcenter - . weight="-phi*phi*phi" report=t)
	#c8=$(hackforce p.dat - fcells=3 | snapcenter - . weight="phi*phi" report=t)
	#c9=$(hackforce p.dat - fcells=3 | snapcenter - . weight="-phi" report=t)

	# each center has 6 values, 3 position and 3 velocity
	echo "$c1  $c2  $c3  $c6"
   
	#  snapcenterp

	#  c5: center based on the particle with the deepest potential; use snapmnmx | snapprint
	#  c6: center based on some fraction of particles sorted by potential; use snapsort | snapcopy i < nbody/k or so
	#  c7: center based on some fraction of particles sorted by potential and weighted by phi*3
	
	#echo $c1 $c2 $c3 $c4

    else
	rm -rf p.ccd
	snapgrid p.dat p.ccd nx=$npix ny=$npix xrange=-${size}:${size} yrange=-${size}:${size}
	ccdblob p.ccd pos=$npix/2,$npix/2 wcs=f box=$box
	nds9 p.ccd
    fi
done
