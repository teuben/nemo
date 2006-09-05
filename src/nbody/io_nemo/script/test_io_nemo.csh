#!/bin/csh -f
#
#

# Check necessary nemo binaries
foreach i ( hackforce snapmask snapprint mkplummer )
    set x=""`$NEMOSRC/scripts/need $i`
    echo "$i [$x]"
    if ( $x != "" ) then
	$NEMOSRC/scripts/need -m  $i
    endif
end
rehash

set RT=/tmp/.io_nemo_run_test

# Create 'run_test' directory
if ( ! -d ${RT} ) then
    mkdir -p ${RT}
endif

# Set io_nemo test programs path
set IONB=compile/bin

#set selt=( "200:300" "100,3000:15000:143,23467:29999:466"  "all" "0:29999" "0:100:2,1000:4000:3,10000,15000:25000"  )
set selt=( "0:100:2,1000:4000:3,10000,15000:25000"  )
set prog=( float double )
set mask=( _s _d )

# check if test snapshot exist
if (! -f ${RT}/plum.30k ) then
    printf "Building test snapshot, wait....\n"
    mkplummer - 30000 | hackforce - - options=mass,phase,acc,phi|snapmask - ${RT}/plum.30k  >&! /dev/null
endif

@ mm = 1
foreach j ( $prog )
    # compute modulo for SNAPMASK extension's name
    @ mm++
    @ snap = $mm % 2
    @ snap++

    set run=${IONB}/io_nemo_test_${j}
    set snpmsk=${IONB}/snapmask${mask[$snap]} #snapmask program

    printf "\n*** Running [$run] with [$snpmsk] ***\n"

# run tes
    foreach i ( $selt )

	/bin/rm ${RT}/plum.30k.* ${RT}/key.* >& /dev/null

	printf "Testing with seltp = [$i]\n"
	( ${run} ${RT}/plum.30k ${RT}/plum.30k.ionemo select=$i ) >& /dev/null
	( snapprint ${RT}/plum.30k.ionemo  options=m,x,y,z,vx,vy,vz,ax,ay,az,phi,key > ${RT}/plum.30k.res ) >& /dev/null
	($snpmsk   ${RT}/plum.30k - select=$i| snapprint - options=m,x,y,z,vx,vy,vz,ax,ay,az,phi,key > ${RT}/plum.30k.mask) >& ! /dev/null

	@ bad=`diff ${RT}/plum.30k.res ${RT}/plum.30k.mask | wc -l`
	if ( $bad != 0 || -z ${RT}/plum.30k.res || -z ${RT}/plum.30k.mask ) then
	    printf "*ERROR* ${RT}/plum.30k.mask and ${RT}/plum.30k.res differ....\n"
	    exit(1)
	endif
    end
end
#
