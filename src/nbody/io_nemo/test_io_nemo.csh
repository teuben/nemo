#!/bin/csh -f


set selt=( "200:300" "100,3000:15000:143,23467:29999:466"  "all" "0:29999" "0:100:2,1000:4000:3,10000,15000:25000"  )

# check if test snapshot exist
if (! -f plum.30k ) then
    printf "Building test snapshot, wait....\n"
    mkplummer - 30000 | hackforce - plum.30k options=mass,phase,acc,phi >&! /dev/null
endif


# run tes
foreach i ( $selt )

  /bin/rm plum.30k.res plum.30k.mask >& /dev/null

  printf "Testing with seltp = [$i]\n"
  (io_nemo_test plum.30k - select=$i| snapprint - options=m,x,y,z,vx,vy,vz,ax,ay,az,phi > plum.30k.res) >& /dev/null

  (snapmask_d plum.30k - select=$i | snapprint - options=m,x,y,z,vx,vy,vz,ax,ay,az,phi >plum.30k.mask) >& ! /dev/null




  @ bad=`diff plum.30k.res plum.30k.mask | wc -l`
  if ( $bad != 0 ) then
     printf "*ERROR* plum.30k.mask and plum.30k.res differ....\n"
  endif

end

