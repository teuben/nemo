#! /bin/bash
#
#  new V4.1 install, very few bells and whistles, just commented lines
#  to suggest another approach
#
#  the old (csh) script "nemo_install" is still available to guide you for some problems

opt=1
nemo=nemo

date0=$(date)

rm -rf $nemo
git clone https://github.com/teuben/nemo $nemo
cd $nemo


#                        when opt=1 it will first compile these packages in $NEMO/opt and use them
if test $opt = 1; then
  ./configure
  make hdf4
  make hdf5
  make cfitsio
  ./configure --with-opt
fi  

# pick a configure

#./configure                       # ok, this will pick yapp_ps
#./configure --with-yapp=pgplot --disable-falcon
#./configure --enable-debug --with-yapp=pgplot --with-pgplot-prefix=/usr/lib     # ok
#./configure --enable-debug --with-yapp=pgplot --with-pgplot-prefix=/usr/lib   --enable-pedantic
#./configure --enable-debug --with-yapp=pgplot
#./configure --enable-debug --with-yapp=pgplot --with-pgplot-prefix=/usr/lib    --enable-single    # one failure nbody/init
#./configure  --with-yapp=pgplot --with-pgplot-prefix=/usr/lib 
./configure --enable-native
#./configure --disable-falcon
#./configure --enable-cfitsio --with-cfitsio-prefix=/usr
#./configure --with-yapp=plplot --disable-falcon
#./configure --with-opt


# on a mac:
#./configure --without-csh 


# and/or re-configure
#./config.status --recheck

#            make build (=build1,2,3,4)
make build1
make build2
make build3
make build4

make check
make bench

date1=$(date)

echo "Started: $date0"
echo "Ended:   $date1"
