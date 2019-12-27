#! /bin/bash
#
#  new V4.1 install, very few bells and whistles, just commented lines
#  to suggest another approach
#
#  the old (csh) script "nemo_install" is still available to guide you for some problems

opt=1
nemo=nemo
python=0

for arg in $*; do\
   export $arg
done

echo "Using: "
echo "  opt=$opt"
echo "  nemo=$nemo"
echo "  python=$python"
echo ""


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
    opt="--with-opt"
else
    opt=""
fi

if test $python = 1; then
    make python
    source python_start.sh
fi    

# pick a configure (remember to use $opt)

#./configure $opt                       # ok, this will pick yapp_ps
#./configure $opt --with-yapp=pgplot --disable-falcon
#./configure $opt --enable-debug --with-yapp=pgplot --with-pgplot-prefix=/usr/lib     # ok
#./configure $opt --enable-debug --with-yapp=pgplot --with-pgplot-prefix=/usr/lib   --enable-pedantic
#./configure $opt --enable-debug --with-yapp=pgplot
#./configure $opt --enable-debug --with-yapp=pgplot --with-pgplot-prefix=/usr/lib    --enable-single    # one failure nbody/init
#./configure $opt --with-yapp=pgplot --with-pgplot-prefix=/usr/lib 
#./configure $opt --enable-native
#./configure $opt --disable-falcon
#./configure $opt --enable-cfitsio --with-cfitsio-prefix=/usr
#./configure $opt --with-yapp=plplot --disable-falcon
#./configure $opt --with-opt


# on a mac:
#./configure $opt --without-csh 


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
