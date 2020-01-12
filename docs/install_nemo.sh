#! /bin/bash
#
#  new V4.1+ install
#
#  the old (csh) script "nemo_install" is still available to guide you for some problems
#  but this is the recommended more cleaner procedure.   You can also use this script
#  to shortcut network installs in clone off a local nemo.git tree and use caching of
#  optional tar files.
#
#  url=file:///home/teuben/NEMO/redo/nemo.git
#
#  opt=1          hdf4 hdf5 cfitsio fftw wcslib gsl
#  python=1       uses an anaconda3 install within $NEMO
#
#  In full version (opt=1 python=1) this script takes about 13 mins
#  A simple version without falcON and any benchmarks takes about 2 mins

echo "install_nemo.sh -- Version 11-jan-2020"

opt=1
nemo=nemo
python=0
url=https://github.com/teuben/nemo

# local git repo via nemo.git for faster testing
if [ -d nemo.git ]; then
  (cd nemo.git; git pull)
  nemo=nemo42
  url=file://`pwd`/nemo.git
else
  echo "Using a direct git clone. You can can activate local hacking with:"
  echo "      git clone $url nemo.git"
fi  

for arg in $*; do\
   export $arg
done

echo "Using: "
echo "  opt=$opt"
echo "  nemo=$nemo"
echo "  python=$python"
echo "  url=$url"
echo ""

# safety
echo Sleep 5 seconds before removing nemo=$nemo ....;     sleep 5


date0=$(date)

rm -rf $nemo
git clone $url $nemo
cd $nemo


#                        when opt=1 it will first compile these packages in $NEMO/opt and use them
#                        it needs a bootstrap configure/build1
if test $opt = 1; then
    ./configure 
    make build1 MKNEMOS="hdf4 hdf5 cfitsio fftw wcslib gsl"
    opt="--with-opt"
else
    echo Skipping optional installs
    opt=""
fi

if test $python = 1; then
    make python
    source python_start.sh
else
    echo No python install
fi    

# pick a configure

#./configure $opt
./configure $opt --with-yapp=ps
#./configure $opt --with-yapp=pgplot
#./configure $opt --enable-debug --with-yapp=pgplot --with-pgplot-prefix=/usr/lib     # ok
#./configure $opt --enable-debug --with-yapp=pgplot --with-pgplot-prefix=/usr/lib   --enable-pedantic
#./configure $opt --enable-debug --with-yapp=pgplot
#./configure $opt --enable-debug --with-yapp=pgplot --with-pgplot-prefix=/usr/lib    --enable-single    # one failure nbody/init
#./configure $opt --with-yapp=pgplot --with-pgplot-prefix=/usr/lib 
#./configure $opt --enable-native
#./configure $opt --disable-falcon
#./configure $opt --disable-falcon --with-yapp=ps
#./configure $opt --with-openmp --enable-native
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
