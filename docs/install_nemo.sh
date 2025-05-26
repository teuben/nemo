#! /bin/bash
#
#--HELP
#  Install script for NEMO V4.1+ install
#
#  You can also use this script to shortcut network installs by cloning
#  off a local nemo.git tree and use caching of optional tar files.
#
#  The script has a number of key=val optional keywords. Their defaults
#  are listed below
#
#  In full version (opt=1 python=1) this script takes about 14 mins
#  A simple version without falcON and any checks and benchmarks takes about 1-2 mins
#   opt=0 python=0 falcon=0                  1'50"
#   opt=0 python=0 falcon=0 mknemos=cfitsio  2'13"
#   opt=0 python=0                           5'30"
#   opt=0 python=1                           6'15"
#   opt=1 python=1                          14'00"
#

echo "install_nemo.sh:  Version 1.11 -- 25-may-2025"

 opt=0                                                # install some optional mknemos= packages
 nemo=nemo                                            # root directory where NEMO will be installed (. = here)
 branch=master                                        # branch used
 python=0                                             # 1=install anaconda3 python? 2=post-install pip
 url=https://github.com/teuben/nemo                   # git repo adres
 mknemos=hdf4,hdf5,cfitsio,fftw,wcslib,gsl,netcdf4    # optional packages via source (needs opt=1)
 falcon=1                                             # install falcON tools?
 amuse=0                                              # install amuse via python? (needs python=2)
 yapp=auto                                            # yapp driver   (auto, ps, pgplot, pglocal, plplot)
 check=1                                              # "make check" at the end?
 bench=0                                              # "make bench" at the end?
 bench5=1                                             # "make bench5" at the end?
 build=3                                              # highest build level
 v=""                                                 # append version (e.g. "-12" or "12" to gcc/clang)
#--HELP

if [ "$1" == "--help" ] || [ "$1" == "-h" ];then
    set +x
    awk 'BEGIN{s=0} {if ($1=="#--HELP") s=1-s;  else if(s) print $0; }' $0
    exit 0
fi

# indirect git clone via nemo.git for faster testing
# bootstrap with:    git clone https://github.com/teuben/nemo nemo.git
if [ -d nemo.git ]; then
  url=file://`pwd`/nemo.git
  echo Found local nemo.git - using this as $url
  (cd nemo.git; git pull)
else
  echo "Using a direct git clone. You can can activate local hacking with:"
  echo "      git clone $url nemo.git"
fi  

#             simple keyword=value command line parser for bash
for arg in "$@"; do
  export "$arg"
done

echo "Using: "
echo "  nemo=$nemo"
echo "  branch=$branch"
echo "  opt=$opt"
echo "  mknemos=$mknemos"
echo "  python=$python"
echo "  url=$url"
echo "  falcon=$falcon"
echo "  yapp=$yapp"
echo "  check=$check"
echo "  bench=$bench"
echo "  bench5=$bench5"
echo "  v=$v"
echo ""

if [ "$nemo" == "." ]; then
    # use local directory, but one has to be in the NEMO root directory
    if [ ! -e nemo_start.sh.in ]; then
	echo No nemo_start.sh.in here
	exit 0
    fi
    nemo=$(pwd)
else
    if [ -e nemo_start.sh.in ]; then
	echo "Warning: you seem to want to make a new nemo inside an old nemo...."
	sleep 5
    fi
    # safety
    if [ -d $nemo ]; then
	echo Sleep 5 seconds before removing nemo=$nemo ....
	sleep 5
    else
	echo Installing NEMO in a new directory $nemo
    fi    
    date0=$(date)
    rm -rf $nemo
    git clone $url $nemo
    cd $nemo
fi
git checkout $branch
./configure 

#                        when opt=1 it will first compile these packages in $NEMO/opt and use them
#                        it needs a bootstrap configure/build1
if test $opt = 1; then
    mknemos=`echo $mknemos | sed 's/,/ /g'`
    make build1 MKNEMOS="$mknemos"
    opt="--with-opt"
else
    echo "Skipping optional installs (opt=0)"
    opt=""
fi

if test $falcon = 0; then
    opt="$opt --disable-falcon"
fi


if [ $python -gt 0 ]; then
    make build1
    NEMO=`pwd` make python
    source anaconda3/python_start.sh
    source nemo_start.sh
    make nemopy
    if [ $python  -eq 1 ]; then
	mknemo amuse
	(cd $NEMO/local/amuse; pip install -e .; make bhtree.code)
    elif [ $python  -eq 2 ]; then
	src/scripts/install_python_modules amuse=$amuse
    fi
else
    echo No python install
fi

if [ $yapp = "auto" ]; then
    with_yapp=""
elif [ $yapp = "pglocal" ]; then    
    source nemo_start.sh
    with_yapp="--with-yapp=pgplot --enable-png --with-pgplot-prefix=$NEMOLIB"
    src/scripts/pgplot.install png=1
else
    with_yapp="--with-yapp=$yapp"
fi
   

# pick a configure

CC=gcc$v CXX=g++$v F77=gfortran$v ./configure $opt $with_yapp
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
#./configure $opt 


# on a mac:
#./configure $opt --without-csh 


# and/or re-configure
#./config.status --recheck

#            make build (=build1,2,3,4)
make build1
[ $build -gt 1 ] && make build2
[ $build -gt 2 ] && make build3
[ $build -gt 3 ] && make build4
[ $python  = 1 ] && make testpython
[ $check   = 1 ] && make check
[ $bench   = 1 ] && make bench
[ $bench5  = 1 ] && make bench5


date1=$(date)

echo "Started: $date0"
echo "Ended:   $date1"

echo All done installing $nemo
echo ""
fnemo=$(pwd)
echo "(ba)sh users:  source $fnemo/nemo_start.sh  to activate NEMO in your shell"
echo "(t)csh users:  source $fnemo/nemo_start.csh to activate NEMO in your shell"

