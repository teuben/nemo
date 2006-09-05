#! /bin/csh -f
#
#  simple script that can be called from a crontab, some of the variables
#  here are very Maryland specific
#  Requirements:
#    - a baseline nemo.tar.gz file, created from CVS
#    - a timestamp file, created via  "date > LastBuild"
#    - this script, autobuild.csh

cd /chara2/teuben/redo_nemo/autobuild

rm -rf nemo
tar zxf nemo.tar.gz
(echo updating with cvs;cd nemo; cvs -q up -dP) >& LastBuild.cvsup.log
find nemo -type f -newer LastBuild | grep -v CVS/Entries > LastBuild.newfiles
date > LastBuild.new

if (-z LastBuild.newfiler) then
  echo No new source code in nemo since `cat LastBuild`
  exit 0
endif


echo Files updated since `cat LastBuild`
cat LastBuild.newfiles


echo Now rebuilding
tar zcf nemo.tar.gz.new nemo
cp nemo/src/scripts/test_a_new_nemo_cvs .
if (1) then
    # default gcc 3.4.1 doesn't compile falcON, so use 3.4.3 for now
    source /n/astromake2/astromake_start
    astroload gcc
endif
./test_a_new_nemo_cvs reuse=1 >& LastBuild.log

if (! -d logs) mkdir logs
foreach file (LastBuild.log autobuild.csh nemo/install.*log)
     set dest=$file:t
     cp $file logs/$dest.txt
end

 
if (0) then
  #  this will move over the build
  mv nemo.tar.gz.new nemo.tar.gz
  mv LastBuild.new LastBuild  
endif
  

  
  




