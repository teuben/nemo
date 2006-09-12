#! /bin/csh -f
#
#  simple script that can be called from a crontab, some of the variables
#  here are very Maryland specific
#  Requirements:
#    - a baseline nemo.tar.gz file, created from CVS
#    - a timestamp file, created via  "date > LastBuild"
#    - this script, autobuild.csh
#
# $Id$

set sleep=0
set dir=/chara2/teuben/redo_nemo/autobuild

foreach arg ($*)
  set $arg
end

if (! -e $dir) then
  echo Directory dir=$dir does not exist
  exit 1
endif


start:

cd $dir

rm -rf nemo
tar zxf nemo.tar.gz
(echo updating with cvs;cd nemo; cvs -q up -dP) >& LastBuild.cvsup.log
find nemo -type f -newer LastBuild | grep -v CVS/Entries > LastBuild.newfiles
date > LastBuild.new

if (-z LastBuild.newfiles) then
# echo No new source code in nemo since `cat LastBuild`
# echo But building anyways....
  goto last_check
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

grep TESTSUITE: LastBuild.log | grep -v OK > LastBuild.fail
foreach file (LastBuild.log LastBuild.fail LastBuild.newfiles autobuild.csh nemo/install.*log)
     set dest=$file:t
     cp $file logs/$dest.txt
end

if (! -z LastBuild.fail) cat LastBuild.fail

 
if (1) then
  #  this will move over the build
  mv nemo.tar.gz.new nemo.tar.gz
  cp -p LastBuild.new LastBuild  
endif


last_check:


cmp -s autobuild.csh nemo/src/scripts/autobuild.csh
if ($status) then
  echo Warning: local autobuild.csh differs from nemo/src/scripts/autobuild.csh
  echo dir=$dir
endif


if ($sleep) then
  sleep $sleep
  goto start
endif
