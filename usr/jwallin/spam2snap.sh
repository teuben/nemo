#! /usr/bin/env bash
#
#   converts a series of SPAM files
#   (as listed in the files.txt file in the run directory)
#   to a single NEMO snapshot.
#
#   note: Replaces  NaN with 0.0
#    

dir=$1
out=spam.snap

if [ ! -e $dir/files.txt ]; then
    echo No files.txt in dir=$dir
    exit 0
fi

cd $dir

first=$(head -1 files.txt)
nbody=$(cat $first | wc -l)
echo Working with nbody=$nbody

rm -f $out
touch $out

# convert to snapshot (the benchmark files had 3999 particles)
# note in my 2022 version one particle was reported with NaN
# so we replace it with 0.0

t=0
while read next; do
    echo $next
    sed s/NaN/0.0/g $next | tabtos - - block1=pos,vel nbody=3999 times=$t debug=-1 >> $out
    t=$(nemoinp $t+1)
done < files.txt

echo Final output in $dir/$out
echo This current version does not have the time properly encoded
