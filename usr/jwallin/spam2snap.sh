#! /usr/bin/env bash
#
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

# convert to snapshot
# note the 

while read next; do
    echo $next
    sed s/NaN/0.0/g $next | tabtos - - block1=pos,vel nbody=3999 debug=-1 >> $out 
done < files.txt

echo Final output in $dir/$out
echo This current version does not have the time properly encoded
