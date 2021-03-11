#! /bin/bash
#

# set -x


size=1000
tmp=tmp$$
debug=0

for arg in $*; do
  export $arg
done

if [ $debug = 1 ]; then
    set -x
fi

/usr/bin/time -f "CPU %e %U %S" ccdgen "" .      size=$size,$size,$size   > $tmp.log 2>&1
/usr/bin/time -f "CPU %e %U %S" ccdgen "" $tmp.1 size=$size,$size,$size  >> $tmp.log 2>&1

if [ $debug = 1 ]; then
    cat $tmp.log
fi

ts=$(grep ^CPU $tmp.log  | sed s/CPU// | tr -d '\n')
dt=$(grep ^CPU $tmp.log  | sed s/CPU// | tabtrend -  xcol=1 debug=-1 | awk '{print $1}')

echo "$ts " `nemoinp "($size*$size*$size)/$dt/(1000*1000)"`


rm -f $tmp.*
