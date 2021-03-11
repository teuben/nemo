#! /bin/bash
#

# set -x


size=1000
tmp=tmp$$


for arg in $*; do
  export $arg
done


/usr/bin/time -f "CPU %e" ccdgen "" .      size=$size,$size,$size   > $tmp.log 2>&1
/usr/bin/time -f "CPU %e" ccdgen "" $tmp.1 size=$size,$size,$size  >> $tmp.log 2>&1

ts=$(grep ^CPU $tmp.log  | sed s/CPU// | tr -d '\n')
dt=$(grep ^CPU $tmp.log  | sed s/CPU// | tabtrend -  debug=-1 | awk '{print $1}')

echo "$ts " `nemoinp "($size*$size*$size)/$dt/(1000*1000)"`


rm -f $tmp.*
