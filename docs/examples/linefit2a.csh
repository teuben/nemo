#! /bin/csh -f
#

set n=1

foreach arg ($*)
  set $arg
end

#

foreach i (`nemoinp 1:$n`)
  linefit2.csh $* doplot=0 >& linefit2.log
  set a=(`grep 'OLS(Y/X)' linefit2.log | awk '{print $4}'`)
  echo $i $a
end
