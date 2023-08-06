#! /usr/bin/env bash
#
#       here are some basic GUI elements
#
#>  IFILE   in=
#>  OFILE   out=
#>  ENTRY   eps=0.01
#>  RADIO   mode=gauss              gauss,newton,leibniz
#>  CHECK   options=mean,sigma      sum,mean,sigma,skewness,kurtosis
#>  SCALE   n=3.141592              0:10:0.001
#<          -digits 3
#<  


#		METHOD1: parse named arguments

for _a in $*;
   eval 
done

#>  -EXEC
echo "If this is executed, tkrun does not work properly yet"
#>  +EXEC

echo ARGS: in=$in out=$out eps=$eps mode=$mode options=$options n=$n

echo TESTSCRIPT
echo 0 : $0
set count=0

again:
  if ($#argv == 0) exit 0
  @ count++
  echo $count : \"$1\"
  shift argv
  goto again
   
  
