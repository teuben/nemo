#! /bin/csh -f
#


set a=2

set b1=1
set c1=4
set d1=1

set b2=2
set c2=6
set d2=0.2

set x=0:10:0.1
set sig=0.1

#set par=2.1,1.1,4.1,1.1,2.1,6.1,0.11
set par=2,1,4,1,2,6,0.2
set free=1,1,1,1,1,1,1

foreach _arg ($*)
  set $_arg
end

nemoinp $x |\
   tabmath - - "$b1*exp(-(%1-$c1)**2/(2*$d1**2))" |\
   tabmath - - "$b2*exp(-(%1-$c2)**2/(2*$d2**2))" |\
   tabmath - - "%1,%2+%3+$a+rang(0,$sig)"  all  > tab1


make -f $NEMOLIB/Makefile.lib gauss2.so


echo tabnllsqfit tab1 load=gauss2.so fit=gauss2 par=$par free=$free
tabnllsqfit tab1 load=gauss2.so fit=gauss2 par=$par free=$free
