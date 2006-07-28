#! /bin/csh -f
#
#     example benchmark for gadget in 'galaxy' mode
#

set manybody=0




foreach a ($*)
  set $a
end

if ($manybody == 0) then
  if ($?MANYBODY == 0) then
    echo Cannot run benchmark without MANYBODY
    exit 1
  endif
  set manybody=$MANYBODY
  if (! -e $manybody/gadget/Gadget/ICs) then
    echo $manybody/gadget/Gadget/ICs does not exist
    exit 1
  endif
endif


set tmp=tmp

rm -rf $tmp
mkdir $tmp
cp gadget.galaxy.param $tmp
cd $tmp

ln -s $manybody/gadget/Gadget/ICs
mkdir galaxy

lamboot
echo time mpirun -np 2 Gadget2.galaxy gadget.galaxy.param 
echo check `pwd`/gadget.galaxy.log for progress
time mpirun -np 2 Gadget2.galaxy gadget.galaxy.param >& gadget.galaxy.log
lamwipe

