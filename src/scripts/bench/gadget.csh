#! /bin/csh -f
#
#     example benchmark for gadget in 'galaxy' mode
#
#          cpu                                  np (number of cpu's)           comments
#                                                1      2      4
# 2 dual-Intel(R) Xeon(TM) CPU 2.66GHz:		1:44   1:02   0:58           UMD::mixcoatl 2006

set manybody=0
set np=1



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
echo time mpirun -np $np Gadget2.galaxy gadget.galaxy.param 
echo check `pwd`/gadget.galaxy.log for progress
time mpirun -np $np Gadget2.galaxy gadget.galaxy.param >& gadget.galaxy.log
lamwipe

