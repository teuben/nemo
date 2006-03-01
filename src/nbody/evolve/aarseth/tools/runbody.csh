#! /bin/csh -f
#
#  This is an example C-shell script that launches a series of
#  nbodyX to search in a two dimensional parameter space of 
#  a particular problem, in this case Q and ALPHA is varied.
#  Each run is saved in a subdirectory, where the subdirectory
#  encodes these parameters. It is also handy to save them in a
#  file in the run directory, which we call nbody1.par here.
#
#  This particular example does a 5x5 grid on 25 particles for
#  10 crossing times, and takes about 20 secs on my current 
#  laptop (P1.6 GHz)
#
#  After you have run this example, here are some commands that
#  illustrate what you can do
#         grep ok series105/*/nbody1.par
#  At the end of this script you can also see an example shell
#  script style analysis
#
# History:
#  1-mar-2006     Created after an excellent Eagle Pub dinner with Sverre        PJT
#
#

set master=series105
set a_values=(1.0 1.5 2.0 2.5 3.0)
set q_values=(0.5 0.4 0.3 0.2 0.1)

set nbody=25
set tcrit=10.0
set deltat=1.0

set par=nbody1.par
set log=nbody1.log


#                           my poor man's "keyword=value" command line parser 
#                           in principle all, but the array values, can be changed
foreach cmdarg ($*)
    set $cmdarg
end


#                           ! You should not have to change anything below this line !
#=====================================================================================

#                           make a master directory for all the experiments
mkdir -p $master

#                           save the script in the master directory
cp $0 $master

#                           loop over q's and a's
set count=0
foreach q ($q_values)
  foreach a ($a_values)
    # keep a counter (1,2,3,....)
    @ count++
    # set the working directory, and change into it
    set wdir=$master/run_q=${q}_a=${a}
    mkdir -p $wdir
    pushd $wdir

    # create the input file, for nbody1 in this case

#===========================================    nbodyX input file
cat > input << END_OF_INPUT
1 0.5
$nbody 1 200 1
0.01 $deltat $tcrit 2.0E-05 0.02
0 0 1 0 0 1 0 0 0 0 0 0 0 0 0
$a 5.0 1.0
$q 0.0 0.0 1.0 1.0
END_OF_INPUT
#===========================================

    # keep the user entertained about progress
    echo "Run# $count : q=$q a=$a "
    # make a parameter file for easy retrieval later
    echo a = $a  >  $par
    echo q = $q  >> $par
    # run program and record the unix shell status in the par file 
    time nbody1 < input >& $log
    set ok=$status
    echo ok = $ok >> $par
    # here you could analyse some files and record the results
    set nb=`cat $log | grep BINARY | wc -l`
    echo nb = $nb >> $par
    echo "   #binaries: $nb"
    # pop back for the next loop
    popd
  end
end


#  Example post-analysis:
#  Since this is a 2-dim parameter search, lets plot a simple matrix
#  of the number of binaries in the run


echo "==============================================================="
echo "BINARY activity:"
echo " a = IMF slope"
echo " q = Initial virial KE/PE"
echo ""
echo "   a   $a_values"
echo "q     ----------------------"
foreach q ($q_values)
  printf "%4s |" $q
  foreach a ($a_values)
    set wdir=$master/run_q=${q}_a=${a}
    set nb=`cat $wdir/$log | grep BINARY | wc -l`
    printf " %3d" $nb
  end
  echo ""
end
