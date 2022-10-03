#! /bin/bash
#
#   compute m16 and add it to the "xv" table

# define the run
run=run701

# take all times that were in the run
times=$(tabcols $run/$run.xv.tab 1)
# or cheat for debugging
#times="20 60 80 100"
#times="20 60"

new=$run/m16.tab
rm -f $new

for t in $times ; do
    echo $t
    ./mkmh97.sh run=$run tstop=$t 
    m16=$(nemopars m16 $run/nemopars.rc | grep -v ^#)
    if [ -z "$m16" ]; then
	m16=0
    fi
    xv=$(grep -w ^$t $run/$run.xv.tab)
    echo $xv $m16 >> $new
done


