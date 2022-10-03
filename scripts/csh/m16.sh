#! /bin/bash
#
#   compute m16 and add it to the "xv" table

# define the run
run=run0

# take all times that were in the run
times=$(tabcols $run/$run.xv.tab 1)

#             simple keyword=value command line parser for bash
for arg in $*; do
  export $arg
done


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


