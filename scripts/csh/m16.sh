#! /bin/bash
#
#--HELP
#   compute m16 and add it to the "xv" table -  helper script after mkmh97.sh
#
#   Timing:   Handling 200 time frames of nbody=10000  took about 9 mins.
#                      100 time frames of nbody=100000 took about 7 mins.

# define the run (required)
run=run0

# times?   by default all times from the simulation are used, via run.xv.tab
#          otherwise manually use a comma separated list here
unset times

#--HELP

if [ -z $1 ] || [ "$1" == "--help" ] || [ "$1" == "-h" ];then
    set +x
    awk 'BEGIN{s=0} {if ($1=="#--HELP") s=1-s;  else if(s) print $0; }' $0
    exit 0
fi

# simple keyword=value command line parser for bash
for arg in $*; do
  export $arg
done

# make sure directory exist
if [ ! -d $run ]; then
    echo run directory run=$run does not exist
    exit 0
fi

# take all times that were in the run
if [ -z "$times" ]; then
    times=$(tabcols $run/$run.xv.tab 1)
fi

new=$run/$run.xvm.tab
rm -f $new

for t in $(echo $times | sed 's/,/ /g'); do
    echo $t
    ./mkmh97.sh run=$run tstop=$t  > /dev/null 2>&1
    m16=$(nemopars m16 $run/nemopars.rc | grep -v ^#)
    if [ -z "$m16" ]; then
	m16=0
    fi
    xv=$(grep  ^"$t " $run/$run.xv.tab)
    echo $xv $m16 >> $new
done


