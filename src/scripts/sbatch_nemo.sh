#! /bin/bash
#
#--HELP
#  NEMO's simple frontend for sbatch. Only tested on zaratan.umd.edu
#
# Usage
#     sbatch_nemo.sh -x scriptfile [args]
#     sbatch_nemo.sh -r runfile.txt [args]
#
# If a file sbatch_nemo.rc is present in the current directory,
# it will be sourced to override any parameters used here.
# A file $SBATCH_TEMPLATE is also allowed after this.
#
# SLURM cheat list 
#     sinfo
#     sbatch run_12345.sh              (repeat a run)
#     squeue -u $USER                  (also shows your JOBID's)
#     scancel JOBID
#     sbalance
#
version=21-sep-2022
#--HELP


# give help
if [ $# -eq 0 ] || [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    set +x
    awk 'BEGIN{s=0} {if ($1=="#--HELP") s=1-s;  else if(s) print $0; }' $0
    exit 0
fi

#                                -l runfile.txt
if [ "$1" == "-r" ]; then
    shift
    # catch the single argument batch call first
    if [ -e "$1" ]; then
	echo Processing lines from $1 line by line
	while IFS= read -r line; do
            echo "LINE: $line"
            sbatch_nemo.sh -x $line
	done < $1
	exit 1
    fi
    exit 0
fi

if [ "$1" == "-x" ]; then
    shift
fi



# processing CLI when key=var
for arg in $*; do
    if [[ "$arg" =~ "=" ]]; then
	export $arg
    fi
done

#                                        probably don't change this
runid=$$
#                                        prefix to run
prefix="/usr/bin/time xvfb-run -a"
#                                        sbatch run file
run=run_$runid.sh
#                                        max sbatch time 
tmax=04:00:00
#                                        max memory used
mem=16G


if [ "$(which sbatch)" != "/usr/bin/sbatch" ]; then
    echo "$0 version=$version"    
    echo "run=$run"
    echo "ERROR:  No sbatch system here on $(hostname)"
    exit 0
fi

if [ -e sbatch_nemo.rc ]; then
    source sbatch_nemo.rc
fi

if [ "$SBATCH_TEMPLATE" != "" ]; then
    source $SBATCH_TEMPLATE
fi

cat <<EOF > $run
#! /bin/bash
#   Probably should not edit this file, it has been created by $0
#   $0 version=$version
#
#SBATCH -o slurm-%j-%x.out
#SBATCH -t $tmax
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=$mem

$prefix $*

EOF


chmod +x $run
echo "$run      - use scancel JOBID to kill this one, JOBID is:"
sbatch $run
