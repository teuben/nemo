#! /bin/bash
#
#  NEMO's simple frontend for sbatch 
#
#  SLURM cheat list 
#     sinfo
#     sbatch run_12345.sh               (this example)
#     squeue -u $USER                  (also shows your JOBID's)
#     scancel JOBID
#     srun -n 1 -c 4 --mem=16G -p toltec-cpu --x11 --pty bash
#

# Usage
#     sbatch_nemo.sh -x scriptfile [args]
#     sbatch_nemo.sh cmdline.txt

#                                --help
if [ "$1" == "--help" ];then
    echo Usage:...
    exit 0
fi

#                                -l cmds.txt
if [ "$1" == "-l" ];then
    shift
    # catch the single argument batch call first
    if [ -e "$1" ]; then
	echo Processing lines from $1 line by line
	while IFS= read -r line; do
            echo "LINE: $line"
            sbatch_nemo.sh $line
	done < $1
	exit 1
    fi

    exit 0
fi




# processing CLI when key=var
for arg in $*; do
    if [[ "$arg" =~ "=" ]]; then
	export $arg
    fi
done

#
runid=$$

#                                        version
version="14-may-2022"

#                                        prefix to run
prefix="/usr/bin/time xvfb-run -a"

#                                        sbatch run file
run=run_$runid.sh

#                                        max sbatch time 
tmax=04:00:00



if [ "$(which sbatch)" != "/usr/bin/sbatch" ]; then
    echo "$0 version=$version"    
    echo "run=$run"
    echo "ERROR:  No sbatch system here on $(hostname)"
    #exit 0
fi


cat <<EOF > $run
#! /bin/bash
#
#   $0 version=$version
#
#SBATCH -o slurm-%j-%x.out
#SBATCH -t $tmax
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

$prefix $*

EOF


chmod +x $run
echo "$run      - use scancel JOBID to kill this one, JOBID is:"
sbatch $run
