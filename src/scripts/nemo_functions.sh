#! /usr/bin/env bash
#
#   some functions to share for nemo bash scripts

_nemo_version="7-feb-2023"

echo "NEMO>> READING nemo_functions $_nemo_version via $0"

function nemo_version {
    v=$(cat $NEMO/VERSION)
    d=$(date -u +%Y-%m-%dT%H:%M:%S)
    g=$(cd $NEMO; git rev-list --count HEAD)
    h=$(uname -a)
    echo "$v  $g  $d  $h"
}

function printf_red {
    # could also use the tput command?
    RED='\033[0;31m'
    NC='\033[0m' # No Color
    echo -e "${RED}$*${NC}"
}

function printf_green {
    # could also use the tput command?
    RED='\033[0;32m'
    NC='\033[0m' # No Color
    echo -e "${RED}$*${NC}"
}

function show_vars {
    # helper function to show value of shell variables using bash dynamic variables
    for _arg in "$@"; do
	echo "${_arg}=${!_arg}"
    done
    
}

