#! /usr/bin/env bash
#
#     can also try:    /bin/bash
#

#--HELP
#       set defaults
aaa=1
bbb=2
debug=0
#--HELP

#       give help?
if [ "$1" == "--help" ];then
    set +x
    awk 'BEGIN{s=0} {if ($1=="#--HELP") s=1-s;  else if(s) print $0; }' $0
    exit 0
fi

#       parse the commandline
for arg in "$@"; do
    export "$arg"
done

#       some debugging?
if [ $debug != 0 ]; then
    set -x
    set -e
    set -u
fi

#       do the work
echo "Computing with aaa=${aaa} and bbb=${bbb}"
echo "And a possibly non-existing ccc=${ccc}"
