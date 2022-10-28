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
if [ "$1" == "--help" ] || [ "$1" == "-h" ];then
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
echo "Computing with aaa=\"${aaa}\" and bbb=\"${bbb}\""
echo "And a possibly non-existing ccc=${ccc}"
echo "    - this should trigger an error if debug=1"
echo "All done."

#--HELP
#  this script will handle spaces in keyword values if properly quoted
#  but it is up to the script itself to properly handle it.
#  In generally it is not adviced to use spaces. 
#--HELP
