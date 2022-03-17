#! /usr/bin/env bash
#
#     can also try:    /bin/bash
#

#       set defaults
aaa=1
bbb=2
debug=0

#       parse the commandline
for arg in $*; do
    export $arg
done

#       some debugging
if [ $debug != 0 ]; then
    set -x
    set -e
fi

#       do the work
echo "Computing with aaa=${aaa} and bbb=${bbb}"
echo "And ccc=${ccc}"
