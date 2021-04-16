#! /usr/bin/env bash
#! /usr/local/bin/bash
#! /bin/bash

#  System Integrity Protection:
#  bash from the system won't import them, /usr/local/bin is ok 

echo SHELL=$SHELL
echo NEMO=$NEMO
echo DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH

# on mac with a native bash 3.x on mac 10.15.7 rthis fails
# brew install bash  installed a mor modern 5.1 in /usr/local/bin
# this script works ok with invoked with /usr/local/bin/bash
#
# export DYLD_LIBRARY_PATH=/usr/lib

export NEMO=test1
export DYLD_LIBRARY_PATH=/usr/lib
