#
#! /usr/local/bin/bash
#! /usr/bin/env bash
#! /bin/bash

#  my frustating experience with SIP on a mac once my rogue laptop was updated to 10.15.7

echo SHELL=$SHELL
echo DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH

#  if you source this script, the DYLD variable should now be visible in your parent shell.
#  If you then run the script
#             ./sip.sh
#  the outcome depends on the first line in the script....
#
#      1.   #                                   you are fine (i don't yet understand this one)
#      2.   #! /usr/local/bin/bash              you are fine
#      3.   #! /usr/bin/env bash                you are SIP'd because a system tool filters is
#      4.   #! /bin/bash                        you are SIP'd because a system tool filters is
#
export DYLD_LIBRARY_PATH=/usr/lib


#  This is for bash, the results for a csh style script is the same.
