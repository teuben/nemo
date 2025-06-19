#! /bin/bash
#
#   runs configure and all three build phases, outputting results into one log file.
#
#   although this works, this is not normally how you want to install NEMO, instead
#   look at the docs/nemo_install.sh script for more flexibilty. 

# ============================================================================================

#   Run configure step and build1+2
./configure
make build1
make build2

#   Move config.log contents into install.log.
cat config.log install.log > install.tmp && mv install.tmp install.log

#   Finally, run build3 and check.
make build3
make check | tee -a install.log
