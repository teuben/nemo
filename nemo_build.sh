#! /bin/bash
#
#   runs configure and all three build phases, outputting results into one log file.

# ============================================================================================

#    Poor man's command line parser
for arg in $*; do\
  export $arg
done

#   Run configure step and build1+2
./configure
make build1
make build2

#   Move config.log contents into install.log.
cat config.log install.log > install.tmp && mv install.tmp install.log

#   Finally, run build3 and check.
make build3
make check | tee -a install.log
