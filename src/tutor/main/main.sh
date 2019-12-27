#! /bin/bash
#
#   our poor man's command line parser in bash

a=1
b=2
c=3

for arg in $*; do\
   export $arg
done

echo a=$a b=$b c=$c
