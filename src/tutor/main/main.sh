#! /bin/bash
#
#   our poor man's command line parser in bash

a=1
b=2
c=3
yapp=xs

for arg in $*; do\
   export $arg
done

# a useful YAPP_PGPLOT and YAPP_PS device function
yapp() {
  if test $yapp = "xs"; then
     echo $1/$yapp
  elif test $yapp = "_ps"; then
     echo fig$1.ps
  else
     echo fig$1.$yapp/$yapp
  fi
}


echo a=$a b=$b c=$c yapp=$(yapp 10)
