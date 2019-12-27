#! /bin/csh -f
#
#   our poor man's command line parser in csh
#

set a=1
set b=2
set c=3

foreach arg ($*)
   set $arg
end

echo a=$a b=$b c=$c
