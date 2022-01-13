#! /bin/csh -f 
#
#     can also try:    /bin/csh -f
#

set aaa = 1
set bbb = 2


foreach arg ($*)
    set $arg
end

echo "Computing with aaa=${aaa} and bbb=${bbb}"
