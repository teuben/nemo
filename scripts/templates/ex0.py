#! /usr/bin/env python
#
#

import sys


aaa = 1
bbb = 2

#   parse
for arg in sys.argv[1:]:
    exec(arg)

print(aaa,bbb)    
