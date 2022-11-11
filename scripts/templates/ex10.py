#! /usr/bin/env python
#
# NEMO's getparam

import math
from nemopy import getparam

keyval = [
    "aaa=1\n       The aaa (int) variable",
    "bbb=10.0\n    The bbb (float) variable",
    "VERSION=1\n   11-nov-2022 PJT",
    ]

usage = "showing off getparam in python"

if __name__ == '__main__':

    p = getparam.Param(keyval, usage)
    print('aaa:',p.get("aaa"))
    print('bbb:',p.get("bbb"))


#  alternatively we can hide the objects even more?
#        nemo.nemo_main(keyval,usage)
#        nemo.getparam("aaa")
#        nemo.getrparam("bbb")

