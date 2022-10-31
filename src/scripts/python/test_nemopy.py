#! /usr/bin/env python
#
import os, sys

from nemopy import getparam

keydef = [
    "a=1\n         the a parameter",
    "b=2\n         the b parameter",
    "VERSION=1\n   30-oct-2022 PJT",
    ]

usage = """
  this is a test program for nemopy.

  Is has just two parameters, a= and b=.
  """

p = getparam.Param(keydef,usage)

print("a as string:",p.get("a"))
print("b as list:  ",p.listf("b"))


