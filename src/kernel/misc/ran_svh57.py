#! /usr/bin/env python

#   random number generator according to SvH 1957,
#   as explained in vH1960. In NEMO also as ran_svh57() in xrand.c
#   and use xrandom to test the C version
#
#  Comparing the python and C version:
#      ./ran_svh57.py   0.580128 1001 > tab1
#      xrandom n=1000 seed57=0.580128 > tab2



import math
import sys

# pick number between 0.57 and 0.91

xlo = 0.57
xhi = 0.91

seed = 0.580128
n    = 10
if True:
    seed = float(sys.argv[1])
    n    = int(sys.argv[2])


x = seed
for n in range(n):
    x2 = x*x
    # eq.(44) in 1960ZA.....50..184V
    if x2 < xlo:
        xnew = x2 + 0.32
    else:
        xnew = x2
    # Von jeder Zahl x_n dieser Folge ist 
    # das erste viertel ihrer ziffern zu streichen,
    # der Rest ist als Zufallzahl zu benutzen.
    #
    # My interpretation:   the Siemens 2002 https://dl.acm.org/doi/10.1145/1458043.1458076
    # has a word length of 12 decimals, so lets cut out the first 3 and shift the decimal
    # point by that amount. It means the resulting random number has 9 digits
    # e.g. 0.123456789012 becomes 0.456789012
    if True:
        ran = int(('%.11f' % xnew)[5:])/1e8        # true chopping of 11 digits
        ran = int(('%.12f' % xnew)[5:])/1e9        # true chopping of 12 digits
    else:
        ran = math.modf(xnew*1000)[0]              # full math
    print(ran)
    x = xnew
