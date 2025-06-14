#! /usr/bin/env python
#
#  tabplot in python

import os
import sys
import argparse
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

#import pdb
#pdb.set_trace()

isScatter = False
numlines = 0

argv = sys.argv[1:]
inputs = [i for i, x in enumerate(argv) if x == '-i' or x == '--input']

# loop over the input files
nf = len(inputs)
if nf == 1:
    inputs2 = [len(argv)]
else:
    inputs2 = inputs[1:] + [len(argv)]

# Defining parser
# add_argument attaches individual assignment specifications to the parser

# Add help
parser = argparse.ArgumentParser(
    usage = """ 
    -i --input ???       input fits cube\n
    -x --xcol 1       X column(s) (1=first column)\n
    -y --ycol 2       Y column(s)\n
    -l --line type,width        TBD\n
    -c --color                  TBD\n
    -p --point type,size       TBD\n
    -VERSION 0.1     3-apr-2023 PJT
    """,

    description = "plots a table like tabplot"
)

parser.add_argument('-i', '--input', required=True) # -i/--input argument
parser.add_argument('-x', '--xcol', default='1', type=str)
parser.add_argument('-y', '--ycol', default='2', type=str)
parser.add_argument('-c', '--color', default="", type=str)
parser.add_argument('-l', "--line", default="", type=str)     #Type, width
parser.add_argument('-p', "--point", default="", type=str)  #Type, size

# Store parsed results
avs = []
for i in range(len(inputs)):
    section_args = argv[inputs[i]:inputs2[i]]
    av = parser.parse_args(section_args)
    avs.append(av)

#Show help if input not provided
if(len(avs) == 0):
    parser.print_help()
    sys.exit()

p = vars(avs[0])

#Declare cli params
infile = p['input']
xinput = p['xcol']
yinput = p['ycol']
colorinput = p['color']
lineinput = p['line']
pointinput = p['point']




# get indexes for xcol and ycol
try:
    xcols = [int(x) for x in xinput.split(',')]
    ycols = [int(y) for y in yinput.split(',')]
except ValueError:
    print("ERROR: excpected 'ints' or 'ints,ints")
    exit(1)
if(len(xcols) > len(ycols)):
    numlines = len(xcols)
else:
    numlines = len(ycols)

# get CLI parameters
colors = colorinput.split(',')
linewidths = lineinput.split(',')
pointsizes = [1] * numlines
try:
    pointsizes = [int(p) for p in pointinput.split(',')]
except ValueError:
    pointsizes = [1] * numlines

# if no xcol or ycol specified, assume col 1 is x and col 2 is y
if len(xcols) == 0:
    xcols = [1]
if len(ycols) == 0:
    ycols = [2]

# if no color or width specified, set color to black and width 1
if colors[0] == "":
    colors = ["black"] * numlines
if linewidths[0] == "":
    linewidths = [1] * numlines

# gather data
data = np.loadtxt(infile).T

xdata = [0] * len(xcols)
ydata = [0] * len(ycols)

for i in range(len(xcols)): # Read x data
    xdata[i] = data[xcols[i]-1]

for i in range(len(ycols)): # Read y data
    ydata[i] = data[ycols[i]-1]

# Plotting
plt.figure()
if len(xdata) == 1: # Case: only 1 xcol, plot each y against the only x
    for i in range (len(ydata)):
        plt.plot(xdata[0], ydata[i], color=colors[i], linewidth=linewidths[i], marker = '.', markersize=pointsizes[i])
elif len(ydata) == 1: # Case: only 1 ycol, plot each x against the only y
    for i in range (len(xdata)):
        plt.plot(xdata[i], ydata[0], color=colors[i], linewidth=linewidths[i], marker = '.', markersize=pointsizes[i])
else:               # Case: more than 1 xcol and ycol, plot each ycol with its xcol at the same index
    for i in range (len(xdata)):
        try:
            plt.plot(xdata[i], ydata[i], color=colors[i], linewidth=linewidths[i], marker = '.', markersize=pointsizes[i])
        except IndexError:
            print("ERROR: xcols and ycols mismatch")
            sys.exit()

plt.show()
