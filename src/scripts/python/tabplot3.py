#! /usr/bin/env python
#
#  tabplot in python

import os
import sys
import argparse
import numpy as np
import numpy.ma as ma
import pandas as pd
import json
import matplotlib.pyplot as plt
from io import StringIO

#import pdb
#pdb.set_trace()

isScatter = False
numlines = 0

argv = sys.argv[1:]


# Defining parser
# add_argument attaches individual assignment specifications to the parser



# Add help
parser = argparse.ArgumentParser(
    description = "plots a table like tabplot", 
     formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=40, width=100) #Widens the available width so it all prints on the same line
)

parser.add_argument('-i', '--input', default = "", help = "inupt fits cube") # -i/--input argument

parser.add_argument('-x', '--xcol', default='1', type=str, help = "X column(s)")
parser.add_argument('-y', '--ycol', default='2', type=str, help = "Y column(s)")
parser.add_argument('-c', '--color', default="black,", type=str, help = "Line Width(s)")
parser.add_argument('-l', "--line", default="0,", type=str, help = "Color(s)")   
parser.add_argument('-p', "--point", default="1,", type=str, help = "Point size(s)")  
parser.add_argument('-lb', "--labels", default=",", type = str, help = "Labels for each plot (required for legend)")
parser.add_argument('-lg', "--legend", action = "store_true", help = "Add a legend")
parser.add_argument('-t', '--title', default = "", type = str, help = "Add a title")
parser.add_argument('-o', '--out', default='', type=str, help = "Save to file")



#Gets defaults - needed for when there's errors with extracting arguments
def get_defaults():
    defaults = {}
    for action in parser._actions:
        if action.dest != 'help':  # skip help
            defaults[action.dest] = action.default
    return defaults



# Store parsed results + handle pipes for missing results
# Manually splitting here so -i or --input can be set to "" when needed for piping
avs = []
sublist = [argv[0]]
for el in argv[1:]:
    if el == '-i' or el == '--input':
        av = parser.parse_args(sublist)
        avs.append(av)
        sublist = []
    
    sublist.append(el)

avs.append(parser.parse_args(sublist))



plt.figure()
out = False; outpath = None; title = None


#Creating a separate plot for each varset
for varset in avs:
    p = vars(varset) 

    #Declare cli params
    infile = p['input']
    xinput = p['xcol']
    yinput = p['ycol']
    colorinput = p['color']
    lineinput = p['line']
    pointinput = p['point']
    labels = p['labels']
    legend = p['legend']
    out_arg = p['out']
    ttl = p['title']


    #Handling pipe by making infile a readble stream if not detected
    if not infile or infile == "":
        if not sys.stdin.isatty():
            data_str = sys.stdin.read()
            infile = StringIO(data_str)
        else:
            parser.print_help()
            sys.exit()


    #Saving output path if we want it (title too)
    if out_arg != "" and not out:
        out = True; outpath = out_arg
    
    if title == None:
        title = ttl



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
    colors = [i for i in colorinput.split(',') if i != ""]
    linewidths = [i for i in lineinput.split(',') if i != ""]
    labels = [i for i in labels.split(",") if i != ""]

    pointsizes = [1] * numlines
    try:
        pointsizes = pointinput.split(',')
        pointsizes = pointsizes[:numlines] if len(pointsizes) > numlines else pointsizes + [1 for _ in range(numlines - len(pointsizes))]
        pointsizes = [int(i) for i in pointsizes]
    except ValueError:
        pointsizes = [1] * numlines

    # if no xcol or ycol specified, assume col 1 is x and col 2 is y
    if len(xcols) == 0:
        xcols = [1]
    if len(ycols) == 0:
        ycols = [2]

    # set non-specifics and trim the length if needed
    colors = colors[:numlines] if len(colors) > numlines else colors + ["black" for _ in range(numlines - len(colors))]
    linewidths = linewidths[:numlines] if len(linewidths) > numlines else linewidths + [0 for _ in range(numlines - len(linewidths))]
    labels = labels[:numlines] if len(labels) > numlines else labels + ["" for _ in range(numlines - len(labels))]


    # gather data 

    df = pd.read_csv(infile, sep=None, engine='python', header=None)
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.select_dtypes(include="number").dropna(axis=1, how='all')
    data = df.values.T


    #data = np.loadtxt(infile).T  #w np

    #print(f"\nRead the following:\n{data}\n")

    xdata = [0] * len(xcols)
    ydata = [0] * len(ycols)

    for i in range(len(xcols)): # Read x data
        xdata[i] = data[xcols[i]-1]

    for i in range(len(ycols)): # Read y data
        ydata[i] = data[ycols[i]-1]

    # Plotting

    if len(xdata) == 1: # Case: only 1 xcol, plot each y against the only x
        for i in range (len(ydata)):
            plt.plot(xdata[0], ydata[i], color=colors[i], linewidth=linewidths[i], marker = '.', markersize=pointsizes[i], label = labels[i])

    elif len(ydata) == 1: # Case: only 1 ycol, plot each x against the only y
        for i in range (len(xdata)):
            plt.plot(xdata[i], ydata[0], color=colors[i], linewidth=linewidths[i], marker = '.', markersize=pointsizes[i], label = labels[i])

    else:               # Case: more than 1 xcol and ycol, plot each ycol with its xcol at the same index
        for i in range (len(xdata)):
            try:
                plt.plot(xdata[i], ydata[i], color=colors[i], linewidth=linewidths[i], marker = '.', markersize=pointsizes[i], label = labels[i])
            except IndexError:
                print("ERROR: xcols and ycols mismatch")
                sys.exit()

    if legend:
        plt.legend()

#Showing all figures at once (or saving)
if title is not None:
    plt.title(title)

if out:
    plt.savefig(outpath)

else:
    plt.show()
