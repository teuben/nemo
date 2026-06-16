#! /usr/bin/env python
#
#    plot one or more RSR spectra, optionally a band or piece of the spectrum
#

import sys
import numpy as np
from docopt import docopt
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons

_version = "14-mar-2024"

_help = """Usage: tab_plot.py [options] TABLE1 [TABLE2...]

Options:
  -y PLOTFILE             Save plotfile instead of default interactive. Optional.
  -t --title TITLE        Title of plot. Optional.
  --xscale XSCALE         Scale factor to apply to X axis [Default: 1.0]
  --yscale YSCALE         Scale factor to apply to Y axis [Default: 1.0]
  --xlab XLAB             X-axis label [Default: X]
  --ylab YLAB             Y-axis label [Default: Y]
  --xrange XMIN,XMAX      Min and Max of X-coordinate. Optional.
  --yrange YMIN,YMAX      Min and Max of Y-coordinate. Optional.
  --irange IMIN,IMAX      Min and Max of X-coordinates on top axis (twiny). Optional.
  --size SIZE             Plotsize in inches. [Default: 8.0]
  --dpi DPI               DPI of plot. [Default: 100]
  --ycoord YCOORD         Dashed line at YCOORD. Optional.
  --boxes LL,UR           4-Tuples of lower-left upper-right boxes. Optional
  -h --help               This help
  -d --debug              More debugging output
  -v --version            The script version

One or more ASCII tables are needed, with columns 1 and 2 designating the
X and Y coordinates.
"""

args = docopt(_help, options_first=True, version=f'tab_plot.py {_version}')
print("Arguments:", args)

# Extract command-line arguments
xscale = float(args.get('--xscale', 1.0))
yscale = float(args.get('--yscale', 1.0))
size = float(args.get('--size', 8.0))
dpi = int(args.get('--dpi', 100))
title = args.get('--title', 'Spectral Plot')
xlab = args.get('--xlab', 'X')
ylab = args.get('--ylab', 'Y')
xmin = xmax = ymin = ymax = None
ycoord = None
plotfile = args.get('-y', None)
imin = imax = None

if args['--xrange']:
    xmin, xmax = map(float, args['--xrange'].split(','))
if args['--yrange']:
    ymin, ymax = map(float, args['--yrange'].split(','))
if args['--ycoord']:
    ycoord = float(args['--ycoord'])
if args['--irange']:
    imin, imax = map(float, args['--irange'].split(','))

# Load data
spectra = [args['TABLE1']] + args.get('TABLE2', [])
data = []
labels = []
lines = []

for filename in spectra:
    d = np.loadtxt(filename).T
    data.append(d)
    labels.append(filename)

# Create the plot
fig, ax = plt.subplots(figsize=(size, size), dpi=dpi)
plt.subplots_adjust(left=0.3)

for d, label in zip(data, labels):
    line, = ax.step(d[0] * xscale, d[1] * yscale, label=label, where='mid', linestyle=':')
    lines.append(line)

if xmin is not None and xmax is not None:
    ax.set_xlim(xmin, xmax)
if ymin is not None and ymax is not None:
    ax.set_ylim(ymin, ymax)

ax.set_xlabel(xlab)
ax.set_ylabel(ylab)
ax.set_title(title)

# Interactive CheckButtons
# Adjust the position to avoid overlapping with the main plot area
rax = plt.axes([0.01, 0.5, 0.2, 0.3], facecolor='aliceblue')  

# Styling the CheckButtons
check = CheckButtons(rax, labels, [True] * len(labels))

# Change colors and style of checkboxes
for rec, text in zip(check.rectangles, check.labels):
    rec.set_facecolor('lightblue')       # Color of the checkbox
    rec.set_edgecolor('darkblue')        # Edge color of the checkbox
    text.set_color('darkblue')           # Text color
    text.set_fontsize(10)                # Font size
    text.set_fontweight('bold')          # Font weight

# Function to toggle visibility of the data lines
def func(label):
    index = labels.index(label)
    lines[index].set_visible(not lines[index].get_visible())
    plt.draw()

check.on_clicked(func)


# Optional dashed line at y-coordinate
if ycoord:
    ax.axhline(y=ycoord, color='gray', linestyle='--')

# Save or show the plot
if plotfile:
    plt.savefig(plotfile)
else:
    plt.show()
