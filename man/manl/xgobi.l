.TH XGobi 1 "12 April 1994"
.SH NAME
xgobi \- interactive dynamic graphics program for data visualization
.SH SYNOPSIS
.B xgobi
[ X options ]
[
.B \-std
mmx|msd|mmd
] [
.B \-dev
std_deviation
] 
filename
.SH DESCRIPTION
.I xgobi
is an interactive dynamic graphics program for data
visualization in the X Window System.  It is especially
designed for the exploration of multivariate data.
Its basic plot is a scatterplot, and these are some
of the tools available for scatterplot display and
manipulation:
.RS .25i
.LP
Cycling rapidly through two-variable scatter plots.
.LP
Three-dimensional rotation, including trackball control.
.LP
Grand tour rotation: a smooth randomized sequence of
two-dimensional projections in order to explore a
higher-dimensional point cloud.
.LP
Brushing: the ability to change the color or plotting
character of a point or a group of points.
.LP
Scaling: moving and reshaping the plot interactively.
.LP
Identification: displaying an label next to a point
in the plot window.
.LP
Linking:  Brushing, identification and touring are linked;
that is, actions in the window of one XGobi process are
immediately reflected in another XGobi window displaying
the same data.
.LP
Line editing: Any pair of points can be connected
with a line segment.
.LP
Transformation:  A menu of data transformations is
available.
.LP
Postscript plotting:  The contents of the plot windows
can be written out to postscript files for high-quality
plotting.
.LP
Online help:  Most of the XGobi buttons and windows
have associated help files that can be displayed during
an XGobi session.  Click on the Help button toward the
top right of the window for instructions.
.RE
.LP
XGobi has a direct manipulation interface, and all the above
actions are performed using the mouse.
.LP
XGobi can be used in conjunction with S (1).  Execute
.IR help(xgobi) 
within S for information.
.SH OPTIONS
.sp
.B \-std
mmx|msd|mmd
.LP
.RS .25i
By default, the data are scaled into the plotting window
using the minimum and maximum values of each variable or
variable group, in such a way that the midpoint of the
variable is at the center of the plotting window and no
points fall outside the window.  Instead,
to scale using mean and standard deviation, specify 
.IR \-std
.IR msd;
to scale using the median and median absolute deviation,
specify
.IR \-std
.IR mmd.
.RE
.sp
.B \-dev
x
.RS .25i
.LP
If you have specified 
.IR -std
.IR msd
or
.IR -std
.IR mmd,
then you can also specify the number of standard deviations
(or median absolute deviations) from the mean (or median)
to be contained within the plotting window, using the argument
.IR -dev
.IR x,
where 
.IR x
is a real number between 0 and 100.  The default is 2.
.RE
.sp
X options
.sp
.RS .25i
The standard X command line options can be used with XGobi.
These include 
.IR "-display machinename:0",
used when running an X program on one machine and displaying
its output on another, and
.IR "-title Title",
where
.IR Title
is a string you want to appear in the window manager titlebar.
.RE
.SH FILES
XGobi accepts standard input, but is most often used with
files, partly because of the additional plot control that can
be achieved using a set of files.  The data input file
should be an ASCII file
with the data matrix arranged in rows and columns; rows must
be distinguished by carriage returns, and columns can be
separated by any amount of white space.  (The input file
can also be a binary file, which can be produced within
XGobi once the ASCII data has been read in.)
XGobi accepts
other input about the display of the data from files as well.
.P
If the data is in a file named
.B filename
or
.B filename.dat
(either of which must be an ASCII file), or
.B filename.bin
(the binary version of the data), then the other files are as follows:
'''
.sp
.nf
.B filename.row
.B filename.rowlab
.B filename.case
.fi
.RS .25i
.LP
Row or case labels: a label for each row of the data matrix,
which is displayed in the identification mode.  The file
should contain one label per line.
.LP
.RE
'''
.sp
.nf
.B filename.col
.B filename.collab
.B filename.column
.B filename.var
.fi
.RS .25i
.LP
Column or variable labels: a label for each column of the data
matrix, which becomes part of the XGobi variable selection panel.
The file should contain one label per line.
.RE
'''
.sp
.B filename.colors
.RS .25i
.LP
Brushing colors: a color for each point in the plot,
representing a row or case of the data.
The file should contain one color per line.  (It is probably
best if the colors correspond to the colors used in
brushing; see the later section on resources.)
.RE
'''
.sp
.B filename.glyphs
.RS .25i
.LP
Brushing glyphs: a glyph type for each point in the plot,
representing a row or case of the data.
The file should contain one glyph type per line.
The glyph types are as follows:
.RS .75i
.sp
.nf
1 through 5: Five sizes of '+'
6 through 10: Five sizes of 'X'
11 through 15: Five sizes of open rectangle
16 through 20: Five sizes of filled rectangle
21 through 25: Five sizes of open circle
26 through 30: Five sizes of filled circle
31: A single-pixel point
.fi
.RE
.RE
'''
.sp
.B filename.erase
.RS .25i
.LP
Erase: a column of 1s (to have a point erased on startup) and 0s
(to have the point plotted).  There should be one value per line
and as many lines as there are rows in the data.
.RE
'''
.sp
.B filename.lines
.RS .25i
.LP
Line segments: specifications for the pattern of
line segments which connect pairs of points.
The file should contain two numbers per line.  The pair
of numbers represents the row numbers of the two
points that should be connected.
.RE
'''
.sp
.B filename.linecolors
.RS .25i
.LP
Line colors: a color for each line in the .lines file.
The file should contain one color per line.  (It is
best if the colors correspond to the colors used in
brushing; see the later section on resources.)
.RE
'''
.sp
.B filename.nlinkable
.RS .25i
.LP
The number of rows to be linked for brushing and
identification.  By default, nlinkable is equal to the
number of rows in the data.  This feature can be used to link
ordinary scatterplots with plots that have some decorations
requiring additional points, such as clustering trees.
.RE
'''
.sp
.B filename.vgroups
.RS .25i
.LP
Variable groups: an integer for each column in the data.
each set of columns that is represented by the same integer 
will grouped together for scaling and transformation.
The file is just one long line of integers.
For example, an input file with four columns could have a
.vgroups file containing the line 1 2 2 3.  The second
and third columns are then grouped.  The range
of their plotting axes is be the same, and if column
2 is transformed, column 3 is transformed at the
same time.
.RE
'''
.sp
.B filename.resources
.RS .25i
.LP
Resources: a set of datafile-specific XGobi resources,
which specify the size of the plotting window and some
user-selection option settings.
The file is in the format of a standard X resource file.
It can be directly edited so that other resources
can be specified.
See the later section on resource files for more
information.
.RE
'''
.LP
All of the above files can be created outside of XGobi,
using an editor or other UNIX utilities, and several
of them (glyphs and colors, line segments and line color,
resources)
can be stored during an XGobi session, in which case they
represent the results of interactions performed during
that session.
.SH "X11 RESOURCE FILES"
X11 resource files are used to control size, color, layout,
labels, and various other windowing features.  The XGobi
resources you are most likely to find useful set the sizes
of the plotting window and the variable selection panel,
the colors available for brushing, the font types and
sizes, and the title appearing in the title bar of the X
window.  Some examples follow:
.LP
In this example, the first 5 of the 10 brushing colors are
specified.  These colors should be drawn from the set of
colors available on your workstation or X terminal, which
you can probably determine using the X command 'showrgb.'
.RS .5i
.nf
*brushColor0: white
*brushColor1: red
*brushColor2: green
*brushColor3: skyblue
*brushColor4: yellow
.fi
.RE
.LP
Here, the sizes of the plotting window and variable selection
panel are set.  All units are in pixels.
This file was created interactively during
an XGobi session.
.RS .5i
.nf
*XGobi*PlotWindow.height: 120
*XGobi*PlotWindow.width: 490
*XGobi*VarPanel.width: 390
.fi
.RE
.LP
Here, the title and font have been adjusted for making a
videotape.  An explanatory title and a large font are used.
.RS .5i
.nf
*XGobi.title:    XGobi: Experimental Laser Data
*XGobi.iconName: XGobi: Laser
*plotFont: -*-lucida-medium-r-normal-sans-20-*-*-*-*-*-*-*
.fi
.RE
.LP
Several logical variables controlling XGobi's behavior
can be defined in a resource file.  These are also
represented by buttons, most of which are found on the
Option menu accessed from the main panel.  Descriptions
can be found in the help window for each button and in
the full-length user's manual.
.RS .5i
.nf
*showAxes: True
*showPoints: True
*showLines: True
*linkBrush: True
*linkLineBrush: True
*linkIdentify: True
*linkTour: True
*plotSquare: True
*jumpBrush: True
*reshapeBrush: True
*carryVars: True
*glyphType: 6
*glyphSize: 1
*defaultPrintCmd: lpr -Pps1
.fi
.RE
.SH AUTHORS
Deborah Swayne (dfs@bellcore.com), Dianne Cook
(dcook@iastate.edu), and Andreas Buja
(andreas@bellcore.com)
.SH CONTACT
Deborah Swayne (dfs@bellcore.com)
.SH SEE ALSO
.LP
Sample data files are included in the XGobi distribution.
.LP
Deborah Swayne, Dianne Cook and Andreas Buja (1991),
\f1User's Manual for XGobi, a Dynamic Graphics Program for
Data Analysis Implemented in the X Window System (Version 2)\fP,
Bellcore Technical Memorandum.
.SH BUGS
Probably too numerous to mention, but you would do the
authors a favor if you report the ones you notice.

