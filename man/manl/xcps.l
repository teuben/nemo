.TH XCPS 1 "Release 4" "X Version 11"
.SH NAME
xcps \- Colour PostScript dump facility for X11.4
.SH SYNOPSIS
.B xcps
[-option ...]
.SH DESCRIPTION
The \fIxcps\fP program is a window dumping facility for the X Window System.
It dumps windows directly in Colour or B/W PostScript format to standard
output. By default, \fIxcps\fP creates Colour-PS. The window to be dumped
is defined by clicking the mouse in that window, which can also be the root.
A beep indicates that dumping has started. Until this is fininshed, indicated
by two beeps, the window should not be changed. The resulting PS-image
will contain the window title, if available, the time of creation and the
name of the user. For more info check the options.
.SH OPTIONS
.TP 8
.B \-help
This causes \fIxcps\fP to print out a verbose message describing its options.
.TP 8
.B \-bw
This option indicates that the output should be grayscale-PS for 
B/W-PostScript printers instead of the default Colour-PS.
.TP 8
.B \-land
This option indicates that the output should be in LandScape mode instead
of the default Portrait mode.
.TP 8
.B \-curs
This option indicates that the user wants to specify a rectangular area
anywhere on the screen to be dumped instead of a complete window.
The first corner is defined by pressing a mouse button. While keeping the
button pressed, move to the opposing corner. A dragging rectangle, defining
the area to be dumped will appear. By releasing the mouse button the second
corner is defined.
.SH OPTIONS
.TP 8
.B \-display \fIdisplay\fP
This argument allows you to specify the server to connect to;
see \fIX(1)\fP.
.TP 8
.B \-x \fIx\fP
Specify x-dimension of dump. 
.TP 8
.B \-y \fIy\fP
Specify y-dimension of dump. 
.TP 8
.B \-bdrs
If a single window is dumped, include the X-window border.
.TP 8
.B \-frame
If a single window is dumped, include the window manager title bar.
.TP 8
.B \-noax
Do not draw a frame around the dumped image.  
.TP 8
.B \-text \fIcomment string\fP
This option specifies a comment string, to be printed above the dumped
image. If the comment contains spaces, enclose it within double quotes.
.TP 8
.B \-notext 
This option specifies that al text (user, date, window) should be omitted.
.TP 8
.B \-name \fIwindow name\fP
This option indicates the name of the window to be dumped. This avoids the
use of the pointer.
.TP 8
.B \-debug
Print some debugging information.

.SH "EXAMPLE"
\fIxcps -bw -curs > file.ps\fP creates a greyscale PS file. The area to be 
dumped will be defined by the user.

.SH Error Messages
Possible errors messages deal with the simultaneous specification of eg. 
-curs and -frame or -brds, or dimensions exceeding the physical paper size.
.SH "BUGS"
None noted sofar, but according to Murphy's law there definitely will be some.
.SH "NOTES"
The core-code of \fIxcps\fP was created by merging parts of \fIxwd\fP and
the PostScript library from \fIMIDAS\fP ((c) European Southern Observatory).
.SH "AUTHORS"
Edwin Huizinga, Kapteyn Lab Groningen (The Netherlands)

xwd: Tony Della Fera, Digital Equipment Corp., MIT Project Athena
William F. Wyatt, Smithsonian Astrophysical Observatory

MIDAS library: Preben Grosbol, ESO/IPG
     

