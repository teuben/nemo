.TH X3D 1 "15 August 1995"
.SH NAME
x3d \- interactive (X-windows based) 3d viewer
.SH OPTIONS
.nf
Options:
-display  <display>        X display name
-flicker                   Make window flicker
-geometry <geometry>       Window geometry
-help                      Print this message
-menu                      Turn off menu
-stats                     Turn off stats
-stereo                    Turn off stereo
-controls                  Switch controls style
-scale <amount>            Scale object by amount
-print <filename>          Postscript File for
                           writing lines into
-hpgl                      Write an HPGL plotter
                           file instead
.fi
.SH X3D DATA FORMAT
The conversion utility called con that takes a list of segments
.nf
x1 y1 z1 x2 y2 z2
. . .
x1 y1 z1 x2 y2 z2
.fi
and converts them into the .obj format (see Q_and_A for more information)
.SH CONTROL
.nf

  QUIT                    UPPERCASE AND LOWERCASE Q
  TOGGLE CONTROLS MENU    UPPERCASE AND LOWERCASE M
  TOGGLE POSITION DISPLAY UPPERCASE AND LOWERCASE P
  TOGGLE STEREO DISPLAY   UPPERCASE AND LOWERCASE S
  TOGGLE FLICKER          UPPERCASE AND LOWERCASE F
  ROTATE OBJECT ABOUT X   VERTICAL   POINTER MOVEMENT WHILE HOLDING LEFT BUTTON
  ROTATE OBJECT ABOUT Z   HORIZONTAL POINTER MOVEMENT WHILE HOLDING LEFT BUTTON
  MOVE OBJECT BACKWARD    UPPERCASE AND LOWERCASE J
  MOVE OBJECT FOREWARD    UPPERCASE AND LOWERCASE K
  MOVE OBJECT RIGHT       UPPERCASE AND LOWERCASE H
  MOVE OBJECT LEFT        UPPERCASE AND LOWERCASE L
  MOVE OBJECT UP          UPPERCASE AND LOWERCASE I
  MOVE OBJECT DOWN        UPPERCASE AND LOWERCASE U
  ROTATE ABOUT X          UPPERCASE AND LOWERCASE X A
  ROTATE ABOUT Y          UPPERCASE AND LOWERCASE Y B
  ROTATE ABOUT Z          UPPERCASE AND LOWERCASE Z C
  AUTOROTATE ABOUT X      1 2 3
  AUTOROTATE ABOUT Y      4 5 6
  AUTOROTATE ABOUT Z      7 8 9
  ADJUST FOCUS            SQUARE AND CURLY BRACKETS
  ADJUST FOCUS            SQUARE AND CURLY BRACKETS
  WRITE TO FILE           UPPERCASE AND LOWERCASE D
.fi
.SH BUGS
If the program dies,  try
.nf
	x3d -menu -stats -geometry 100x100+0+0 saturn.obj
.fi
and see if it might work now.
.SH SEE ALSO
snap3dv(1NEMO)
.SH AUTHOR
Mark Spychalla (spy@castlab.engr.wisc.edu)
University of Wisconsin - Madison
.SH HISTORY
.nf
15-aug-95	Made up a (NEMO) manual page
.fi
