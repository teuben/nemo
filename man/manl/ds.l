.\" @(#)ds.1 1.0 87/11/13 WLS;
.\" $Header$
.TH DS 1 "13 Nov 1987"
.SH NAME
ds \- a window and mouse-based utility for displaying a picture
.SH SYNOPSIS
.B ds 
[
.B \-debug
]
[
.B \-d
]
[
.B \-mono
]
[
.B \-m
]
[
.B \-xflip
]
[
.B \-x
]
[
.B \-yflip
]
[
.B \-y
]
[
.B \-xyflip
]
[
.B \-xy
]
[
.B \-zlow
.I lowpixvalue
]
[
.B \-l
.I lowpixvalue
]
[
.B \-zhigh
.I highpixvalue
]
[
.B \-h
.I highpixvalue
]
[
.B \-zoom
.I zoomvalue
]
[
.B \-z
.I zoomvalue
]
[
.B \-range
.I colormapsize
]
[
.B \-r
.I colormapsize
]
[
.B \-color
.I colormapselection
]
[
.B \-c
.I colormapselection
]
.B filnam
[
.B lowpixvalue
]
[
.B highpixvalue
]
.SH DESCRIPTION
.SH OPTIONS
.LP
Most options can be specified in two forms, a short form and a long form.
.LP
The file named \fBfilnam\fP
is displayed.  This name may possibly be a name meaniful only to some package,
like an \s-2AIPS\s0 file name (See Section on File names).
If no \fBfilnam\fP argument is specified then a name of BLANK512x512 is
assumed which is a dummy file which containing a 512 by 512 array of zeroed
pixels.
.LP
If the \fB-zlow\fP or (alternately) \fB-l\fP option
is not present then the next positional argument is interpreted as
\fBlowpixvalue\fP which has the same meaning as the \fB-zlow\fP option (see
below).
.LP
Likewise if the \fB-zhigh\fP or (alternately) \fB-h\fP option
is not present then the next positional argument is interpreted as
\fBhighpixvalue\fP which has the same meaning as the \fB-zhigh\fP option (see
below).
.\"
.IP \fB-debug\fP
.IP \fB-d\fP
Turn on the debugging flag.  This generates volumnious information on the
standard error output that is sometimes useful for debugging \fIds\fP.
.\"
.IP \fB-mono\fP
.IP \fB-m\fP
Force the display into monochrome mode, even if \fIds\fP is being run from
a color display.  In monochrome mode a dithered array of black or white pixels
is generated to represent a picture and all color facilities are disabled other
than switching between a positive or a negative (inverted) display.
.\"
.IP \fB-xflip\fP
.IP \fB-x\fP
Flip the picture horizontally. This changes the origin from its default
position on the left of the window to instead be on the right.
.\"
.IP \fB-yflip\fP
.IP \fB-y\fP
Flip the picture vertically. This changes the origin from its default position
on the top of the window to instead be on the bottom.
.\"
.IP \fB-xyflip\fP
.IP \fB-xy\fP
Interchange the vertical and horizontal axes.  The location of the origin on
the window is unchanged.
.\"
.IP "\fB-zlow\fP \fIlowpixvalue\fP"
.IP "\fB-l\fP \fIlowpixvalue\fP"
Set the low pixel value for scaling to \fIlowpixvalue\fP. Exponential form
(3.4.e-2 for example) is allowed.  \fIlowpixvalue\fP and \fIhighpixvalue\fP
together set the scaling between the displayed pixels and the pixels contained
in the original picture.  If \fBcolor\fP is \fIpositive\fP then
\fIlowpixvalue\fP corresponds to a black display and \fIhighpixvalue\fP
corresponds to a white display.  \fIhighpixvalue\fP can be lower than
\fIlowpixvalue\fP in which case a negative display will be obtained.
.\"
.IP "\fB-zhigh\fP \fIhighpixvalue\fP"
.IP "\fB-h\fP \fIhighpixvalue\fP"
Set the high pixel value for scaling to \fIhighpixvalue\fP.  Exponential form
(3.4.e-2 for example) is allowed.  \fIlowpixvalue\fP and \fIhighpixvalue\fP
together set the scaling between the displayed pixels and the pixels contained
in the original picture.  If \fBcolor\fP is \fIpositive\fP then
\fIlowpixvalue\fP corresponds to a black display and \fIhighpixvalue\fP
corresponds to a white display.  \fIhighpixvalue\fP can be lower than
\fIlowpixvalue\fP in which case a negative display will be obtained.
.\"
.IP "\fB-zoom\fP \fIzoomvalue\fP"
.IP "\fB-z\fP \fIzoomvalue\fP"
Set the initial value of the zoom factor of the display.  If the zoom factor
is equal to 1 then one pixel on the displayed picture corresponds to one
pixel on the display. A zoom factor of 2 produces a times 2 magnification with
one pixel in the picture interpreted as a 2 by 2 array of pixels on the
display. Magnification is currently done by pixel duplication.
Negative zoom factors are interpreted as being the reciprical of their
abolute value.  For example a zoom factor of -2 produces a times two reduction
in size with a 2 by 2 array of pixels in the picture corresponding to 1 pixel
on the display. Reduction is currently done by pixel averaging with the pixel
on the display being the average of the corresponding pixels in the picture.
Currently the initial zoom factor must be a value for which the picture will
fit on the screen. Otherwise the zoom value will be decreased to such a value.
The default initial zoom value is 1.
.\"
.IP "\fB-range\fP \fIcolormapsize\fP"
.IP "\fB-r\fP \fIcolormapsize\fP"
Set the number of colors used to display the picture on a color or grey/level
display.  This must be a power of two and be less than or equal to 256.  The
default value is 128.  A \fIcolormapsize\fP of 1 is equivalent to the
\fB-mono\fP option, setting the display into monochrome mode.
.\"
.IP "\fB-color\fP \fIcolormapselection\fP"
.IP "\fB-c\fP \fIcolormapselection\fP"
.\"
.LP
\fIds\fP also accepts all of the generic tool arguments discussed in
\fIsuntools\fP(1).
.SH SUBWINDOWS
.LP
\fIds\fP consists of up to four subwindows.
.IP "control panel" 16
.\"
.IP "picture canvas" 16
.\"
.IP "ramp canvas" 16
.\"
.IP "blanking window" 16
.\"
.SH "CONTROL PANEL ITEMS"
.\"
.IP "File Name Text Item"
This item shows the name of the file currently displayed and can be used to
select another file to display.
To select this item: point the cursor at it and press the
left mouse button.
A triangular blinking cursor will appear
to the right of the file name.
To prevent accidental changing of the file file, the
triangular cursor will not appear in this item unless it has been
explicitly selected. With the mouse cursor then still inside the
control panel the keyboard can be used to change the file name. The old file
name can be cleared out with your line kill character (usually control U or
control X) or your delete character (usually Delete or backspace).
The file name entered here may possibly be a name meaniful only to some
package, like an \s-2AIPS\s0 file name (See Section on File names).
When carriage return is pressed
.B ds
will be notified to load the new file. If the new file can not be found
the old file will continue to be displayed.  When the new file is being loaded
the cursor will be replaced by a cursor shaped like an hourglass, to indicate
that
.B ds
is busy and the message "loading file..." will appear in the message item.
.\"
.IP "Z Low Text Item"
.IP "Z High Text Item"
These are the two Z boundary items.
The value of these two items determine the
mapping between pixel values in the original picture file and the brightnesses
and colors shown on the screen.
Ordinarily value of the Z Low Text Item corresponds to pixels that are
displayed as black and the value of the Z High Text Item corresponds to pixels
thatat are displayed as white, yielding a positive display
(currently the only exception to this statement is
when the Colormap Cycle Item has a value of Negative. Then the sense of
the statement is reversed: Z Low pixels are white and Z High Pixels are black).
However a Z Low Text Item value higher that the
value of the Z High Text Item is legal, yielding a negative or color inverted
display.
Exponential form (3.4.e-2 for example) is allowed and is useful if the picture
contains floating point pixels.
The keyboard can be used to change the value of the item if the mouse cursor is
inside the control panel and a triangular blinking cursor is present in the
item.  If the triangular cursor is not there already, it can be brought to the
item by selecting the item (pointing to it with the mouse cursor and pressing
the left mouse button).
The old value can be cleared with your line kill character (usually control
U or control X) or your delete character (usually Delete or backspace).
Then the new value can be entered.
Pressing carriage return
notifies
.B ds
of your new value and sends the triangular cursor over to the other Z boundary
Text item (the triangular cursor can only be moved to the File Name Text Item
by selecting it explicitly).
The new value will not be applied to the display of
the picture until the Repaint Button Item on the control panel is invoked.
.\"
.IP "Zoom Factor Slider Item"
This item controls the value of the zoom factor. Position the cursor inside
the rectangular box (the "slider rectangle") and press the left mouse button.
Move the cursor to the left or the right until the zoom factor (given in square
brackets to the left of the slider rectangle) has the value desired.  Then
release the left mouse button and the display will be readjusted to the new
zoom factor.  The operation can be aborted and the previous zoom factor
retained by moving the cursor out of the slider rectangle before releasing the
left mouse button.
.\"
.IP "Colormap Cycle Item"
.\"
.IP "X flip toggle"
Point cursor at this item and press left mouse button to change the value of
the X flip toggle.
Changing the value of the X flip toggle flips the picture horizontally.
When the X flip toggle is off the origin is on the left side of the window.
When the X flip toggle is on the origin is on the right side of the window.
.\"
.IP "Y flip toggle"
Point cursor at this item and press left mouse button to change the value of
the Y flip toggle.
Changing the value of the Y flip toggle flips the picture vertically.
When the Y flip toggle is off the origin is on the top side of the window.
When the Y flip toggle is on the origin is on the bottom side of the window.
.\"
.IP "XY flip toggle"
Point cursor at this item and press left mouse button to change the value of
the XY flip toggle.
Changing the value of the XY flip toggle flips the picture on a diagonal line
through the picture, interchanging the horizontal and vertical axes.
The location of the origin on the window is unchanged.  The XY flip is applied
before the X and Y flips.  Thus an X flip will always interchange the left and
right hand sides of the display and a Y flip will always interchange the top
and bottom of the display.
When the XY flip toggle is off the pixel rows in the picture run horizontally.
When the XY flip toggle is on the pixel rows in the picture run vertically.
.\"
.IP "Blank Button"
Point cursor at this "button" and press left mouse button to invoke the screen
blanking feature intended for use with hard copy devices.  A window the color
of the background (usually black) appears covering the whole screen. In its
center is a window on the picture being displayed at the current pan and
zoom.  For more information see the subsection on "Blanking Feature".
.\"
.IP "Repaint Button"
Point cursor at this "button" and press left mouse button to repaint the
image in the picture canvas.  This is generally used after the brightness
scaling factors "Z Low" or "Z High" have been changed to show the image with
the new scaling.
.\"
.IP "Quit Button"
Point cursor at this "button" and press left mouse button to exit ds. A
confirmer window will appear.  Press the left mouse button a second time to
confirm the quit or press the right mouse button to cancel the quit.
.\"
.IP "Help Button"
Point cursor at this "button" and press left mouse button for help on the
various features in ds. At the moment all that happens is the message appears
in the message item "help not yet implemented". Pressing the help button will
eventually bring up popup menu of items for which help is requested.
.\"
.IP "Message Item"
.\"
.SH "COMMANDS IN PICTURE CANVAS"
.SH "BLANKING FEATURE"
.SH "FILE NAMES"
.PP
\fIDs\fP can display files that have several different formats.  Each format
has associated with it two routines, an examine and a load routine.  Each
examine routine is called in turn to look at the file.  The examine routine
can return 3 different statuses.  If it returns NOT_OURS it goes on to the
next examine routine in the list.  If it returns BAD the search is aborted. If
it returns OK the load routine for that format is called.  If BAD status was
returned by the examine routine or if the load routine returns bad status
an attempt is made to reload the previous picture.
A status of NOT_OURS returned by the examine routine
may mean that the file didn't exist under that name, as for some of the formats
it may be stored under a modified version of that name.  For example: iraf
files may have a .imh tacked onto the end of the name or AIPS files may require
a search of an AIPS catalog to translate the name into the long funny
name where the image is actually stored. A name that begins with the characters
BLANK is not expected to be on disk at all but results in a image containing
zeroed pixels.  
.PP
Once a picture is loaded then it is present in memory and no further reference
to it is made by the \fIds\fP program.  The file formats recognized by the
load routine are:
.RS
.IP "1) Iraf:"
These files are generated by the IRAF image processing package.
If filnam is a string of characters (uppercase/lowercase
significant) then filnam or filnam.imh (with a .imh extension)
are recognized.  The accompanying .pix file (which must exist)
is located by reading the header in the .imh file.  Hostname
prefixes in the .pix file name are currently ignored.
.IP "2) Fits:"
These files are disk images in the FITS tape transmission
standard. The file name recognized is the same as the name of
the disk file used to store the images.
.IP "3) Sun Raster:"
These files are in the standard form for picture storage used
by Sun Microsystems software. The file name recognized is the
same as the name of the disk file used to store the images.
Either 8 bit or 1 bit monochrome images are recognized.  Monochrome images
result in pixels that either contain 0 or 255.
.IP "4) Aips:"
These files are generated by the AIPS image processing package.
The name can either be the name of the actual image head (of
the form CBBxxxxx.xxx;1 ) or the AIPS catalog name.
The AIPS catalog name  is constructed in the form
aipsname`class`seqno`aipsid .  A backquote separates the name
sections.  A doubled backquote is interpreted as a single
backquote in a name section.  The "aipsid" name section is
optional and if omitted it defaults to the user's uid number
from the /etc/passwd file.  If the "seqno" section is omitted
then the picture with an otherwise matching name that has with
the highest seqno will be used.  A name or class of * will
match any name or class. An omitted class is equivalent to a
class of *.
the search for an aips catalog (file name beginning with the
letters "CAA") begins in the current directory.  If that fails
the search continues at directories pointed to by environment
variables DA01, DA01, DA02, ... DA0n.  The search terminates
when the file is found or when an environment variable does
not exist.
.IP "5) NMR picture:"
These files are used to store NMR images from the Laurie
Imaging center.  The file name recognized is the same as the
name of the disk file used to store the images.  The file
consists of 256 x 256 2 byte pixels offset by a 2048 byte
header. Pixels are in byte-swapped order, i.e. the less
significant byte is first.
.IP "6) blank file:"
Generates a picture filled with zero pixels. This picture can
be used as a simple default or as a blank background for
graphics.  The file name consists of the letters BLANK followed
optionally by a single number or by two numbers separated by
the letter "x". The name BLANK generates a 512 by 512 pixel
blank picture.  The name BLANKn where "n" is a number generates
an n by n pixel blank picture.  The name BLANKnxm where n and
m are numbers generates an n by m pixel blank picture.
.RE
.SH BUGS
.LP
