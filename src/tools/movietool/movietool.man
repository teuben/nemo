.TH movietool L "20 November 1990"
.SH NAME
Movietool \- a display tool for sequences of
.I Sun rasterfiles
providing "video animation" in a Sun window.
.SH SYNOPSIS
.I movietool
.I [SunView Generic Args]
.I [-encoded]
.I [-zoom factor] 
.br
.I 	[-sub x y w h]
.I [-ident]
.I [-XOR]
.I [-background]
.br
.I 	[-mono_panel]
.I [-audio]
.I [-CDROM]
.I  rasterfile[s]
.SH DESCRIPTION
.I Movietool
displays a sequence of color or monochrome raster files (see "man 5 rasterfile")
in real-time or near-real-time.
On a Sun 3/50M with 4 MB of RAM, monochrome images of 600 by 800 pixels
may be played at about 16 frames per second.
.I Movietool
also allows previewing single or multiple raster files,
even if you're not interested in animation.
.SH ARGUMENTS
After the SunView generic arguments (see generic_args(1)),
Movietool takes the following arguments:
.IP "-e[ncoded]"
With this option,
the rasterfiles' data are stored as byte-encoded data
(provided the rasterfile format was
.I RT_BYTE_ENCODED).
This could typically reduce memory consumption by as much as
.I 50-80 %,
providing much smoother animation, or longer sequences,
with a given amount of available RAM memory.
The decoding process of course requires more of the CPU
than a simple memory copy operation,
but when you are lacking RAM memory this tradeoff is a big win.
.IP "-z[oom] factor"
does raster replication to magnify your rasterfiles.
The
.I factor
must be between 2 and 10,
just for safety.
This option also works with the -encoded option,
where rasters are zoomed during the decoding process.
.IP "-s[ub] x y w h"
picks a 
.I subregion
of your rasterfiles starting at pixel (x,y) of width=w and height=h.
This option has not been implemented for the -encoded option yet,
so Movietool will refuse the -s together with the -e flag.
.IP "-i[dent]"
Use 
.I identical colormaps
for all (color) frames,
thereby reducing memory and time consumption.
The color map of the first frame is used throughout.
.IP "-X[OR]"
For a sequence of rasterfiles, store the
.I XOR
of successive rasters.
When playing a sequence, the updating takes place by means of XOR-ing
what is on the screen.
Since most of the XOR'ed pixels are presumably zero,
this provides the possibility for very fast updating of the screen.
However, the simple frame buffers
(I have tried /dev/bwtwo1, /dev/cgfour0 and /dev/cgtwo0 with GP2)
do not appear to do XOR's efficiently.
This option does not work with the -encoded option.
.IP "-b[ackground]"
The canvas background will be painted the color of the upper left,
or first, pixel of the first image,
rather than the default SunView color.
The background is the area of the canvas not covered by the
raster image currently displayed.
.IP "-m[ono_panel]"
With this flag,
the Movietool control panel will use the
default monochrome colors of your screen.
This may be useful for images whose colormaps make the control buttons
invisible (or at least hard to distinguish).
.IP "-a[udio]"
Play an audio file on the system's speaker.
The audio player filename can only be changed by recompiling Movietool.
.IP "-C[DROM]"
Pop up a slightly modified version of the Sun CD-ROM "cdplayer" program,
intended for making Movietool a "multimedia" tool.
The audio CD-ROM will start and pause along with the movie.
The CDROM player filename can only be changed by recompiling Movietool.
.IP "rasterfile[s]"
argument[s] must be in the
.I rasterfile(5)
format.
Files not conforming to this format are just ignored.
.SH "CONTROL PANEL"
There are
.I Play
(animate the sequence),
.I Step
(single step through the movie),
.I Clear
(clear display canvas), and
.I Eject
(finish up) buttons on the control panel.
If only a single rasterfile is to be displayed,
the control panel consists of only the Eject button.
.LP
While the movie is playing,
the mouse cursor shows an hour-glass in the control panel,
and the movie may only be interrupted by pressing the
.I STOP
button located at the top left-hand corner of the keyboard
(\fIL1\fP on many keyboards).
Since the display-canvas is non-retained and attempts to use the overlay
plane on the framebuffer (if available),
damage to the window is not repaired while the movie is playing.
To repair, interrupt the movie and restart it.
.LP
If the
.I -XOR
option has been selected,
the first line of the panel also features a switch between
XOR-ing and a normal copy,
as well as the possibility to view the XOR'ed raster.
.LP
Second line features a 
.I Frame
slider indicating framenumber as the movie proceeds,
which in addition may be used to position at a particular frame.
The file name corresponding to the current frame is displayed
(except when Playing).
.LP
Third line shows,
and allows selection of,
play 
.I Direction: forward
or
.I backward.
.PP
Fourth line controls repetition of the movie:
The
.I Repeat
button allows repetition as indicated by Direction.
If
.I Auto-reverse
is flipped on,
the movie turns around at either end.
In this case, the Repeat button is disabled.
.PP
If you toggle the
.I Clipping
button "Off",
an attempt will be made to decode
the rasterfiles directly onto the screen pixrect,
when pressing the "Play" button.
This is mostly useful when playing RT_BYTE_ENCODED raster images,
where you avoid the internal copying to a memory pixrect
before displaying on the screen.
It is your own
.I responsibility
to enlarge Movietool sufficiently so that the Clipping Off option does
not damage other windows !
If you disable clipping,
it is recommended to zoom up Movietool by using the Resize=>FullScreen popup.
However, even if you disable clipping,
Movietool may not always be able to provide faster display,
and clipping may occur correctly anyway.
Also, when single-stepping through rasterfiles,
clipping will be performed because speed is not an issue.
.LP
Finally, the
.I Frames/sec
slider lets you select the speed of the movie.
By choosing the maximum value, you get whatever maximum performance
your system will deliver.
At lower values delays will be imposed,
but if the actual speed is lower than the value selected,
the slider will flick to the current number of frames per second.
Note that the system's time-of-day clock is so coarse that at above 16
frames/sec it is hard to provide proper timing.
.SH PICTURE DISPLAY
The raster image will always be centered on the display canvas,
provided the canvas is large enough.
In case your image is larger than the available canvas,
the scroll bars allow you to select part of the image.
.SH CAPACITY
On a 3/60C, the rate at which rasters can be displayed on the screen
seems to be about 6-8 Mega-pixels per second,
which may easily be 25-50 frames per second for small images.
Monochrome images generally play 8 times faster than color images of
the same size.
Also, if the Movietool is shrunk to display just the interesting region of
the rasterfiles (or when using the
.I -sub
argument),
the display may become significantly faster.
The most important limitation to long image sequences is the amount
of RAM available to hold the memory-pixrects containing the images.
When paging to the disk sets in, 
you will notice immediately !
.SH AUTHOR
Movietool was written by:
.nf
Ole H. Nielsen
Lab of Applied Physics, Bygn. 307
Technical University of Denmark, DK-2800 Lyngby, Denmark.
E-mail: ohnielse@ltf.dth.dk
.fi
.LP
You may redistribute Movietool freely,
provided you maintain the Copyright notice of the author,
and provided you don't attempt to sell Movietool for profit.
.LP
The development of
Movietool was inspired by a similar tool called VCR by 
Hsuan Chang (hsc@vuse.vanderbilt.edu),
but Movietool has rather improved functionality and performance.
