.TH xrastool n "11 May 1993" "Ver: 1.0"
.SH NAME
\fIxrastool\fP \- animate rasterfiles in an \fIXView\fR window
.SH SYNOPSIS
\fIxrastool\fR [generic_options] [other_options] [filename...]
.SH DESCRIPTION
\fIxrastool\fR is an animation package written in XView.  Currently
only the 8-bit Sun Rasterfile (RLE or standard) image format is
supported, but many other formats can be converted to 8-bit rasterfiles
using the \fIPBMplus\fR image utilities or \fIxv\fR.  \fIxrastool\fR
provides the user with an \fIXView\fR panel interface for convenient
access to animation, display, and colormap functions.  There are a
number of options to help maximize the display rate.  On an unloaded
Sparc IPX running \fItwm\fR under X11R5, a display rate of 50
frames/sec can be achieved with 400x400 8-bit images.  These speeds are
obtained by first loading the rasterfiles into client-resident Ximages
then transferring them to server-resident pixmaps.  Consequently, the
maximum number of frames displayable (and hence the duration of the
animation) is limited only by the core memory available to the server.
.SH LAYOUT
\fIxrastool\fP consists of two panels, a canvas frame for displaying
images, and a pop-up info window.  When invoked without any options or
filenames specified, \fIxrastool\fR displays the Main Panel in its
default state.  If any images were specified and at least one was
loaded successfully, the first will be displayed in the Canvas frame,
also in the default state.  The cursor is moved automatically to point
to the \fBMore...\fP button (see below) on the Main Panel for
convenience.  This also ensures that the colormap is displayed
correctly for the first image by forcing the window manager to load
\fIxrastool\fP's colormap segment.  The Sub Panel and Info window can
be called up from the Main Panel.  The following two sections describe
the control panels in detail.
.SH Section 1: The Main Panel
The Main Panel contains all the controls needed to prepare and run an
animation.  Features include image loading facilities, sizing controls,
timer controls, and cycling buttons.  All the controls are described in
detail in what follows.  Note that any equivalent command-line options
are given in parentheses beside the control name (see Section 4).
.IP "\fBImage #\fP (-\fIstart\fP)"
The current image number is shown here and may be modified, either
using the text field or the incremental buttons.  Only numbers in the
range 1 to N, where N is the number of images loaded, are accepted as
input, although the value N + 1 may be displayed prior to loading a new
image (see \fBNew\fP button).
.IP "\fBNew\fP"
Click on this button to clear the image filename field and set the
image number to N + 1 in preparation for a new image.  Any displayed
image will be replaced with the current background.
.IP "\fBLoad\fP"
Clicking this button will cause \fIxrastool\fR to load the image shown
in the \fBImage\fP field.  The image will be assigned the image number
shown in the \fBImage #\fP field, so \fBLoad\fP can be used to replace
existing images if desired.  Note that if there is insufficient server
memory to store the image, the functions XCreateImage() or
XCreatePixmap() function may abort with an error.
.IP "\fBLogical Functions\fP"
The buttons at the top right of the main panel are used to change the
logical function of the main Graphics Context (GC).  The logical
function is applied to the next and all subsequent images to be
displayed on the Canvas.  Use \fBRefresh\fP to apply the function to
the current image.  Choices are:
.RS
.IP "\fB+\fP (-\fIlogical\fP SRC)"
This applies the GXcopy function to the images so that any existing
image (destination) is entirely replaced with each new image (source).
This is the default behavior.
.IP "\fB|\fP (-\fIlogical\fP OR)"
This applies GXor between the source and destination images,
effectively merging them together (color behavior may be surprising).
.IP "\fB&\fP (-\fIlogical\fP AND)"
Applies GXand between the source and destination.
.IP "\fBx\fP (-\fIlogical\fP XOR)"
Applies GXxor.  Two invocations of GXxor on the same image returns the
original image.
.IP "\fBe\fP (-\fIlogical\fP ERA)"
Applies GXandInverted.
.IP "\fB-\fP (-\fIlogical\fP NEG)"
Applies GXcopyInverted, giving the negative of the image.
.RE
.IP "\fBImage\fP"
The filename of the current image is displayed in the text field.
Typing a new filename and pressing RETURN or clicking on \fBLoad\fP
instructs \fIxrastool\fP to replace the current image with the contents
of the specified rasterfile.  A new image can be appended to the
existing set by first clicking on \fBNew\fP, then supplying the new
filename.
.IP "\fBW\fP (-\fIw\fP (width))"
This is the image width field.  Changing this number either using the
text field or the incremental buttons causes the current \fBSizing\fP
option to change to \fBFixed\fP.  Any existing image will have its frame
adjusted to accommodate the new size.
.IP "\fBH\fp (-\fIh\fP (height))"
The height of the Canvas frame can be adjusted with this control,
in the same way as the width (\fBW\fP).
.IP "\fBSizing\fP"
There are three sizing options that can be used for displaying
images.
.RS
.IP "\fBAuto\fP (-\fIauto\fP)"
With this option, images are displayed up to a specified maximum size,
by default 640x480.  The maximum size can be changed by first setting
\fBW\fp and \fBH\fp, then clicking on \fBAuto\fP.  The display frame
will shrink to accommodate images with dimensions smaller than the
maximum dimensions.  For larger images, the scrollbars can be used to
view hidden portions of the image.
.IP "\fBFull\fP (-\fIfull\fP)"
With this option, the entire image is displayed as long as it can be
accommodated on the default screen.  The \fBW\fP and \fBH\fP fields
show the displayed image dimensions.
.IP "\fBFixed\fP (-\fIfixed\fP)"
This option fixes the image size to the current values displayed in the
\fBW\fP and \fBH\fP fields.  Fixed canvases are the fastest to draw,
so this is the preferred option when maximizing display rates.  Hidden
portions of images can still be viewed with the scrollbars.  Images that
are smaller than the fixed size are padded with the current background
(see the \fBBackdrops\fP item in the Sub Panel) and may be centered in
the frame by choosing \fBCentering\fP from the Sub Panel.
.RE
.IP "\fBSet\fP"
This button is used to set the \fBW\fP and \fBH\fP fields to the
currently displayed size.  This means the frame size can be adjusted by
first using the mouse (e.g. dragging the resize corners if the window
manager is olwm) and then clicking on \fBSet\fP.  This button also
causes \fBFixed\fP to be selected automatically.
.IP "\fBTimer Scale\fP"
The MINIMUM delay between displaying one frame and the next during animation
is set by the \fBTimer Scale\fP and \fBTimer Slider\fP.  The former is a
coarse adjustment and the latter a fine adjustment.  It is possible to
select delay ranges between 1 microsecond and 999 seconds (about 15 min),
though 20 milliseconds (i.e. 50 frames/sec) is about the best achievable
on a Sparc IPX.  The practical lower limit is set by hardware CPU and I/O
rating.  Note that at the fastest speeds, the \fIxrastool\fP controls may
become very sluggish -- be patient! Ideally, all display options should
be set PRIOR to animating at top speed.
.RS
.IP "\fBus\fP (-\fItimerscale\fP us)"
This control selects the microseconds timer scale.
.IP "\fBms\fP (-\fItimerscale\fP ms)"
Use this for the milliseconds timer scale, which is the practical range
for animations.
.IP "\fBsec\fP (-\fItimerscale\fP sec)"
The seconds scale can be used for slideshows.
.RE
.IP "\fBTimer Slider\fP (-\fItimervalue\fP (timer))"
This slider is used for fine adjustments of the delay timer for
animations.  The range is 1 to 999 in the units selected on the
\fBTimer Scale\fP.  The value can be set using the text field or by
dragging the slider control.  Unit increments or decrements can be
performed by clicking either side of the control.
.IP "\fbFast\fP (-\fIfast\fP)"
Fast mode can be selected by pressing this button.  This sets
\fBFixed\fP sizing, \fBus\fP timer scale, unit timer value, and many of
the options on the Sub Panel.  In this mode, \fIxrastool\fP will
display frames at the fastest possible rate.  Any default fast setting
can be overridden by adjusting the appropriate control in the usual
way.
.IP "\fbDirection\fP (-\fIfwd\fP -\fIrev\fP)"
Cycling direction can be toggled between forward (FWD) and reverse
(REV) using this control.  \fBBlink\fP mode and the \fBLoop Back\fP and
\fBOne Way\fP cycling options will take over direction control when
invoked.
.IP "\fbStep\fP"
Clicking on this button will display either the next or previous image
depending on the current \fbDirection\fP.  The step function is called
internally when animating.
.IP "\fBBlink\fP (+|-\fIblink\fP)"
This button toggles blinking on and off.  The current image and the next
or previous image are blinked at the current display rate.
.IP "\fBCycle\fP (+|-\fIcycle\fP)"
This button displays the stored images in sequence, subject to the
current Direction, \fBCycling Option\fP, and display rate.  The
animation can be stopped by clicking on \fBCycle\fP again.
.IP "\fBCycling Option\fP"
There are three choices of cycling option:
.RS
.IP "\fBLoop\fP (-\fIloop\fP)"
This is the default.  When the end of the current sequence of images
is reached, cycling continues from the beginning ad infinitum.
.IP "\fBLoop Back\fP (-\fIloopback\fP)"
When the end of the current sequence of images is reached, the
current \fBDirection\fP is reversed and cycling continues.
.IP "\fBOne Way\fP (-\fIoneway\fP)"
When the end of the current sequence of images is reached, cycling
is terminated.  Clicking on \fBCycle\fP again will repeat the sequence,
starting with the first or last image depending on the current
\fBDirection\fP.
.RE
.IP "\fBRefresh\fP"
Clicking here redisplays the current image.  This can be used on an
image to apply a \fBLogical Function\fP to itself.  It is also useful if
the window manager has corrupted the image display for whatever
reason.
.IP "\fBMore...\fP (+|-\fIsubpanel\fP)"
This pops up the Sub Panel.
.IP "\fBInfo...\fP"
Use this button for info relating to the author and copyright.
.IP "\fBQuit\fP"
Clicking here terminates the \fIxrastool\fP session.
.SH Section 2: The Sub Panel
When selected using \fBMore...\fP on the Main Panel or by including the
-\fIsubpanel\fP option on the command line, the Sub Panel is displayed.
The user can perform various fine-tuning operations, including colormap
manipulation, using the Sub Panel.
.IP "\fBNo Scrolling\fP (+|-\fIscrolling\fP)"
Selecting this removes the scrollbars from the image.  Deselecting
replaces them.
.IP "\fBNo Resizing\fP (+|-\fIsizing\fP)"
This option overrides all sizing requests.  The frame size can still be
adjusted manually using the mouse, but a new size cannot be \fBSet\fP
until \fBNo Resizing\fP is deselected.
.IP "\fBNo Moving\fP (+|-\fImoving\fP)"
When an image is to be displayed, \fIxrastool\fP attempts to reposition
the frame so that no part of the image will be off-screen.  Selecting
\fBNo Moving\fP disables this behavior.
.IP "\fBNo Scaling\fP (+|-\fIscaling\fP)"
Selecting this option disables all color \fBScaling\fP effects.
.IP "\fBNo Updates\fP (+|-\fIupdates\fP)"
Normally when cycling, the \fBImage #\fP and \fBImage\fP fields, as
well the size fields and Main Panel footer, are updated between images.
Selecting \fBNo Updates\fP disables this behavior when cycling.
Updates when using \fBStep\fP manually cannot be suppressed.
.IP "\fBNo Backdrops\fP (+|-\fIbackdrops\fP)"
By default when cycling, images smaller than a \fBFixed\fP frame will
have undefined regions painted with the current backdrop.  \fBNo
Backdrops\fP disables this feature when cycling.
.IP "\fBLock Colors\fP (+|-\fIlockcolors\fP)"
Selecting this locks the current image colormap so that all
subsequently displayed images use the locked colormap rather than their
own.  This feature results in interesting effects between images with
very different colormaps, but significantly improves the display rate
for images with identical or nearly identical colormaps.  If the command
line option -\fIonecmap\fP is used, the \fBLock Colors\fP control is
made redundant.
.IP "\fBCentering\fP (+|-\fIcentering\fP)"
When selected, and if the \fBSizing\fP is \fBFixed\fP, the current image
and all subsequently displayed images will be centered in the display
frame, assuming the frame is larger than the image.  Any gaps will be
filled with the current backdrop.
.IP "\fBLive Cursor\fP (+|-\fIlivecursor\fP)"
With a live cursor, the mouse can be used to examine the current color
\fBScaling\fP of an image.  Simply moving the mouse across the image
will display the cursor coordinates with respect to the top left corner
of the image, the value of the pixel at that location, the original
color RGB triad to which the pixel value corresponded (after
compression -- see below), and the current color at that location.  For
example:

.RS
.IP "(122,36) = 53 (96,96,96) --> (255,0,0)"
.RE
.IP
shows that the pixel at x = 122, y = 36 on the image has a value of 53,
which originally corresponded to the gray-scale value (96,96,96) but
has since been mapped into the color red (255,0,0).  The coordinates
and mapping information are displayed in the footer of the Canvas frame.
.IP "\fBBackdrops\fP (-\fIbackdrop\fP (index))"
Various bitmap backdrops can be selected for filling gaps left when an
image is smaller than its frame.  The patterns available are displayed
on the choice buttons (blank means white).
.IP "\fBScaling\fP"
There are six color scaling options available with \fIxrastool\fP.  The
RGB color sliders below the options allow fine tuning of the scaling.
If the code was compiled with LIVE_SLIDERS, changing one color slider
may instantaneously affect the other sliders, depending on the option.
Otherwise the other sliders will be updated after the slider control of
the first is released.  The \fBx\fP and \fB/\fP buttons to the right of
each slider are only used with the \fBColor\fP option.  Each scaling
option has a default slider setting built in which is used when the
option is first selected.
.RS
.IP "\fBBW\fP (-\fIbw\fP)"
The slider controls a fractional cutoff: colormap entries with a total
intensity less than (slider/100) * 255 are mapped to black while the
remainder are mapped to white.  The default setting is 50.  Any slider
can be used to change the setting.
.IP "\fBGray\fP (-\fIgray\fP)"
The sliders control the fractional contributions of each color in the
image to a gray scale map R_new = (R/100) * R_old + (G/100) * G_old +
(B/100) * B_old, B_new = G_new = R_new.  The sum of the slider values R
+ G + B is fixed at 100 unless any two sliders have a zero setting.  The
default setting is R = 30, G = 59, B = 11 which is the map used by
\fIppmtopgm\fP in the \fIPBMplus\fP package.
.IP "\fBColor\fP (-\fIcolor\fP)"
This is the default scaling option.  With this option, the sliders move
independently, controlling the fractional intensity of their respective
colors in the image.  Selecting \fBx\fP or \fB/\fP for each slider
determines whether the color contribution is diminished (former) or
enhanced (latter).
.IP "\fBCutoff\fP (-\fIcutoff\fP)"
With this option, the sliders control the pixel value above
which no colors are displayed.  Hence the first color in the colormap
can be displayed by setting a slider to 1.  All colors are displayed
for a setting of 100 (the default).  In a future release, the colormap
will be sorted in order of usage, so that the most or least "important"
colors can be screened out using the \fBCutoff\fP option.
.IP "\fBContour\fP (-\fIcontour\fP)"
This option maps the image colors by intensity into 10 repeating
colors, with the total number of contours controlled by the sliders
(from 1 to 100, default 25).  The contour colors have been selected in
order of intensity so that the brightest colors in the image get the
brightest contours.
.IP "\fBRandom\fP (-\fIrandom\fP)"
This option simply maps the image colors into random values.  The values
can be changed with the sliders.
.SH Section 3: Keyboard Accelerators
Limited keyboard support is provided on the Canvas window.  Pressing the
spacebar while the cursor is on the Canvas is equivalent to pushing the
\fBStep\fP button.  Pressing RETURN performs a \fBSet\fp operation.
Press 'c' to toggle cycling, and 'q' to quit the program.
.SH Section 4: Command Line Options
Most of the panel controls have command-line equivalents so that an
\fIxrastool\fP session can be started with various options already
set.  The usage is:
.IP "\fIxrastool\fP"
[-fast] [+|-onecmap] [+|-subpanel] [-start (index)]
.ti +0.5
[-private|-hist|-pair] [-logical (function)]
.ti +0.5
[-w (width)] [-h (height)] [-auto|-full|-fixed]
.ti +0.5
[-timerscale us|ms|sec] [-timervalue (timer)]
.ti +0.5
[-fwd|-rev] [-blink|-cycle] [-loop|-loopback|-oneway]
.ti +0.5
[+|-scrolling] [+|-resizing] [+|-moving] [+|-scaling]
.ti +0.5
[+|-updates] [+|-backdrops] [+|-lockcolors]
.ti +0.5
[+|-centering] [+|-livecursor] [-backdrop (index)]
.ti +0.5
[-bw|-gray|-color|-cutoff|-contour|-random]
.ti +0.5
[-red (value)] [-green (value)] [-blue (value)]
.ti +0.5
[filename...].

Most of these options have been described in the previous sections.
Note that the "+" prefix is used to override defaults or previous
values.  For example, -\fIfast\fP invokes many of the Sub Panel options;
if it was desired to retain scrollbars, then +\fIscrolling\fP should be
placed after -\fIfast\fP on the command line.  Options appearing last
take precedence over those appearing first.  The options that have not
yet been described are:
.IP "-\fIonecmap\fP"
Instructs \fIxrastool\fP to store only the first image's colormap and
automatically \fBLock Colors\fP to it.  This speeds up the image loading
considerably as it bypasses any colormap compression (see below).  Of
course, the colormaps of all the images loaded should be similar, or,
ideally, the same.  The -\fIfast\fP option automatically invokes
-\fIonecmap\fP.
.IP "-\fIstart\fP (index)"
Instructs \fIxrastool\fP to display image (index) first.  If image
(index) does not exist, the first image will be displayed.
.IP "-\fIprivate\fP"
Instructs \fIxrastool\fP to use a private colormap for each image
so that no attempt at colormap compression is made.  Window colors,
including those of the \fIxrastool\fP control panels, may be severely
affected when viewing these images.
.IP "-\fIhist\fP"
Constructs a histogram of color usage for each image (or just the
first in the case of -\fIonecmap\fP), and eliminates those colors
that are least used until the image colormap can be melded with
the existing USED portion of the default display colormap.  Currently
\fIxrastool\fP reserves 30 colors for this purpose, so at most 226
colors from any image can be displayed.  If more colors are added to
the window session after xrastool has started up, they may be
"switched out" when viewing the images.  (NOTE: under twm, the panels
and canvas "black out" when they no longer have the keyboard focus;
under olwm, all colors can be seen at all times).
.IP "-\fIpair\fP"
Determines which colors are most alike in pair-wise fashion for each
image until 30 colors are eliminated.  This method is much slower than
-\fIhist\fP, but it has the advantage that it does not rely on the same
colors being used in EACH image for a multi-image animation (under
-\fIhist\fP it is possible that a color in one image may not be used in
a second -- even though they have the same colormaps -- and hence be
eliminated from the latter and not the former, resulting in noticeable
flickering when animating).  The -\fIpair\fP option is recommended when
using -\fIonecmap\fP.
.IP "-\fIred\fP (value) -\fIgreen\fP (value) -\fIblue\fP (value)"
These options allow the user to override the initial default color
\fBScaling\fP settings.  For \fBGray\fP scale, the images may be "burned
out": it is up to the user to ensure R + G + B = 100 if desired.
Changing the color sliders, however, will restore the constraint.  For
\fBColor\fP, the \fB+\fP and \fB-\fP choices can be selected by
supplying positive or negative color values, respectively.
.SH HINTS
For best results, never exceed the core memory of the server when
loading images.  Swapping will make the animator painfully slow.  Even
with more careful use of memory, \fIxrastool\fB may need to perform one
or even two cycles through the images before the animation becomes
smooth.  Once everything is up to speed, the animator controls should
respond better as well.  Also, to load animations that consist of a
large number of files, it is best to use the Unix wildcard facility
(e.g. \fIxrastool\fP -fast movie*.ras), so be sure that the image
sequence is the same chronologically as alphabetically.
.SH HISTORY
\fIxrastool\fB replaces its Sunview predecessor \fIrastool\fB, which
was never officially released.  These tools were designed to aid in
visualization of gravitational N-body problems carried out at the
Institute of Astronomy, Cambridge in partial fulfillment of the
author's Ph.D. thesis.
.SH COPYRIGHT
\fIxrastool\fP is Copyright (C) 1993 Derek C. Richardson under the
terms and conditions of the GNU General Public License (included
with the source distribution).
.SH AUTHOR
Derek C. Richardson, Institute of Astronomy, Cambridge U.K.
.SH BUGS
All bug reports, comments, praise, and criticism should be e-mailed to
dcr@mail.ast.cam.ac.uk.  Enjoy!
