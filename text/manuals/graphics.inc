%% This file is to be included by latex in nemo.tex
%
% Chapter: Graphics and Image Display
%
\myfile{graphics.inc} \mylabel{c:graphics}

NEMO programs also need to display their data of course.
We shall make
a distinction between {\it graphics} and {\it image} data.  A simple but
flexible {\it graphics} interface has been defined in NEMO and is used
extensively in programs.  
% Many plots used to get data into publication 
% quality format are still obtained by means of outside NEMO programs.
To display {\it image} data we currently rely on borrowware. Often
images need to be copied to a FITS file.


\section{The YAPP graphics interface}

\mylabel{s:yapp}
The programs in NEMO which use graphics are rather simple and allow no
interactive processing, except perhaps for a simple 'hit-the-return-key'
or 'push-a-mouse-button' between successive plots or actions.  A very
simple interface was defined ({\bf yapp}, Yet Another Plotting Package)
with basic plot functions.  \index{yapp}
There are currently a few yapp-interfaces
available, such as a postscript-device, suntools/windows, mongo and
pgplot.  If your output device is not supported by the ones available
in the current yapp directory
({\tt \$NEMO/src/kernel/yapp}), you have to write a new one!  A
reasonably experienced programmer writes a functional yapp-interface in
an afternoon, or so.

Although this method results in a flexible graphics interface, a
program can currently only be linked with one yapp-interface.  This will in
general result in the existence of at least two versions of the same
program, each for another graphics output device.  We use the 
convention that the ones for a
postscript printer have a {\tt \_ps} appended to their original name: the 
program which has the original name is the one whose display is the current
screen,
Hence we may see program names such as {\tt snapplot} (general Sun screen
within suntools), {\tt snapplot\_ps} (postscript), \index{snapplot, families}
{\tt snapplot\_cg} (color Sun screen) and {\tt snapplot\_sv}
(variable size sunview pixwindow for storing small images for movies).
Again: actual names may differ on your system.

If programs are linked with the multiplexing libraries
{\it yapp\_mongo}\index{mongo}
or {\it yapp\_pgplot}\index{pgplot}
interface, several device drivers are transparently present through
mongo, and the hidden system keyword {\tt yapp=} is then used to select
a device (a default can be set by using the {\bf YAPP} environment
variable).  \index{yapp=, system keyword} \index{YAPP, environment} See
also Appendix~\ref{a:iface}. 

\section{General Graphics Display}

Another convenient way to present data in graphical form is by using
the table format. We have already encountered the {\it tables} created by
many NEMO programs. These tables can be used by NEMO programs
such as {\it tabplot(1)}, {\it tabhist(1)}, and other packages
such as {\it mongo(1L)}, {\it sm(1L)} and 
{\it gnuplot(1L)},
\index{mongo} but \index{sm} even by completely foreign
packages such as Lotus 1-2-3.  \index{Lotus 1-2-3} 
Binary tables need to be converted to ASCII of course.

\section{Image Display Interface}

Data in {\it image(5NEMO)} format \index{image(5)} can be transferred in
{\it fits(5NEMO)} format and subsequently displayed and analyzed within
any astronomical image processing system.  They are generally much
better equipped to display and manipulate data of this kind of format. 
A number of standalone display programs can also understand FITS
format.  \index{fits(5)} An excellent example of this is 
{\it ds(1L)}\footnote{{\tt ds} is an X-windows specific 
program written by Bill Sebok, and is capable of auto-detecting
a number of different file formats}
\index{Sebok W.L.}, which is also available \index{ds} through NEMO.  If
compiled properly, {\tt ds} can even understand {\it image(5NEMO)} files
directly, without the need to go through {\it fits(5NEMO)} format 
using\index{ccdfits(1)}
{\tt ccdfits}.  Another powerful standalone display program which
can be used is {\it saoimage(1L)}, which runs as a server under
X-windows, and is used within the IRAF package.  \index{X-windows}
\index{SAOIMAGE} It also 
understands FITS format files directly.
% although the filename must end in {\tt ".fits"}.  %% NOT QUITE TRUE

\begin{verbatim}
    % nemosao
    % ds
    % saoimage -fits
\end{verbatim}

%Although the title of this section is short for {\tt IDI},
%\index{IDI}
%it seems unlikely that NEMO will use, or adapt to, an 
%image interface definition as {\it e.g.} IDI 
%\footnote{A set of image interface definitions, cooked up by
%the astronomical community (ESO/STI/etc.)} in the near future.
