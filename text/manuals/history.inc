%% This file is to be included by latex in nemo.tex
%
% chapter: Future, Present and Past
\myfile{history.inc}
\mylabel{a:history}

This Appendix contains a number of fairly random remarks on some
remaining difficulties in NEMO, and what we may do to resolve them, and
other things which may happen to NEMO over the next
years\footnote{The rumors of my death have been greatly
exaggerated...}.  A proper history is given of its past, and
can be found in itemized form in \ref{as:history}

\section{Some present problems}

\subsection{Graphics}

Yapp can only be linked in once, there's no dispatcher, and
different versions of the program exist for different graphics output
devices.

It would be nice to have more general dispatcher, which fully utilizes the 
{\tt yapp=} system keyword, may be used such that only one compiled
version of the program is needed.


{\it E.g.} on the SUN we could use:
\begin{description}
\item{yapp\_core: }   generic yapp using suntools (deprecated)
\item{yapp\_cg: }       extended suntools version with color support (deprecated)
\item{yapp\_ps: }       output to a postscript file, 
        for Postscript device. (color not yet implemented)
\item{yapp\_mongo: } interface which connects to mongo-87. (deprecated)
\item{yapp\_pgplot: } interface which connects to pgplot package. This
implementation also has
the great advantage of being able to handle a variety of
terminals and printers fairly transparent. This is the one we use mostly.
\item{yapp\_plplot: } interface which connects to plplot package. 
\end{description}

We use the environment variable YAPPLIB\index{YAPPLIB,
environment} in Makefile's.  This environment variable is also
set by your local NEMO guru in NEMORC.local \index{NEMORC, .local file}
file (normally
initialized by reading it by a users {\tt .cshrc} file) or you have to add
them specifically to your {\tt .cshrc} file.  This makes it easier to install
a new version of NEMO where a different graphics package will
been used. 

In some older versions of NEMO some of the Makefile's may have not
been modified to have these flexible setups, or the variable.

Sometimes, the need occurs for a specific YAPP interface, besides
default one.
Programs which specifically address a graphics device
have the base name
({\it e.g.} {\tt snapplot}, appended with an underscore and
the device ID, {\it e.g.} {\tt snapplot\_ps} or {\tt snapplot\_cg}. 
Best is to check the \$NEMOBIN directory.

The installer should carefully read the {\tt Makefile} and/or {\tt
README} file in the directory 
{\tt \$NEMO/src/kernel/yapp} for instructions regarding
specific installations of yapp interfaces. 

Programs would have to be recompiled manually, as in the following
example, because most Makefile's have hardwired graphics library
names in them:
\begin{verbatim}
    % make snapplot YAPPLIB="$YAPP_PS"
    % mv snapplot $NEMOBIN/snapplot_ps
\end{verbatim}

\subsection{System independent file structure}

\mylabel{sifs}

Currently the file structure is tied in with the operating system
routines {\it fread(3)} and {\it fwrite(3)}, and binary files cannot be
guarenteed\footnote{for IEEE and twos-complement data 
automatic byte-swap detection this problem has been solved} to be 
used across \index{data, format problem} machines with different data
types (size/floating point conventions etc.).  This problem is only
partially solved by using programs such as {\it tsf(1NEMO)} and {\it
rsf(1NEMO)}: it still requires physical data modification, transport and
again modification.  (UNICOS) cannot be read on a SUN, however datafiles
on e.g.  an Alliant, Multiflow and SUN are binary compatible because of
IEEE floating point \index{IEEE, floating} numbers and the proper size
twos complement integers.  Support exists now for byte-swapping, such
that files on Sun and Dec can be read and written either way.
However, machines like the CRAY supercomputer with its
deviating size and floating point format will have to convert their
data as is exemplified in Section \ref{s:exch-data}.

A possible solution is the way data is written to disk in 
a package\index{miriad, dataformat}
as {\tt miriad}: the layer just before the {\it fread/fwrite}
packs the data in some predefined standard format
(IEEE floating point and twos complement integers seem
an obvious choice at this time). This causes a small overhead
on some machines, and on other machines it is nothing more than
a copy operation or even passing of pointers.

\subsection{File size - float vs. double}

For a really large number of bodies (to take the example
of particle pushers) filesize becomes important for
many analysis programs which become I/O bound:

\begin{description}

\item{$\bullet$}
    It is not always necessary to keep information in double
precision. Images are also stored in doubles, in good
faith with old C, where all math is intrinsically done 
in double precision. This would save a factor of 2 in space.
There exist data i/o routines with force float/double
conversions ({\it get\_data\_coerced()}).

\item{$\bullet$}
    It is not always needed to keep all 6 phase space
coordinates, besides the data structure of snapshots has
phase coordinates rather well tied in, it's not easy yet
to separate positions and velocities, and only store
positions in a datafile. This would save another factor 
of 2.

\item{$\bullet$}
Images: totally unneeded to have them in double precision.
Should become {\tt 'float'}. Since we have used the {\tt real} type
\index{real}

\end{description}



\subsection{Some shortcuts and hints}

\begin{description}
\item{$\bullet$}
        You can set {\bf \$NEMO} to a 'universal' path as "/usr/nemo"
        and make a symbolic link of this file to the actual 
        physical location of nemo. In a file cluster system the
        actual location of NEMO might be on a common fileserver.
        Normally you have to become superuser to make the link
\begin{verbatim}        
                # cd /usr
                # ln -s /usr/guinness/nemo
\end{verbatim}
        This means that your {\tt .cshrc} can always keep the same
\begin{verbatim}        
              setenv NEMO /usr/nemo
\end{verbatim}
        and have the symbolic link take over the work for you.
 
\item{$\bullet$}
        The use of FLOAT\_OPTION is encouraged on the SUN3, because it allows
        a flexible change from a system with different floating
        point hardware.
        Perhaps the CFLAGS in most Makefiles should have the
        inclusion CFLAGS=\$(FLOAT\_OPTION) for compatibility
        reasons with other-than-SUN systems??? Beware that only the
        NEMO variable is exported to a Makefile (but see below)

\item{$\bullet$}
        When, while trying to compile, the cc compiler does not seem
        to find the NEMO include file, it is probably an indication
        that the {\tt cc} in {\tt \$NEMOBIN} has not the proper flags. In
        particular, some versions of the {\tt cc}-compiler do not support
        the -L flag (e.g. Ultrix, Sun UNIX 4.2 Release 3.1FCS)
        Make sure the '{\tt cc}' and '{\tt make}' are properly placed in the
        {\tt \$NEMOBIN}.

\item{$\bullet$}
On a SUN3 floating point intensive programs will run a lot faster when
the inline floating point libraries are used.  Instead linking with
the standard math library ({\tt -lm}), link it with {\tt
/usr/lib/f68881.il} or {\tt /usr/lib/ffpa.il}. The fortran TREECODE
only gains about 5\% in speed, but floating point intensive
programs can gain up to 30\% in speed.

\item{$\bullet$}
The best run-time performance from SUN4 compilers for
compute-bound applications is usually obtained from some combination of the
following compile-time options:
\begin{verbatim}     
     Fortran 1.3.1:
             -O4 -cg89 -libmil -dalign -fnonstd -Bstatic
     C 1.0:
             -O4 -cg89 -libmil -dalign -fnonstd -Bstatic -fsingle
\end{verbatim}         
These are discussed in the Numerical Computation Guide which accompanies C 1.0
and Fortran 1.3.1.  Also the default {\tt swap} and {\tt /tmp}
partitions supplied by
SunOS are often insufficient to fully optimize some large programs.  Use
{\it swapon(8)} in the first instance and {\tt -temp=...} compile option, 
described in {\it cc(1)} and {\it f77(1)}, in the second instance.
%           (dgh@validgh.com)


\end{description}

\section{Future}

A wishlist and what may be forthcoming in some future release of NEMO:

\begin{description}

\item{$\bullet$}
graphics: yapp\_server to work across machines. Not  necessarily
X11 server - but likely so. This finally means
a full implementation of the yapp=  system keyword and also would
make executables a lot smaller. Also a {\it plot(3)} interface.

\item{$\bullet$}
graphics: yapp\_x: An honest X11 graphics server.

\item{$\bullet$}
image display: more support for display facilities, e.g. {\it ds(1L)}.

\item{$\bullet$}
loadobj should understand constant expressions

\item{$\bullet$}
SPH : utilities. 
% Makes Lars happy.

\item{$\bullet$}
The everlasting expansion of the manual: more examples,
tutorials for course work, Figures, tables etc etc.

\item{$\bullet$}
loadobj for COFF (SUN386i, most SYS5 implementations, Convex).
(partially done, 3b1 version works) \index{COFF}

\item{$\bullet$}
???shared libraries for the latest SUN OS 4.1 system???
Is a rather laborious thing, and difficult to maintain
in a environment where the library is frequently upgraded.
\index{shared libraries}

\item{$\bullet$}
dynamic object loader for N-body diagnostics
\index{diagnostics, N-body}

\item{$\bullet$}
options.h through stdinc.h or nemo.h?

\item{$\bullet$}
handle multi-snapshot files more efficiently

\item{$\bullet$}
loadobj for MF, UNICOS, Alliant and Convex do not work!!
The only reliable implementation we have is BSD (SUN OS) and Ultrix
and a SYSV (3b1). SPARC also seems to be stable, and MIPS COFF is
not quite done yet. \index{MIPS COFF}

\item{$\bullet$}
loadobj in yacc?, nemoinp in yacc?  (cf. sm)

\item{$\bullet$}
Utility for chaining programs in (nsh?) shell scripts, automatic passing
of in/out files from on to the next. \index{chaining programs}
See {\tt pipe} shell script for example.

\item{$\bullet$}
Install the official NBODY1 and NBODY2 programs from Aarseth
\index{Aarseth S.} with NEMO's user interface and file structure,
as has been done for NBODY0

\item{$\bullet$}
Adaptation of all C programs to the ANSI standard. Usage of the 
GNU gcc compiler\index{gcc}\index{gnu, cc}
recommended for portability. This process is 
now underway (march 90 - PJT).
It is also likely that the GNU make program will be used on the long run,
and making small updates to the package will be a lot easier.
(more portable  {\tt make}?)\index{gnu, make}

\end{description}

\newpage

\section{New Features}

This section is not updated frequently, for more timely information it is probably
better to consult
{\tt http://www.astro.umd.edu/nemo/whatsnew.html}.

\subsection{Release 3.3}
To be released around the 2nd NbodySchool (Amsterdam, July 23-30 2005). Improved
support to help installing various ancillary packages needed for the summer school
(dubbed ``manybody'' in {\tt \$NEMO/usr/manybody}). The I/O library was
updated with the ability to handle blocked I/O.

\subsection{Release 3.2}
Was released April 11, 2003, just after the Strasbourg N-body school. The multi-CPU
directory tree often made the bootstrap installation with libraries such as
pgplot, gsl, cfitsio etc. harder. Also the installation on MacOSX 10.3 was now
streamlined, though still not perfect "out of the box" like it still does on 
most solaris and linux boxes.

\begin{description}

\item{$\bullet$}
back to a top level {\tt make bins}, which creates "all" binaries. It currently creates
about twice as many binaries as the old more reliable/robust 
{\tt src/scripts/testsuite -b} method.

\end{description}

\subsection{Release 3.0}
Was released April 1, 2001. Certainly not a joke. Source code has been released within
CVS for subsquent development with the partiview and starlab modules, also under
CVS.
Installation has now completely been done using configure, an
autoconf product. 3.0.0 was never released, 3.0.1 was the first official release that
worked under linux. Solaris and SGI.

\begin{description}

\item{$\bullet$}
configure support

\item{$\bullet$}
loadobj support now using .so files, not .o files.


\end{description}


\subsection{Release MD-2.5}
Was released December 1999, with initial support for configure. This series had 5 subreleased
though 2.5.5, but 2.5.6 was never released in favor of 3.0.0. No major changes through this
release.

\subsection{Release MD-2.4}
Was released April 1, 1997. 

\subsection{Release MD-2.3}

\begin{description}

\item{$\bullet$}
We are now using an official versioning scheme, with major,
minor and patchlevel. 

\item{$\bullet$}
The dreaded dynamic object loader\index{loadobj} has been
made to work on the DEC Alpha (OSF1 V3.2) and SGI (IRIS 5.2)
using a new system\index{ldso} utility {\tt ldso}. This has increased
the portability, but not taken away the fundamental difficulty
in installing NEMO on a new operating system.

\end{description}

\subsection{Release MD-2.2}
\begin{description}

\item{$\bullet$}
Full support for Solaris 2.x. (tested on 2.3) 
The default graphics {\tt yapp}
\index{yapp, default install} device should now be generic X
and Postscript. loadobj\index{loadobj, installation} method 
is now supported by the operating system.

\item{$\bullet$}
Compiling with {\bf -DSINGLEPREC} actually caused lots of programs,
include {\tt hackcode1}, not to run properly. Most of the bugs 
associated with this have been removed, but for example LINUX
is very sensitive to such errors and still a lot of programs will
crash in this mode.

\end{description}

\subsection{Release MD-2.1}
\begin{description}

\item{$\bullet$}
{\it Literate Programming:} trying CTEX\index{ctex}
embedded comments in e.g. the potential\index{potential, ctex}
descriptors, and
{\tt anisot.c}\index{anisot, ctex}

\item{$\bullet$}
Properly documented and advertised use of using binary
structured files in pipes.

\item{$\bullet$}
Some support for Starlab in the form of NEMO++. Translation
programs to convert from a {\tt dyn} to a {\tt snapshot} have
been written. Apart from properly ANSI-coding the NEMO kernel,
the routine {\tt nemomain.c} needs to be present in C++ format
too: {\tt nemomain.C}.

\item{$\bullet$}
Potentials are now supposed to return a pattern speed (even if
it was not changed) into the first argument of the potpars array.
This to deal with rotating potentials.

\end{description}

\subsection{Release MD-2.0}

With this major release upgrade the directory structure has been modified
away from 'user' oriented to 'topic' oriented.  The 'user' oriented stuff
is now under {\tt \$NEMO/usr}, whereas the more stable 'topic' oriented
under  {\tt \$NEMO/src}. Some files live in both, in which case the
{\tt \$NEMO/src} version should take preference.

\begin{description}
\item{$\bullet$} {\it dlopen()} version of loadobj, as well as 
{\it ldl} (``gnu''). Still some problems, but they are Sun bugs,
not ours. Certain complicated expressions fail. \index{loadobj}
{\it dlopen} may also work on NeXT.

\item{$\bullet$} {\it yapp} and {\it loadobj} separated out and cleaned up

\item{$\bullet$} small subset of numerical recipes maintained

\item{$\bullet$} stories in snapshots \index{story, in snapshot}

\item{$\bullet$} random access additions to filestruct, and support
\index{random access}
for little and big-endian machines using auto-byte-swapping. Hence
no Macro's used.  Slow?

\item{$\bullet$} {\it bodytrans} using {\it fie}.
\index{bodytrans(3)}

\end{description}

\subsection{Release MD-1.4}
Yet another non-release - summer 1990.
\begin{description}
\item{$\bullet$}
Miriad shell implemented as {\tt nemo} through an alias. \index{nemo,shell}
The old {\tt nemo} program renamed to
{\tt nemoshow}. \index{nemoshow} The {\tt mirtool} can also be compiled
to {\tt nemotool}.

\item{$\bullet$}
Sault's FITS I/O routine replace Werong. No real need for a Fortran-C
interface anymore. New FITS routines in image.

\item{$\bullet$}
Experimental Micro-NEMO in {\tt src/nemo/micro}.

\item{$\bullet$}
Manual now reaches about 100 pages.

\end{description}

\subsection{Release MD-1.3}
Newly released option until February 1990.

\begin{description}

\item{$\bullet$}
potential(5NEMO) has extra time parameter, relevant programs have been
updated.

\item{$\bullet$}
hackcode3 \index{hackcode3} is an
experimental version which allows an extra external potential through
the potential(5) format. It also has the option of keeping the first
{\tt nrigid} particles rigid, and more silly options are expected.
It should be noted that hackcode2, an experimental toy, has never
been released.

\item{$\bullet$}
The user interface {\tt getparam()} now supports reading the value(s)
of a keyword from a file using the {\tt key=@file} construct.

\end{description}
\subsection{Release MD-1.2}

Newly released option until November 1989
\begin{description}

\item{$\bullet$}
Introduction of the HOSTTYPE environment variable in the NEMORC file.
This meant that
new environment variables such as NEMOBIN, NEMOLIB and NEMOOBJ
are derived from NEMO and HOSTTYPE. See the NEMORC file.
It also meant that basically all Makefile's had to be updated.

\item{$\bullet$}
Libraries are a bit screwed up now, and best is to include
both libT.a and libJ.a in the minimum list of libraries. Makefiles
are being updated for this.
One can also use libNEMO.a, and use utilities such as
mklib, addlib and mkbin.

\item{$\bullet$}
bodytrans enhanced, has a proper database of expressions, which
is dynamically updated when bodytrans(3) is used. The file
{\tt \$NEMOOBJ/bodytrans/BTNAMES} contains a list of extra 
expressions currently understood. Make sure those directories are
write permissible by the world.

\item{$\bullet$}
Loadobj is now also functional for Sparc (sun4) and COFF (System V UNIX).
The SPARC version requires loading with {\bf -Dsun4 -Bstatic} since
often required symbols, such as integer multiplication, are in
shared libraries. This unfortunately makes those binaries larger
than they could be.

\end{description}

\subsection{Release MD-1.1}

A description of the major differences between the existing IAS version 
and the newly released MD version 1.1 (summer 1989):

\begin{description}

\item{$\bullet$}
 Improved user interface: suntools menu interface at help=8.
 The sophistication of the user interface is determined at
compile time through a 
number of
 compile switches in {\tt getparam.c}

\item{$\bullet$}
Standardized usage of some standard interactive facilities using
    {\tt setparam} instead of and with help of {\tt getparam}.

\item{$\bullet$}
Installation through Makefiles is more flexible. The installation
and various administrative
utilities are more streamlined to make porting to
non-SUN's easier.

\item{$\bullet$}
Aarseth's straight N-body code - as published in BT87 is in
nbody0 and an Ahmad-Cohen version in
nbody2 (this last one has not been tested out well enough)
These programs have a NEMO interface, which also handles some
Fortran-C interface questions. \index{Aarseth S.}

\item{$\bullet$}
filestruct\_old is now obscure: \index{filestruct, old}
all of Piet's clib programs have been
converted (snaplist, snapenter, snapdist, mkplummer) \index{Hut P.}

\item{$\bullet$}
various new yapp's, the \$NEMO/src/nemo/yapp directory has been
cleaned up and documentation has been updated.

\item{$\bullet$}
bodytrans saving: the program bodytrans can be tested and also save
files ({\bf -DSAVE\_OBJ} compile switch) - same for library routines.

\item{$\bullet$}
snapplot and movies work better together; snapplotedit.

\item{$\bullet$}
potcode

\item{$\bullet$}
ds, image display on image(5) or fits(5) files.

\item{$\bullet$}
fortran interface


\item{$\bullet$}
GRAVSIM added to nemo/src tree \index{Bellon M.}

\item{$\bullet$}
lars/treecode has been added, but needs a good  \index{Hernquist L.}
recursive fortran compiler. NEMO interfaces built in. Not been able
to test well - Does not work on Suns.. \index{Hernquist L.}

\item{$\bullet$}
The NEMO startup file is now called {\tt NEMORC}. 


\item{$\bullet$}
This expanded manual.

\end{description}

\newpage

\section{HISTORY of NEMO}
\mylabel{as:history}

\begin{description}

\item[Oct 86 - Jun 87]
Initial development on a network of SUN3 workstations
at the Institute for Advanced Study, Princeton, NJ by
Barnes, Hut and Teuben.
\index{IAS}

\item[Jun 87]
export version, for easy installing on BSD4.2, we call it
alpha version 1.0a.

\item[Jul 87]
test phase for installing on an VAX 8300 running
Ultrix 1.2: 
mods: Makefile, new cc and make (cc has no -L flag on ultrix) tsf.c
casting of pointers and advancing pointers fixed

\item[-Mar 88]
        various installations at MIT, Drexel U. in Philadelphia,
        U.of Illinois at Urbana (Convex, Alliant, Gould) tested and
        done. New programs and updated
        programs keep coming in at a steady rate.

\item[Mar 88]
        Many UIUC changes: getparam() has a few new system keywords,
        (debug, yapp, host).
        History mechanism in data file
        I/O means that now get\_hist() and put\_hist() must be
        called, optionally add\_hist(). 
        dprintf() added to getparam.c for user debugging aid, can now
        be called in user programs. 
        Yapp\_mongo used. Environment variables
        YAPPLINT and YAPPLINP are now encouraged in user Makefiles.
        Documentation into one big TeX file for users as well as programmers.

\item[Jun 88] IAS and UIUC versions have been merged again. Improved
        filestructure (Josh), user interface, data history mechanism,
        yapp interface (Peter). Documentation end manuals significantly 
        improved. Still a number of items in the 'problems' area not resolved.

\item[Nov 88] IAS and UIUC versions have been merged again.

\item[Spring 89] Manual updated - working on class/course problems - lots
        of working examples added to manual. MD version is now
        slowly diverging from IAS version.

\item[Summer 89] Groningen version installed  on a combined
    SUN/Alliant network - shows difficulty of maintaining 
    a shared disk environment with different versions of
the binaries (binX, libX, objX, datX) - 

\item[August 1989] Version prepared for official beta release 1.1.
(Maryland)

\item[December 1989] Experimental multi-CPU release (1.2)
(Maryland)

\item[February 1990] Minor upgrades, mainly SUN4 and multi-CPU
bugs. (1.3). 
Toronto's version has been dubbed ZENO now. \index{ZENO}

\item[May 1990] Slightly expanded manual for the Pittsburgh Workshop
and some minor upgrades every here and there.

\item[Summer/Fall 1990] An attempt to merge Starlab and Nemo - added
some extra
functionality to filestruct and merged story concept.
Total directory structure overhaul:
all code related to a particular topic is in its own directory in
{\tt \$NEMO/src}. The tree starting at {\tt \$NEMO/usr} will now be used 
for user contributed software. In particular, {\tt \$NEMO/src/kernel}
contains a small core of NEMO which can be used without any of
the parallel branches.

\item[Fall 1991] Slowly progressing the {\tt \$NEMO/src} tree. A few new
sites for export maintained.

\item[Summer 1993] Some support for C++ and Starlab.

\item[January 1994] Solaris 2.x support, started WWW.

\item[Winter 1995] Linux, Dec Alpha and SGI support added, since they allow
dynamic object loading. Manual in html (latex2html)

\item[Spring 2001] Installation converted to configure/autoconf, also 
using CVS for source code control now.

\item[April 2004] Slight directory change for the directories created
during the installation, this ends the era of the multi-CPU tree.
Also removed a last environment variable that was used in Makefiles,
now they are all inherited from the ones created by configure.
Released between the Strasbourg Nbody-school and the College Park FAM04
tutorial/workshop weeks.

\item[Summer 2005] Second summerschool, at MODEST-5c, In Amsterdam. Manual
expanded with more examples.

\end{description}
