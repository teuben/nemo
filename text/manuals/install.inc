%    -*- LaTeX -*-     for inclusion in nemo.tex
%
% chapter Installing NEMO, updates and exporting
%
\myfile{install.inc}
\mylabel{a:install}

This Appendix describes installation and maintenance of NEMO. 
NEMO currently consists of almost 60MB of source,
documentation, datafiles etc.
obtained either as a compressed tar file (around 6MB) or, preferably 
actually, via git\index{git}. As of December 2017
we no longer maintain the CVS \index{CVS} release.
In NEMO V4 the installation has been simplified using autoconf, and
users are encouraged to use git to install and update  their version
of NEMO.

The first Section \ref{as:suninst} walks you through an installation 
using the new configure script.
For the most it is likely to be close to the
procedure on any standard UNIX system, although a number of known
deviations are described in Section \ref{as:otherinst}. Exporting a
working NEMO system is described in Section \ref{as:export}. Incremental
updates (import as well as export) are described in Section
\ref{as:update}. Finally, various maintenance issues are discussed in
Section \ref{as:maint}.

The installation examples shown below assume you are using something
like the {\tt csh} shell; you would likely have to replace 
commands such as {\tt setenv}, {\tt source}, {\tt rehash}, {\tt set} etc.  
by the appropriate ones for another shell.  We mereley leave this as
an excersize for the reader (:-).

\section{git - NEMO V4.0}

Most of this is now described in more detail in the file 
{\tt README.install}, but
briefly here are the steps of a typical current NEMO install, using PGPLOT,
in csh (there still is no proper sh support):

\begin{verbatim}
  git clone https://github.com/teuben/nemo
  cd nemo
  (mkdir local; cd local; git clone https://github.com/teuben/pgplot)

  ./configure --with-yapp=pgplot --with-pgplot-prefix=`pwd`/lib
  source nemo_start.csh
  make postconfig
  rehash

  make pgplot
  make vogl
  make libs
  make bins
 
\end{verbatim}

This process takes a few (2-3) minutes 
and produces about 200 programs and occupies about 125MB of diskspace

\section{Bootstrap}


\section{Linux workstations - NEMO V3.0 and V3.1}
\mylabel{as:suninst}

We are going to install NEMO in the directory designated
by the environment variable {\tt\$NEMO} on a Linux
workstation. \index{NEMO, environment} The {\bf NEMO} 
environment variable will be needed (see also 
Appendix~\ref{a:setup}),
and will generally be set 
in your {\tt .cshrc} file.  We will do this a little
later.  Even if your are not on a Sun workstation, please read this
section carefully and the next section on installation procedures on
non-Sun UNIX systems.  \index{.cshrc}

If you haven't done so yet, create NEMO's root directory now
(or make at least sure that the directory is empty):

\begin{verbatim}
    % mkdir $NEMO                (or whatever $NEMO will be)
    % cd $NEMO
\end{verbatim}

It is perhaps convenient to install a special ``{\tt nemo}'' user.
In case your site has active users that are donating software to NEMO
(the ultimate goal\index{goal, NEMO} of this project)
the sub-directories need to have write permission for them. The use
of UNIX 'groups' may come in handy here. To begin with, {\tt \$NEMO/usr}
is intended for this purpose.

We will now take you through the installation step by 
step\footnote{Remember there is also the {\tt bootstrap} method, which
automates large parts of the following steps}.

\subsection{tar}

A tar ``image'' of the source code
will be either on tape, or is already 
available\footnote{See Section \ref{as:export} how to create
such tar files from an existing NEMO installation}
on disk. Extract it from tape \index{import, tar image}
(the\index{tar, import}
name of the tapedrive in this example is {\tt /dev/rmt0}):

\begin{verbatim}
    % tar xf /dev/rmt0
\end{verbatim}

or extract it from a disk tarfile (called {\tt nemo.tar} here):

\begin{verbatim}
    % tar xf nemo.tar
or:
    % zcat nemo.tar.Z | tar xf -
or:
    % gzip -dc nemo.tar.gz | tar xf -
or:
    ftp> get nemo.tar.gz "|gzip -dc | tar xf -"
\end{verbatim}

depending on the nature and location of the (compressed) tar  file. It
may have been compressed (\verb+.Z+)  with the standard UNIX {\it
compress(1)}\index{gzip} program\index{compress}, or the GNU utility
{\tt gzip}. Instead of using {\it uncompress(1)},  the\index{uncompress}
utility {\it zcat(1)}\index{zcat}  can be piped to {\it tar(1)} to
extract the tar file much more efficiently.   You may even transfer the
file accross the (inter)net  via the {\tt ftp}\index{ftp}\footnote{not
all versions of {\tt ftp} allow this feature though} utility and piping
it's output straight into {\tt tar} to avoid having to store the tar
file on your local disk.

BE SURE to be in the {\tt \$NEMO} directory while doing this! On some
UNIX systems it may happen that {\tt tar} complains about not being able
to create directories or such. You can probably ignore it, since  SYSV
{\tt tar} wants to preserve ownership (stored by a user id, which may be
unknown or even invalid on your target machine). The {\tt tar} program
may have a flag for this.

If you haven't set the  {\bf NEMO} environment\index{NEMO, environment}
variable (cf. Appendix~\ref{a:setup}) now is really a good time to do
it. 

First add the following to your {\tt .cshrc} file, somewhere before the
{\bf PATH} environment\index{PATH, environment} is defined:

\begin{verbatim}
   setenv NEMO /usr/nemo                        (or whatever)
   source $NEMO/NEMORC
\end{verbatim}

%% The {\bf HOSTTYPE} environment variable \index{HOSTTYPE, environment}
%% has been added here too; note the backquotes if this is the way you want
%% to set it.  It is set to a value such as {\tt sun3}, {\tt sun4} or {\tt
%% convex}, to be able to handle multi-machine environments on the same disk. 
%% If you do plan to support only one machine, the {\bf HOSTTYPE} could be
%% set to a blank string ({\tt ""}).  

After having read the {\tt NEMORC} startup file, a number of environment
variables and alias will have been added to your current environment.
Include the {\tt \$NEMOBIN} directory
in your search path before any major system
areas, {\it e.g.}:

\begin{verbatim}
   set path = ( . $NEMOBIN $path)
\end{verbatim}

and make it current:

\begin{verbatim}
   % rehash
\end{verbatim}

Now that all source and documentation is in place, the installation can
start.  As said before, the top level {\tt \$NEMO/Makefile} can
steer the whole\index{makefile, install} installation\index{install,
makefile} process heirarchically, calling {\tt Makefile}'s in lower level
directories and so on and so on.\footnote{At least, this is the goal of
this game}

\subsection{make dirs}

First we need to complete the creation of the NEMO directory
structure\footnote{A full listing was presented in Appendix \ref{a:dirs}}
as needed for this particular {\tt \$NEMOHOST}. 
Also a few files within
the system need to have write permission by all users:

\begin{verbatim}
    % make dirs
\end{verbatim}

You will likely see comments that certain directories and files
exist, which is fine. This install command can be rerun without
any danger when updates come in. If you get some kind of error
message like:

\begin{verbatim}
    Badly placed ()'s.
    *** Error code 1 (ignored)
\end{verbatim}

it is most likely that you are within a strict SYSV-like environment, that 
directly imports the SHELL variable into {\tt make}, 
and since\index{SHELL, environment}
our makefile's are written with the assumption that {\tt /bin/sh}
is the make-shell, one has to use something like:

\begin{verbatim}
    % make dirs SHELL=/bin/sh
\end{verbatim}

By copying the appropriate {\tt make.*} script from 
{\tt \$NEMO/src/scripts/} into {\tt \$NEMOBIN}, this problem
is generally solved.


Also a few files are created that are needed when running NEMO.  These
are summarized in Table \ref{t:instfiles}. If you expect to run
the import/export scripts and be in regular contact with the central
source archive, you need to give them your 
{\bf NEMOSITE}\index{NEMOSITE, environment} name alias,
and edit the {\tt NEMORC.local} file\index{NEMORC.local, file} appropriately.
The {\bf YAPPLIB} probably needs to be set to one of the available
{\bf YAPP\_*} variables from the {\tt NEMORC} file, since the default
is {\bf YAPP\_NULL} (no graphics output!).
The sensible thing to do, is going through the installation of the
library (as discussed in Section \ref{ss:install}),
after which some tests in the directory
{\tt \$NEMO/src/kernel/yapp} will tell you which graphics device
is easiest to install on your machine.

\begin{table}[ht]
\caption[Installation created files]{Installation created files}
\mylabel{t:instfiles}
\centering
\medskip
\begin{tabular}{|l|l|l|}\hline\hline
%
File    & Purpose   & Notes \\ \hline
%
etc/motd      & New announcements & You edit to suit your site \\
adm/TIMESTAMP & installation      & do not hand {\tt touch} this! \\
adm/Usage     & usage log         & needs {\tt chmod a+w} permission \\
NEMORC.local  & import/export     & edit at least NEMOSITE \\
nemo\_start   & easier startup    & not created with {\tt make} \\ 
VERSION       & version id        & must be {\tt major.minor.patch} numbers \\
&& \\ \hline
\end{tabular}
\end{table}

\subsection{make scripts}

The next step is to copy a number of 
important\footnote{only the {\tt cc} script is installed now}
scripts into {\tt \$NEMOBIN}:

\begin{verbatim}
    % make scripts
\end{verbatim}

Currently the NEMO environment needs a search path where the system
C-compiler is preceded by our\index{cc, script} own {\tt cc} script. 
This script then calls the real UNIX C-compiler (often {\tt /bin/cc}). 
On some UNIX machines, a {\tt make} script needs\index{make, script} to
be installed too to precede the system {\it make(1)} utility.  This will
be discussed later. 

If you don't like the choice the installation procedure has given you,
pick another one, and copy it yourself instead, for example, to select
the gnu\index{gnu, cc compiler} compiler, issue:

\begin{verbatim}
    % cd $NEMO/src/scripts
    % cp cc.gnu $NEMOBIN/cc
    % rehash
\end{verbatim}

\subsection{More on scripts}
\mylabel{ass:scripts}

In the previous section only a few very essential scripts were copied
to {\tt \$NEMOBIN}, but the directory {\tt \$NEMO/src/scripts} contains
a number of other useful scripts for users,
as summarized in Table~\ref{t:uscripts}.

\begin{table}[ht]
\caption[User scripts]{Useful user scripts for NEMOBIN}
\mylabel{t:uscripts}
\centering
\medskip
\begin{tabular}{|l|l|l|}\hline\hline
%
Script    & Purpose   & Notes \\ \hline
%
bake\index{bake}	  & make replacement	      & Template in {\tt \$NEMOLIB} \\
mknemo\index{mknemo}	  & Find and install programs & \\
mkpdoc\index{mkpdoc}      & check if doc file needs update & Also called by mknemo \\
manlaser\index{manlaser}  & send manual pages to (laser)printer & Ignores small .so files\\
rmsf\index{rmsf}	  & remove NEMO binary data files & \\ 
&&\\ \hline
\end{tabular}
\end{table}

There are some more scripts, but they are more useful for 
NEMO maintenance, as summarized in Table~\ref{t:mscripts}, and
should perhaps not be copied to {\tt \$NEMOBIN} for direct access
by users.

\begin{table}[ht]
\caption[Maintenance scripts]{Useful maintenance scripts}
\mylabel{t:mscripts}
\centering
\medskip
\begin{tabular}{|l|l|l|}\hline\hline
%
Script    & Purpose   & Notes \\ \hline
%
nmlist\index{nmlist}    & listing of library          & \\
ranlib\index{ranlib}    & fake if none present        & needed on SYS5 \\
tardot\index{tardot}	&   gather all hidden . files & see also submit \index{submit} \\
&&\\ \hline
\end{tabular}
\end{table}

\subsection{make install}
\mylabel{ss:install}

A full installation starts with installing the essential libraries,
after which application programs can be compiled at will.

On a Sun4 workstation the default installation is invoked by:
\begin{verbatim}
    % make -i sun4 >& make.log &
    % tail -f make.log
     .....
      <Control-C>                       #  to abort tail
\end{verbatim}
We have added an (stdout + stderr) redirection of the screen output
to a file here and also put the job in
the background for convenience.
Also, the extra {\tt -i} flag will force the installation through
possible errors.
The log file {\tt make.log} is to check for errors afterwards.

Particularly at this stage, before you compile source code, you may wish
to review if certain system dependant features compiled correctly, and
patch the system.  No good facilities are available from the top level,
but {\tt Makefile}'s in the appropriate directories usually contain
information and are flexible enough to patch the system with a minimum
amount of effort.  A number of standard patches are described in the
next Section. 

\subsection{mknemo}

The previous process only installs the necessary libraries. At this stage
each program needs to be loaded manually.  The script {\tt mknemo}
should be installed for this (see section \ref{ass:scripts}).  To
install the programs {\tt hackcode1, mkplummer} and {\tt snaprstat} you
could issue:\footnote{you may need a {\tt rehash} command here}

\begin{verbatim}
    % mknemo hackcode1 mkplummer snaprstat
\end{verbatim}

The {\tt mknemo} script is a little bit smart about certain different
structures of compiling programs, {\it e.g.} if there is no single
sourcecode found (as with {\tt hackcode1}), it does accept a directory
with that name, after which the local {\tt Makefile} is allowed to take
over.  Clones, such as {\tt hackcode1\_qp} cannot be installed this way
though. 

Don't forget to run the {\tt rehash} command if a command was new.
Everything should be set up to use NEMO! 

\subsection{Documentation}

The documentation of NEMO is manyold: first of all there is 
a {\it Users and Programmers Guide}
(what you are reading now), which is currently in LaTeX 
format.\index{documentation, latex file}\index{LaTeX}\index{makeindex}
Before you format this huge document, be sure to check if
the {\tt nemo.dvi} or even the (compressed) {\tt nemo.ps} postscript
file are not present in the manuals directory:

\begin{verbatim}
    % cd $NEMO/text/manual
    % zcat nemo.ps.Z | lpr
\end{verbatim}

or compile it:

\begin{verbatim}
    % cd $NEMO/text/manual
    % make nemo
\end{verbatim}

and print it:

\begin{verbatim}
    % dvipr nemo    # or whatever your local version is called
\end{verbatim}

The program {\tt makeindex} is needed by LaTeX, it is a public domain 
utility, see the local Makefile for more details. The actual commands
may be slightly different on your host.
If {\tt makeindex} does not compile or work for you,
a dummy zero length file {\tt nemo.ind}
must be created, before LaTeX can be run. You will have to live
without an index in this case, or ask the distributors for
the {\tt nemo.dvi} or proper {\tt nemo.idx} file. 
Be also aware that the manual ({\tt nemo.tex}) uses include 
files ({\tt .inc}) using the\index{inc, files}
\TeX\ ``\verb+\input+'' command, and may also include direct postscript
files using the special {\tt $\backslash$ PSinsert} command. 
There are a few logical ``flags'' using the \TeX\  {\tt $\backslash$ newif}
command, {\it e.g. } in the top of the {\tt nemo.tex} document they are
declared as:

\begin{verbatim}
    \newif\ifnemo         % declaration of the logical "variables"
    \newif\ifdebug
    \newif\ifindex
    ...
    \nemotrue             %  setting the variables to true or false
    \debugtrue
    \indextrue
    ...
    \ifindex              % example use of such a conditional 
       \printindex
    \fi
\end{verbatim}

In particular, the {\tt $\backslash$ ifdebug} variable is used to comment
the margins with all the items that go into the index. This makes 
proofreading the manual very efficient.

The second form of documentation are online manual pages. For all
programs there is (or at least should be) a 
manual page (see {\it man(1)}). These can be individually
printed using a shell-script {\tt manlaser} or:
\index{manlaser}\index{troff}\index{documentation, manual pages}
\begin{verbatim}
    % manlaser program.1
or
    % troff -man -t program.1 | lpr -t 
\end{verbatim}
The command script {\tt \$NEMO/usr/nemo/maint/PrintMan} can\index{PrintMan}
be used
to print out {\bf all} manual pages 1,3,5 and 8, but also check the 
{\tt catman}  \index{catman} target in the toplevel Makefile.
This creates\index{fmt, files}
{\tt fmt} files, that can be sent directly to the printer. Be careful
about redundant information because of the {\tt .SO} troff command. The
{\tt fmt} are normally logical links \index{linked, manual files}
in this case, the script should
take care of this, and only print out proper files and not links
or zero length files.

Third, there are the {\it doc} files. They are needed by various
front-end shell (as described in Appendix \ref{a:iface}).
The script {\tt mkpdoc} will\index{doc, files}
check if an update of this {\it doc} file is needed; this script
is called by {\tt mknemo}, the quick\index{mknemo}
and dirty procedure to install programs. In principle 
the {\it doc} files don't
have to be saved, since they can be generated from the source code.

Fourth, there are the {\it ctex} files. In a number of places the
source code contains mathematical details in TeX format. A small
utility, {\tt ctex}\index{ctex, tex extractor}, extracts the
these from the source code, and standard {\tt tex} utilities can
be used to format them. For example, most of the {\it potential(5NEMO)}
descriptors have an associated CTEX section, from which Appendix
{\ref a:potential} can be generated automatically. Other examples
are {\tt anisot.c} and {\tt mkop.c}.

\subsection{SUN only: mirtool, ds, movie's}

A number of window based utilities are available for Sun 
workstations. To install them, you have to do it after the
general NEMO installation, as some of them need the NEMO 
library. 

See {\tt \$NEMO/Makefile} for details how to install them.


\section{Tailoring}
\mylabel{as:otherinst}

NEMO has also fairly successfully been ported to various other UNIX
machines, such as\index{Convex}\index{Gould}\index{VAX}\index{Sun} a
CONVEX supermini computer, several GOULD machines, a VAX workstation
running Ultrix 3.0, a Multiflow Trace-14  and Cray Unicos.
They all required a few simple modifications, about which this
section reports. 

The first few (small) steps in the top-level Makefile
require some attention, in specific the cc-compiler and
the make-utility have been assumed to have too many options
which are generally not true on all versions
of UNIX. Note that these 'inflexibilities' may disappear with
time when BSD and SYSV UNIX merge
(expected release SUN OS 5.0 1992?). 
% Gould has recently added the
% '-L' and '-I' flags to their {\it cc(1)} compiler, as do more vendors
% these days.

\subsection{Replacement scripts in NEMOBIN}

Before NEMO V3 programs like cc, make and ranlib were sometimes 
replaced with a script in \$NEMOBIN.

\subsection{YAPP graphics device driver}

% \underline{{\tt YAPP graphics package:}} src/kernel/yapp/*

\index{yapp, installation} 

In compiling the major libraries and utilities in NEMO the proper
graphics package should be installed too.  The {\tt
\$NEMO/src/kernel/yapp} directory contains a number of yapp
implementations and your need to decide which yapp interface(s) to use. 

To test a particular graphics device you can use:

\begin{verbatim}
    % cd $NEMO/src/kernel/yapp
    % make yapptest_mongo
\end{verbatim}

where the trailing {\tt \_mongo} part tests, in this case, the mongo
device driver. 

This usually means that the {\bf YAPPLIB} environment variable in
\index{YAPPLIB, environment} the startup script
{\tt \$NEMO/NEMORC.local} has
to\index{NEMORC.local} be adjusted accordingly. 


\subsection{Math libraries} 

There are one or two programs which use {\it Numerical Recipes in
C} \index{Numerical Recipes} routines, and the library {\tt libnumrec.a}
should be present in {\tt \$NEMOLIB}, or be symbolically/logically
linked to the real one if references remain unresolved during linking.
Note we might want to use double instead of
float. Some other routines use IMSL, and we also supply a few 
emulated IMSL\index{IMSL, emulated}
routines in {\tt \$NEMO/src/kernel/misc/imsl.c}


\subsection{LOADOBJ dynamic object loader}

\index{dynamic loader, loadobj} \index{loadobj, installation}

This package is still a major problem on most non-BSD UNIX
implementations, and unfortunately many useful programs depend on it's
functionality.  The BSD version (Sun3, Sun4, Ultrix) as well as SystemV
COFF version are fairly stable, but versions for complicated cpu's or
non-standard object code loaders, such as Alliant, \index{Multiflow}
Convex, Unicos and Multiflow have \index{Convex} given \index{Alliant}
problems.  \index{Unicos} See source code documentation in {\tt
loadobj.c} must contain an entry for the proper {\tt loadobjXXX.c} file. 
Good news is that the new SYSV4 revision, as well as Sun OS 4.1 include
the {\it dlopen(3x)} library utilities for dynamic linking. 
\index{dlopen} \index{dynamic loader, dlopen}

\subsection{Tailoring the NEMO kernel}

It can be handy to know what a minimum subset of library routines
will be in order to run major portions of NEMO.  Some of the comments
below have never been tried out, so beware! Also some comments apply how to
tailor certain routines for less functionality, if
your local hostmachine is not as sophisticated. It is our intention
to have them defined in a file {\tt options.h}, although this 
may complicate multi-cpu 
shared source code. In such cases they may need to be defined
via the {\it cc} script {\bf -D} command line option, just as we currently
use the {\bf -I} switch.

The source code of most of the barebone NEMO library routines can be
found in {\tt \$NEMO/src/kernel/} and below.

{\bf loadobj.c:} This has been discussed before, a proper symbol
must be defined, in order for {\tt loadobj.c} to include the correct
loadobj source code type. Note that the trigger \index{trigger, cpp}
\index{cpp, trigger} may be used for you if no symbol has been
supplied. usually {\it cpp} supplies a trigger, such as {\it sun, sparc,
unicos, mips} etc.

{\bf bodytrans.c:} This routine is normally compiled with the {\bf
-DSAVE\_OBJ} compilation switch \index{bodytrans(1)} in which case
bodytrans functions, which are not yet found, will be saved to the
standard repository directory {\tt \$NEMOOBJ/bodytrans}.  
this option is now standard defined in {\tt options.h}.

{\bf getparam.c:} This file contains many {\tt \#define}'s
which allow you to add a number of extensions to the
user interface which have been discussed in Appendix \ref{a:iface}: 
data history
mechanism, interactive keyword and menu prompting, parsing of
expressions of {\it getXparam()} variables, interrupt to the review
section, remote machine execution, etc.  For more documentation see the
source code. none of these run through {\tt options.h} though.

{\bf image.c:} This routine handles all image i/o
({\it e.g.} read\_image, write\_image, create\_image), and can be 
compiled in {\bf -DCDEF} or
{\bf -DFORDEF} mode, depending how your matrices in c should be
stored in memory. The fits interface may not work
properly when the wrong def is used???

\subsection{isolation}

If certain portions of NEMO's subroutine library appeal to you,
yet you don't want to drag the whole NEMO library along (this may
not always be possible though) you can isolate them, with probably
the following modifications:

\begin{itemize}
\item{$\bullet$} A call to {\it error} must be replaced with 
a simple {\it printf} followed by {\it exit}, since these are
supposed to be fatal errors, and the program should abort.

\item{$\bullet$} A call to {\it warning} must be replaced with 
a simple {\it printf}. no need to exit the program though.

\item{$\bullet$} A call to {\it dprintf} must be replaced with 
a simple {\it printf}, by leaving of the first (integer) argument,
or by totally leaving it out of the code. It is merely a debug-type 
printf, and could also be used with an {\tt \#ifdef debug}.

\end{itemize}

These are all small printf-like routines.

\section{Exporting}
\mylabel{as:export}

There are a few automated ways to make full and partial
(incremental) export versions of NEMO. Full exports are
described here, partial and incremental backups in the
next section. \index{export, full}
Again, it's best to be in NEMO's root directory:
\begin{verbatim}
    % cd $NEMO
\end{verbatim}
The first possibility is to make a tar-file of the 
barebone source:\footnote{To exclude hidden dot-files, issue a
``{\tt make purge}'' before creating exporting tarfiles}
\index{tar, export}
\begin{verbatim}
    % make tarfile [FILE=nemo.tar] [SRC="src usr"]
\end{verbatim}
         
With the optional FILE keyword (note the upper case!) you can save
the tar-file in a directory in which you have write permission, just
in case your name is not {\tt nemo}.  The default 
for {\tt FILE} is {\tt nemo.tar} in the
current directory. 
It is currently about 10 Mb, if both the
{\tt src} and {\tt usr} tree are exported.

It is adviced to compress\index{compress} the tar file, since  
in this case it saves a considerable fraction:
\begin{verbatim}
    % compress nemo.tar
\end{verbatim}

The second possibility:
\begin{verbatim}
    % make tarbackup [BACKUP=nemo-fullbck.tar]
\end{verbatim}
creates a complete tar backup image on a diskfile. this diskfile can
be any device (file/tape) using the 
optional backup  keyword. Note that this tar image
will be very large, and is only handy to copy NEMO between identical
systems, or make a backup of the whole system, since also
binary files are backed up.

\section{Small updates: tar import and export}
\mylabel{as:update}

The directory {\tt \$NEMO/adm} is the working directory from where
incremental updates are processed.  A nightly script can be set up by
the NEMO system manager to export new files to all other relevant sites
which run NEMO on a regular basis.  All new files are bundled together
in a tar file, to preserve directory structure and time stamps.  These
tar files are then sent to the central site or other sites.
Similarly, other sites may
export their tar files with updated files to your local site.  These
imported tarfiles are normally located in {\tt \$NEMO/adm/import}. 
Scripts {\tt import.csh} and {\tt export.csh} handle the traffic between
tar files and the local directory.  See also {\it import(8NEMO)} and
{\it export(8NEMO)}.  \index{export, partial} Note that inter-site
collisions are not always handled properly.

\subsection{make tarfile}

If you know what to export, you can also use the {\tt tarfile} target in the
root install {\tt Makefile}. 
For example, in the case where NEMO has to be ported 
to a system with not that much filespace, you could start by copying the
kernel first on the originator (as shown by the {\tt 1\%} prompt below) 
into a tarfile,
the tarfile is then copied to the destination (as shown by the {\tt 2\%}
prompt below) and
extracted out there, after which installation can start:

\begin{verbatim}
    1% cd $NEMO
    1% make tarfile ASCIIDIRS="inc src/kernel"
   ....
    2% cd $NEMO
    2% tar xv nemo.tar
    2% cd $NEMO/src/kernel
    2% make install
\end{verbatim}
This would copy all default ASCIIFILES, as specified in the root
{\tt Makefile}, and your choice of ASCIIDIRS into the tar file
(this bare minimum NEMO kernel is currently about 2Mb).

After having gone through installation, as described earlier, the next
step could be to add the N-body package (or any other one
listed in the {\tt \$NEMO/src} directory) on top of it:

\begin{verbatim}
    1% cd $NEMO
    1% make tarfile ASCIIDIRS=src/nbody ASCIIFILES=
   ....
    2% cd $NEMO
    2% tar xv nemo.tar
    2% cd $NEMO/src/nbody
    2% make install
\end{verbatim}

This tarfile is currently about 1.5Mb.


Typically the installation procedure adds relevant routines to the
NEMO library, after which {\tt mknemo} should be able to install
requested programs from that package.
\begin{table}[ht]
\caption[Small export tar files]{Small export tar file}
\mylabel{t:tarfiles}
\centering
\medskip
\begin{tabular}{|l|l|l|}\hline\hline
%
ASCIIDIRS	& ASCIIFILES	&  Description \\ \hline
%
"inc src/kernel src/scripts" &   -           &  Bootstrap kernel \\
src/tools	& ""            & Various utilities \\
src/nbody	&  ""		&  Nbody integration \\
src/image	&  ""		&  Image utilities \\
src/orbit	&  ""		&  Orbit utilities \\
&& \\ \hline
\end{tabular}
\end{table}

\subsection{WEB maintenance}

\section{Maintenance} 
\mylabel{as:maint}


The directory {\tt \$NEMO/usr/nemo/maint} contains a number of shell scripts
and C-programs which may be useful to do a few consistency checks, cross
correlate two versions of nemo etc. see local documentation, as this is
still part of the old NEMO V1.x release.

In the directory {\tt \$NEMO/adm} a few administrative files
will grow over time. They may have to be cleaned up
on regular intervals, or their features need to be shut down.
see also Section ~\ref{as:update}.

\begin{itemize}
\item the file {\tt NEMO.LOG} \index{NEMO.LOG}
is the accumulated log of the nightly script, which exports new
files to remote sites. 

\item the file {\tt Usage} \index{Usage} is a usage
list of all NEMO programs. It can be used to create various statistics
about usage by user, program, machine, time etc. Use at your own
discretion. the {\tt Usage} file is created by the user interface,
{\it getparam(3NEMO)}, and a logical defined flag {\bf USAGELOG}
\index{USAGELOG} in {\tt getparam.c}
needs to turned off. All of NEMO's programs would need to be recompiled
of course.
\end{itemize}




\section{Other Libraries}
\mylabel{as:otherlibs}

Although NEMO can be installed without any additional libraries, a few
common ones 

Many libraries use the {\tt configure} scripts, and notably with the 
{\tt --prefix=} command line to install their include files, libraries
and other ancillary data in standard locations such as
{\tt /usr/local}. For this purpose NEMO now supports installing them
within the NEMO hierarchy, to ensure ease of porting with the same
compiler and environment:

\begin{verbatim}
    configure --prefix=$NEMO/opt
    make
    make install
\end{verbatim}


\begin{center}
\begin{table}[h]
\caption{Optional Libraries NEMO can use}
\begin{tabular}{||l|l|l||}
\hline
library	  &  version       &   comments \\
\hline
pgplot      & 5.22 (carma)      &  yapp graphics library \\
plplot      & 5.x               &  yapp graphics library \\
hdf         & 4.x               &  data I/O library (e.g. hdfgrid, sdsfits) \\
gsl         & 1.4               & Gnu Scientific Library - optional \\
cfitsio     & 5.x               & fits I/O library - optional \\

\hline


\end{tabular}
\label{t:libraries}
\end{table}
\end{center}

