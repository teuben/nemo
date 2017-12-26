%    -*- LaTeX -*-     for inclusion in nemo.tex
%
% Appendix: The user interface - references
\myfile{aiface.inc}
\mylabel{a:iface}

This Appendix overviews the basic command-line interface of NEMO programs. 

Every NEMO program accepts input through a user supplied parameter list
of {\it 'keyword=value'} arguments. In addition to these 
\index{keywords, program} {\bf program keywords},\index{program keywords}
there are a number of globally defined {\bf system 
keywords},\index{system keywords}
known to every NEMO program. \index{keywords, system}

\section{Program keywords}

Program keywords are unique to a program, and need to be
looked up in the online manual page or by using the 
{\tt help=} system keyword (dubbed the inline
help). Parsing of
``values'' is usually done, though sometimes primitive. Program
keywords also have the ability to read the value(s) of a keyword from a file
through the {\tt keyword=@file} construct. This is called the 
{\bf include keyword file}, and is very handy for long keyword values,
not having to escape shell characters etc.
\index{include, keyword file} \index{@, include keyword file} 
% Note that include keyword files do not nest, and that newlines are replaced
% by blanks.

\section{System keywords}

The 'hidden' system keywords, although overridden by 
any program defined counterpart, can also be set by an 
equivalent environment variable (in upper case), and are:

\begin{description}
\item[help=] Sets the help level to a program. As with all
system keywords, their value can be fixed for
a session by setting the appropriate
environment variable in upper case, {\it e.g.} {\tt "setenv HELP 5"}.
\index{HELP, environment} \index{help=, system keyword}
By using the keyword form, the value of the environment variable 
will be ignored.

The individual help levels are numeric and add up to combine
functionality, and are hence powers of 2:
\begin{itemize}
\item[{\bf 1 }] Remembers previous usage of a program, by
maintaining a keyword file from program to program.  These files are
normally stored in the current directory, but can optionally be stored
in one common directory if the environment variable 
{\bf NEMODEF}\footnote{mirtool also uses this environment variable} is
set.  The keyword files have the name {{\it "progname"}{\bf.def}},
{\it e.g.} {\tt snapshot.def}\footnote{This may result in long
filenames, Unix SYS5 allows only 14 characters - a different solution is
needed here}.  When using this lowest help-level it is still possible to
use UNIX I/O redirection.  This help level
reads, as well as writes the keyword file during the program execution;
hence the user needs both 
read and write permission in the keyword directory.  As can also
be seen, programs cannot run in parallel while using this help-level: they
might compete for the same keyword file.
Within the simple commandline interface it is not possible
to maintain a global keyword database, as is {\it e.g.}  the case in AIPS;
you would have to use the {\tt miriad} shell.

\item[{\bf 2 }] prompts the user for a (new) value for every
keyword; it shows the default (old) value on the prompt line, which can
then be edited.  
It is not possible to combine this level with UNIX I/O redirection. 
By combining the previous helplevel with this one, previous
values and modified ones are maintained in a keyword file.

\item[{\bf 4 }] provides a simple fullscreen menu interface, by having
the user edit the keyword file. The environment variable
{\bf EDITOR} can be used to set any other editor than good old 
{\it vi(1)}.
It is not possible to combine this level with UNIX I/O redirection. 

\item[{\bf 8,16,... }] although not processed, it is reserved for the
next levels of menu interface.


Example: ``{\tt help=3}'' will remember old keywords in a local keyword file,
prompt you with new values, and puts the new values in the keyword file
for the next time.  The ``{\tt help=5}'' option happen to be 
somewhat similar to the way {\tt AIPS} and {\tt IRAF} appear to the user. 

Help levels can also include an alpha-string, which generally display
the values of the keyword, their default values or their help strings.

\item[{\bf ? }]
    lists all these options, as a reminder. It also displays the
    version \index{version, user interface} of the 
    {\tt getparam} user interface package.

\item[{\bf h }]
    list all the keywords, plus a help string what the keywords does/expects.
This is really what we call the inline manual or inline 
help. \index{inline, help} \index{manual, inline} \index{help, inline}

\item[{\bf a }]
    list all arguments in the form {\it keyword=value}.

\item[{\bf p,k }]
    list parameters (keywords) of all arguments in the form {\it keyword}.

\item[{\bf d,v }]
    list defaults (values) of all arguments in the form {\it value}.

\item[{\bf n }]
    add a newline to every {\it keyword/value} string on output.
In this way a keyword file could be build manually by redirecting this
output.

\item[{\bf t }]
    output a documentation file according to the
\%N,\%A specifications \index{mirtool} of 
{\tt miriad}\footnote{Both {\tt mirtool} and {\tt miriad} need such a doc-file
\index{doc file, miriad} to lookup keywords and supply help}.
Is mainly intended to be used by scripts such as {\tt mktool}. 
The procedure in NEMO to update a {\tt .doc} file would be: \index{miriad}
\begin{verbatim}
         % program help=t > $NEMODOC/program.doc
    or:
         % mktool program 
\end{verbatim}
if the script {\tt mktool} has been\index{mktool}
installed\footnote{Obviously this is priviliges to NEMO superusers only}

\item[{\bf q }]
    quit, do not start program. Useful when the helpstring contains
options to print.
% is this really correct, perhaps we need a 'g' options, to force 'go'

Example: {\bf key=val help=1q} redefines a keyword in the keywordfile,
but does not run the program. This is also a way to 'repair' a keyword
file, when the program has been updated with new keywords.
{\bf key=val help=1aq} redefines the keyword,
shows the results but does still not run the program. 
Finally, {\bf key=val help=1a} redefines a keyword, shows
the result and then runs the program.
\end{itemize}

\item[host=] Runs the program on a remote host.  \index{remote execution}
It depends on the \index{host=, system keyword}
implementation of software on local as well as remote host if and how
this option works.  Among SUN systems the {\tt rsh} command is used, and
assumes a shared disk with the same absolute pathname (NFS).  Future
implementations will have to use more sophisticated RPC (Remote
Procedure Call) or X11 interfaces for distributed networking.  No
environment variable is used here. 

\item[debug=] Changes the debug output level.  
\index{debug=, system keyword} \index{DEBUG, environment}
The higher the debug
level, the more output can appear on the standard error output device
{\tt stderr}.  The default value is either 0 or the value set by the
{\bf DEBUG} environment variable.  The use of the {\tt debug=} keyword
will override your default setting.  A value of '0' for debug 
may still show some warning messages.  Setting debug to 
-1 will prevent even those warning/debug messages.  Legal values are 0
through 9.  Values of {\bf DEBUG} higher than 9 are for system experts
usage only.  You may get some weird screen output. Values larger than
5 cause an error to coredump, which can then be used with debug utilities
like {\it abd(1)} and {\it dbx(1)}.

\item[error=] Specifies how many times the fatal error routine can be
bypassed\index{error=, system keyword}. The {\bf ERROR} environment
variable can also be set for this. The default, if neither of them
present, is 0.

\item[yapp=] Defines the device to which graphics output is send. 
Currently only interpreted for a limited number of yapp devices.  
\index{YAPP, environment} \index{yapp=, system keyword} Some yapp
devices do not even listen to this keyword.  Check {\it yapp(5NEMO)} or
your local NEMO guru which one is installed.  The default device is
either 0 or the value set by the {\bf YAPP} environment variable.  

\subsection{Yapp\_mongo}

\index{yapp, mongo}
Valid devices in {\it yapp\_mongo} are (see also 
MONGO-87\footnote{MONGO is a copyrighted program 
by John Tonry\index{Tonry J.}}
users manual) numbers {\tt 1..6, -6..-1}:

\small\begin{verbatim}
1.  Retrographics 640    -1.  Versatec (portrait)
2.  DEC VT125            -2.  Versatec (landscape)
3.  Tektronix 4010       -3.  Printronix (portrait)
4.  Grinell 270          -4.  Printronix (landscape)
5.  HP 2648A             -5.  Postscript (portrait)
6.  Sun Windows          -6.  Postscript (landscape)
                         -7.  Postscript (square portrait) ***
\end{verbatim}\normalsize

Make sure you have set:
\begin{verbatim}
    setenv MONGOPATH /usr/local/lib/mongo
    setenv MONGOFILES $MONGOPATH/mongofiles.dat
\end{verbatim}
or appropriate private ones, or use whatever your system manager has
provided for.
\subsection{Yapp\_sunview}

\index{yapp,sunview}
For the {\it yapp\_sv} (Sunview) interface the absolute value of the {\tt
yapp=} keyword gives the (square) size in pixels 
of the window created. 
The default is 128. For 
positive value the window will disappear when {\it plclose(3NEMO)}
is called at the termination of the graphics operations (essential for
movie production).  For negative values it remains displayed until a
Carriage Return is given. 

\subsection{Yapp\_pgplot}
\index{yapp, pgplot}
A graphics device in PGPLOT\footnote{PGPLOT is a copyrighted public
domain graphics library written by Tim Pierson\index{Pierson T.}}
is defined by preceding it with a slash, in
the same way as is done in PGPLOT itself.
Optional parameters can be supplied before the slash. The following
list gives an overview of some of the available devices
(your list may be a lot shorter (see {\tt '?'} in list below):

\small\begin{verbatim}
?           Get a list of all currently defined graphics devices
x,y/sun     Sunview, (x,y) = sizes in inches [8,8]
file/ps     Landscape Postscript, file=filename [pgplot.ps]
file/vps    Portrait Postscript, file=filename [pgplot.ps]
/null       Null device
/xwindow    X-windows
/xdisp      pgdisp or figdisp server
/tek4010    Tektronix
/tk4100
/gterm      IRAF Gterm window within sunview
/retro      Retrographics
/file       PGPLOT Metafile
/gf         Graphon
/vt125
/printronix
/peritek
/ws4        Landscape
/tfile
/pk
\end{verbatim}\normalsize

\end{description}

{\it BUGS:} 
When a non-zero help level is used, 
one cannot specify system keywords, other than 
by specifying them through environment variables, 

Also consult manual pages such as {\it getparam(3NEMO)} and {\it
yapp(5NEMO)}.
% in case this description get outdated. 

\subsection{Yapp\_sm}
\index{yapp, sm}

A graphics device in {\it sm}\footnote{SM is a copyrighted program
and subroutine library by Robert Lupton\index{Lupton R.} 
and Patricia Monger\index{Monger P.}} is defined by the same 
text string as you would issue it in interactive mode. The 
{\it sm} command {\tt list devices} will list all currently
compiled device drivers. You need a file {\tt .sm} in your
home directory, or take what your system manager has provided
for.

\newpage
\section{The REVIEW section}

\mylabel{a:review}
By setting the {\bf REVIEW} environment variable a NEMO program is 
always put into the REVIEW \index{REVIEW} \index{interrupting a program}
section just before the start of the actual execution of
the program (the end of the {\it initparam(3NEMO)} routine). 
This functionality is quite similar to using the helplevel
{\tt help=4} (see previous Section).

A NEMO program can also be interrupted, using the quit signal
(see {\it signal(2)}), \index{signal(2)} into the 
REVIEW section, although the program must be adapted to get
keyword information through {\it getparam(3NEMO)} and not through it's
own local database, in order for modified keywords to take effect.
This does not hold for the system keywords, whose new value is always
correctly interpreted by the program.

In the REVIEW section the prompt is {\bf ``REVIEW''} and
the following commands are understood:

\begin{description}
\item{{\bf exit, quit, end}}
  Exit the program (ungracefully).

\item{{\bf stop}}
  Gracefully end the program, but first goes through {\tt finiparam()} (see
{\it getparam(3NEMO)}) to update the keyword file if the
helplevel includes 1.

\item{{\bf set [key=[value]]}}
  Set a new value for a program keyword ({\tt set key=value}), where
{\tt value} may also be blank, or display the contents of a 
program keyword ({\tt set key}).

\item{{\bf show key}}
 Show the value of a program keyword.

\item{{\bf keys}}
 Show the values of all program keywords.

\item{{\bf syskeys}}
 Show the values of all system keywords.

\item{{\bf set syskey[=value]}}
 Set a new value for a system keyword {\tt set syskey=value}
or display its current contents {\tt set syskey}. 

\item{{\bf time}}
 Show the cputime (in minutes) used so far. \index{cputime}

\item{{\bf !cmd}}
 Pass a command {\tt cmd} to the shell.

%\item{{\bf HOTKEY}}
% Abort the program badly, possibly generating a coredump.

\item{{\bf go,continue}}
 Continue execution of the program.

\item{{\bf version}}
 Display version of {\tt initparam()} compiled into program.

\item{{\bf ?, help}}
 Displays all commands and their format.

\end{description}

When the system keyword {\tt debug} is non-zero, the ``REVIEW'' prompt also
includes the process identification number of the process.
