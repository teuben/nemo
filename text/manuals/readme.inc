%% This file is to be included by latex in nemo.tex
%
%  26-mar-94    for V2.2
%
% In the Opening statements: README 

\myfile{readme.inc}

In order to prevent disappointment and frustration this section 
should be carefully read and understood before you endeavor into 
NEMO.

You can use the list below as a checklist, it will help you to decide if NEMO is
really the right way for you to solve your particular problem.  In the
end it may be helpful to have the manual pages {\it programs(8NEMO)}
and/or \index{index, programs} {\it index(1NEMO)} for names of
specific programs and utilities that you may want to use. The items in
the list below are in approximate order of importance, the most
important ones listed first. The further down the list you come, the
easier it will be for you to get along with the NEMO package.

\begin{enumerate}


\item[$\bullet$]

understand NEMO is nothing more than an extension of the UNIX
environment, and that you know how to modify your default UNIX
environment to be able to run NEMO programs (See 
Appendix~\ref{a:setup}). This generally makes it relatively easy 
to use NEMO in other environments and packages.

\item[$\bullet$]

understand the basic user interface: programs have a '{\bf keyword=value}'
argument structure, and that there are '{\bf program}' and '{\bf system}' 
keywords. System keywords can generally be set fixed by an equivalent
environment variable 
(Chapter~\ref{c:iface}, Appendix~\ref{a:iface})

\item[$\bullet$]

realize NEMO is mainly an N-body package, with utilities to create
stellar systems, evolving them, and a large variety of analysis and
display programs.   However, for a number of problems an N-body system
can be effectively and efficiently used to simulate a completely
different physical situation.  There are interfaces to create/convert
data in/to image format, and export them in FITS format
(Chapter~\ref{c:examples}).  A small amount of orbit and table
manipulation utilities are also available.  (See also
Chapter~\ref{c:filestr}). 

\item[$\bullet$]

understand the basic workings of a UNIX shell; how to write shell scripts
in which programs are combined in a modular way, and to automate 
procedures. For this reason NEMO programs are 
also geared towards batch usage. 
(see any introductory UNIX book and Appendix~\ref{a:shell})

\item[$\bullet$]

understand that data is stored in binary format and in a general
hierarchical/ structured file format. The data can always be viewed in
human readable form with the program  {\it tsf(1NEMO)}.  Depending on
the kind of data, the format is structured in a specific way, e.g.
snapshots, images, orbits. Only problem specific programs can read
such datafiles (See Chapter~\ref{c:filestr}).

\item[$\bullet$]

you know about the graphics interface {\tt yapp} (see also\index{yapp}
Chapter~\ref{c:graphics}). Make  sure you understand the local NEMO
implementation how to get graphics  output to for example
screen and printer.
Sometimes they are hidden in one device  driver, and the {\bf yapp=}
system keyword is then used to select the device ({\it e.g.} the {\it
mongo}\index{mongo} and {\it pgplot}\index{pgplot} yapp drivers),  or
programs are linked with different device drivers. In this case, the 
program name will have a underscore and the device name appended, such
as {\it snapplot, snapplot\_cg} and {\it snapplot\_ps}.
(See Chapter~\ref{c:graphics}).

\item[$\bullet$]

realize that NEMO can be extended to suit your own needs. It is  fairly
easy to write new programs or even define a new data structure (as was
done for  snapshots, images and orbits). Do realize however that the
main language  used is C and support for linking FORTRAN and C is
minimal and system dependant. (See Chapter~\ref{c:progr}).

\item[$\bullet$]

you realize what it means to install NEMO if this has not been done yet.
You need to carefully read Appendix~\ref{a:install} in detail; any
experience you have had with the UNIX shell, and utilities like 
``{\it make}''  and ``{\it autoconf}'' will be useful.

\end{enumerate}


\section*{More information}

The manual that you are reading here is in LaTeX format and should be in
{\tt \$NEMO/text/manuals}\footnote{The full NEMO directory structure is
outlined in Appendix~\ref{a:dirs}}.
Other information on NEMO can be found in the following places:

\begin{description}

\item[$\bullet$] manual pages, in {\tt \$NEMO/man} and below, using the
UNIX {\it man} or NEMO {\it nemoman}/{\it manlaser}
commands\index{nemoman}\index{manlaser}. For a graphical
interface programs like\index{manual pages, html}
{\it tkman}\index{tkman} and {\it gman}\index{gman} can be very useful
(hypertext) browsers of manual page. We also automatically maintain a
fully hyperlinked version formatted in html.\index{html, manual pages}

\item[$\bullet$] The {\tt README} file in NEMO's home directory will
always contain ``last minute'' information that has perhaps not found
its way into this manual yet. 

\item[$\bullet$] The file {\tt faq.tex} in {\tt \$NEMO/text} lists {\it
Frequently Asked Questions}, and possibly other files in that directory. 

\item[$\bullet$] A summary sheet {\tt summary.tex} in {\tt \$NEMO/text}
overviews NEMO's usage. 

\item[$\bullet$] {\tt survival.texi}, a short texinfo\index{texinfo}
introduction with basic survival techniques to get around.

\item[$\bullet$] The NEMO home page, currently at
\htmladdnormallink{{\tt http://www.astro.umd.edu/nemo/}}
{http://www.astro.umd.edu/nemo/}, contains many
pointers to a variety of helpful pages, and sometimes
contradicts information written down in this manual.

\item[$\bullet$] The ManyBody compendium, putting together an introduction
to NEMO, starlab and other N-body software, originally written for
the first N-body School (Strasbourg, 2004)\index{manybody}\index{Strasbourg}.



\end{description}


