%% This file is to be included by latex in nemo.tex
%
% Chapter:  Filestructure
\myfile{filestr.inc}
\mylabel{c:filestr}

This chapter gives an overview of the file structure of persistent
data\footnote{Programmers are always free to choose any format they like
in memory - this is usally hidden from the users.  What we mean here is
the disk format. The popular memory (object) models,
and how they interact with persistent data on disk, are discussed in
Chapter~\ref{c:progr}}.
Most of the data handled by NEMO is in the form of a
specially designed\footnote{note this was back in 1986, well before XML was conceived}
XML-like binary format, 
although exceptions like ASCII files/tables will also be discussed. 
Ample examples
illustrate creation, manipulation and data transfer.
We close this chapter with a few examples of function descriptions, 
a dataformat that make use of
the native object file format of your operating system ({\it a.out(5)}
and {\it dlopen(3)}).

\section{Binary Structured Files}

Most of the data files used by NEMO share a common low level binary file
structure, which can be viewed as a sequence of tagged data items.  Special
symbols are defined to group these items hierarchically into sets.  Data items
are typically scalar values or homogeneous arrays constructed from
elementary C data types, but the programmer can also add more complex
structures, such as C's {\tt struct} structure definition, or any user
defined data structure. In this last case tagging by type is not 
possible anymore, and support for a machine independent format
is not guaranteed. Using such constructs is not recommended if the
data needs to be portable accross platforms.

The hierarchical structure of a binary file in this general format can
be viewed in human-readable format at the terminal using a special 
program, {\tt tsf} 
\index{tsf(1)} ("{\it type structured file}").
Its counterpart, {\tt rsf} \index{rsf(1)}
("{\it read structured file}"),
converts such human-readable files (in that special ASCII Structured
File format, or\index{ASF, Ascii Structured File} ASF)
into binary structured files 
(BSF)\index{BSF, Binary Structured File}.
In principle it is hence
possible to transfer data files between different types of computers
using {\tt rsf} and {\tt tsf} (see examples in Section~\ref{s:exch-data}).
\footnote{Note however that currently NEMO's binary files have some limited support
for machine-independancies, e.g. simple endian swap.
Some portability notes about this can be found in Appendix~\ref{sifs}}

Let us start with a small example: With the NEMO
program {\tt mkplummer} we first create an
N-body realization of a spherical Plummer model:\index{mkplummer(1)}
\index{Plummer, model}

\small\begin{verbatim}
    1% mkplummer i001.dat 1024
\end{verbatim}\normalsize

Note that we made use of the shortcut that {\tt out=} and {\tt nbody=}
are the first two program \index{keywords, program} keywords, and they
were assigned their value by position rather than by associated name.
We can now display the contents of the binary file {\tt i001.dat} with {\tt
tsf}\index{tsf}:

\small\begin{verbatim}
    2% tsf i001.dat

    char Headline[33] "set_xrandom: seed used 706921861"
    char History[36] "mkplummer i001.dat 1024 VERSION=2.5"
    set SnapShot                                                            
      set Parameters                                                        
        int Nobj 01750
        double Time 0.00000                                                 
      tes                                                                   
      set Particles                                                         
        int CoordSystem 0201402                                             
        double Mass[1024] 0.00195313 0.00195313 0.00195313
          0.00195313 0.00195313 0.00195313 0.00195313 0.00195313
          0.00195313 0.00195313 0.00195313 0.00195313 0.00195313
          0.00195313 0.00195313 0.00195313 0.00195313 0.00195313
          . . .
        double PhaseSpace[1024][2][3] 4.92932 0.425103 -0.474249
          0.342025 -0.112242 4.60796 -0.00388599 -0.389558 -0.958787
          0.220561 0.213904 3.47561 0.0176012 1.22146 -0.903484
          -0.705422 4.26963 -0.263561 1.04382 -0.199518 -0.480749
          . . .                                                             
      tes                                                                   
    tes                                                                     
\end{verbatim}\normalsize

This is an example of a data-file from the N-body group, and consists of
a single {\it snapshot} \index{snapshot} at time=0.0.  This snapshot,
with 1024 bodies with double precision masses and full 6 dimensional
phase space coordinates, totals 57606 bytes, whereas a straight dump of
only the essential information would have been 57344 bytes, a mere 0.5\%
overhead.  The overhead will be larger with small amounts of data, {\it
e.g.} diagnostics in an N-body simulation, or small N-body snapshots. 

Besides some parameters in the '{\tt Parameters}' set, it consists
of a '{\tt Particles}' set, where (along the type of coordinate system)
all the masses and phase space coordinates of all particles
are defined. Note the convention of integers starting with
a '{\tt 0}' in octal\index{octal}
representation. This is done for portability
reasons.

A remark about online help\index{help, online}:
NEMO also uses the {\it man(5)} format\index{man, manual pages}
for more detailed online help, 
although the inline help (system {\tt help=} keyword)
\index{help, inline} is most of the times sufficient enough
to remind a novice user of the keywords and their meaning.
The {\tt man} command is a last resort, if more detailed information
and examples are needed. 

\begin{verbatim}
    3% man tsf
\end{verbatim}

Note that, since the online manual page is a different file from the
source code, information in the manual page can easily get outdated, and
the \index{help, inline}inline ({\tt help=}) help, although very brief,
is more likely to be up to date since it is generated from the source
code (executable) itself:

\begin{verbatim}
    4% tsf help=h
    in               : input file name  [???]
    maxprec          : print nums with max precision  [false]
    maxline          : max lines per item  [4]
    allline          : print all lines (overrides maxline)  [false]
    indent           : indentation of compound items  [2]
    margin           : righthand margin  [72]
    item             : Select specific item []
    xml              : output data in XML format? (experimental) [f]
    octal            : Force integer output in octal again? [f]
    VERSION          : 29-aug-02 PJT  [3.1]
\end{verbatim}

\section{Pipes}

In the UNIX operating system pipes\index{pipes} can be very
effectively used to pass information from one process to 
another. One of the well known textbook examples is how one
gets a list of misspelled (or unknown) words from a document:

\begin{verbatim}
    % spell file | sort | uniq | more
\end{verbatim}

NEMO programs can also pass data via UNIX pipes, although with a
slightly different syntax: a dataset that is going to be part of a pipe
(either input or output) has to be designated with  the {\tt -}
(``dash'') symbol for their filename.
Also, and this is very important, the receiving task
at the other end of the pipe should get data from only one source
of course. If
the task at the sending end of the pipe wants to send binary data over
that pipe, but in addition the same task would also write ``normal''
standard
output, the pipe would be corrupted with two incompatible sources of
data. An example of this is the program 
{\tt snapcenter}\index{snapcenter}. The keyword {\tt report} must be set to
{\tt false} instead, which is actually the default now.
So, for example, the output of a previous N-body
integration is re-centered on it's center of mass, and subsequently
rectified and stacked into a single image as follows:

\begin{verbatim}
    % snapcenter r001.dat . report=t  | tabplot - 0 1,2,3
    
    % snapcenter r001.dat - report=f       |\
        snaprect - - 'weight=-phi*phi*phi' |\
        snapgrid - r001.sum stack=t
\end{verbatim}

If the keyword {\tt report=f} would not have been set properly,
{\tt snaprect}\index{snaprect}
would not have been able to process it's convoluted
input. Some other examples
are discussed in Section~\ref{ss:data}.

\section{History of Data Reduction}

Most programs\footnote{notable exceptions are basic programs
like {\tt tsf}, {\tt rsf}, {\tt csf} and {\tt hisf}}
in NEMO will automatically keep track of the history of
their data-files in a self-describing and self-documenting
way.\index{history} If a program modifies an input file and produces an
output file, it will prepend the
command-line with which it was invoked to its data history.  The
data history is normally located at the beginning of a data file. 
Comments entered using the frequently used program keyword {\tt
headline=} will also appear in the history section of your data file. 
\index{HISTORY, environment}


A utility, {\tt hisf} \index{hisf(1)}
can be used to display the history of a data-file. 
This utility can also be used to create a pure history file (without any
data) by using the optional {\tt out=} and {\tt text=} keywords.  Of
course {\tt tsf} \index{tsf(1)}
could also be used by scanning its output for the string\index{headline}
{\tt History} or {\tt Headline}:

\begin{verbatim}
    5% tsf r001.dat | grep History
\end{verbatim}

which shows that {\tt tsf}, together with it's counterpart {\tt rsf} has
virtually
the same functionality as {\tt hisf}. 
\footnote{HISTORIC NOTE: To prevent data files having a history written
into them an environment variable {\bf HISTORY} must be set to 0.  This
dates back from older times when not all programs could  properly handle
data files with embedded history items properly. Also note there is no
associated system keyword with {\bf HISTORY}. It is expected that this
feature will disappear, {\it i.e.} history is always forcefully written
into the data files, unless the user interface ({\tt
getparam.o}\index{getparam} in {\tt libnemo.a}\index{libnemo}) was
explicitly compiled with  \index{user interface, getparam} the HISTORY
disabled.}

% Again, use {\it man(1)} to get more information on {\tt hisf}. 



\section{Table format}

\mylabel{s-table}
Many programs are capable of producing standard output in (ASCII)
tabular format.  \index{table, files}
The output can be gathered into a file using
standard UNIX I/O redirection.  In the example \index{radprof(1)}

\begin{verbatim}
    6% radprof r001.dat tab=true > r001.tab
\end{verbatim}

the file {\tt r001.tab} will contain (amongst others) columns with
surface density and radius from the snapshot {\tt r001.dat}.  These
(ASCII) 'table' files can be used by various programs for further
display and analysis.  NEMO also has a few programs for this purpose
available ({\it e.g..} {\tt tabhist} \index{tabhist(1)}
for analysis and histogram
plotting, {\tt tablsqfit} \index{tablsqfit(1)} 
for checking correlations between two columns and
{\tt tabmath} \index{tabmath(1)} for general table handling).  
% Commercial packages such as (S)Mongo \index{mongo} \index{sm}
% and Lotus 1-2-3 \index{Lotus 1-2-3}
% may also come in very handy here and easily read table files.
The manual 
pages of the relevant NEMO programs should inform you how to get nice
tabular output, but sometimes it is also necessary to write a shell/awk
\index{awk, script}\index{shell, script}
script or parser to do the job. {\it Note: the {\tt tab=} keyword
hints at the existence of such features}.

A usefull (public domain) program {\it redir(1NEMO)}\index{redirection}
has been included in NEMO\footnote{see also the {\tt tpipe} tool}
to be\index{tpipe}
able split the two standard UNIX output channels {\it stdout} and
{\it stderr} to separate files. \index{stdout} \index{stderr}
\index{redir}

\begin{verbatim}
    7% redir -e debug.out tsf r001.dat debug=2
\end{verbatim}

would run the {\tt tsf} command, but redirecting the
{\it stderr} standard error output to a file {\tt stderr.out}. There are
ways in the C-shell to do the same thing, but they are
clumsy and hard to remember. In the bourne 
shell ({\tt /bin/sh})\index{bourne shell} this is accomplished much
easier:

\begin{verbatim}
    7$ tsf r001.dat debug=2 2>debug.out
\end{verbatim} 

One last word of caution regarding tables: tables can also be used
very effectively in pipes, for example take the first example,
and pipe the output into {\tt tabplot} to get a quick look 
at the profile:

\begin{verbatim}
    8% snapprint r001.dat r | tabhist - 
\end{verbatim}

If the snapshot contains more than 10,000 points, {\tt tabhist} cannot
read the remainer of the file, since the default maximum number
of libes for reading from pipes\index{pipe, limiting lines}
is set by a keyword {\tt nmax=10000}. To properly read all lines, you
have to know (or estimate) the number of lines. In 
the other case where the input is a regular file, table programs
are always able to find the correct amount to allocate for their
internal buffers by scanning over the file once. For very large tables
this does introduce a little extra overhead.

\section{Dynamically Loadable Functions}

A very peculiar data file format encountered in NEMO is that of the 
function descriptors. They present themselves to the user through
one or more keywords, and in reality point to a compiled\index{loadobj}
piece of code that will get loaded by NEMO (using {\tt loadobj(3NEMO)}).
We currently have 4 of these in NEMO:

\subsection{Potential Descriptors}


The potential
\index{potential descriptor} \index{orbit} descriptor is used in orbit
calculations and a few N-body programs.  These are actually binary
object files (hence extremely system dependent!!), and 
used by the dynamic object loader \index{dynamic loader, loadobj}
during runtime. Potentials are 
supplied to NEMO programs as an input variable ({\it i.e.} a set of 
keywords\footnote{Normally called {\tt potname=, potpars=} 
and {\tt potfile=}, but see also {\tt rotcurves}}).
For this, a mechanism is needed to dynamically load 
the code which calculates the potential. This is done by a
dynamic object loader that comes with NEMO. 
If a program needs a potential, and it is present in the
default repository ({\tt \$POTPATH} 
or \index{POTPATH, environment}
{\tt \$NEMOOBJ/potential}), it is
directly loaded into memory by this dynamic object loader. 
If only a source file is present, {\it
e.g.} in the current directory, it is compiled on the fly \index{compile
on the fly} and then loaded.  The source code can be written
in C or FORTRAN.  Rules and more information
can be found in {\it potential(3NEMO)} and {\it
potential(5NEMO)} The program {\it potlist(1NEMO)} \index{potlist(1)}
can be used to test potential descriptors. See Section~\ref{s:potential}
for examples.

\subsection{Bodytrans Functions}

Another family of object files used by the dynamic
object loader are the {\it bodytrans(5NEMO)} functions. These were
actually the first one of this kind introduced in NEMO.
They are functions generated from expressions containing body-variables
(mass, position, potential, time, ordinal number etc.).  They frequently occur
in programs where it is desirable to have an arbitrary
expression of body variables
{\it e.g.}  plotting and printing programs, sorting program etc.
Expressions which are not in the standard repository (currently 
{\tt \$BTRPATH} 
or \index{BTRPATH, environment}
{\tt \$NEMOOBJ/bodytrans}) will 
be generated on the fly and saved for later use. 
\index{bodytrans(1)} The program {\it bodytrans(1NEMO)} is available
to test and save new expressions. Examples are given in 
Section~\ref{s-dispanal}, a table of the 
precompiled ones are in Table~\ref{t:bodytrans}.

\subsection{Nonlinear Least Squares Fitting Functions}

The program {\tt tabnllsqfit} can fit (linear or non-linear, depending
on the parameters) a function to a set of datapoints from an ASCII table.
The keyword {\tt fit=} describes the model (e.g. a line, plane, gaussian, circle,
etc.), of which a few common ones have been pre-compiled with the program.
In that sense this is different from the previous two function descriptors,
which always get loaded from a directory with precompiled object files.
The keyword {\tt load=} can be used to feed a user defined function to
this program. The manual page has a lot more details.

\subsection{Rotation Curves Fitting Functions}

Very similar to the Nonlinear Least Squares Fitting Functions are the
Rotation Curves Fitting Functions, except they are peculiar to the
1- and 2-dimensional rotation curves one find in galaxies as the 
result of a projected circular streaming model. The program
{\tt rotcurshape} is the only program that uses these functions, the
manual page has a lot more details.