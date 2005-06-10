%% This file is to be included by latex in nemo.tex
%
% Chapter: NEMO Programmers Guide: an Introduction
\myfile{progr.inc}
\mylabel{c:progr}

In this chapter an introduction\footnote{Based on an original report
{\it ``NEMO: Elementary Mechanics Observatory''} by Joshua Barnes} is
given\index{NEMO, the name, another origin} how to write programs within
the NEMO environment. 
 

To the application programmer NEMO consists of a set of macro
definitions and object libraries, designed for numerical work in general
and stellar dynamics in particular.  A basic knowledge how to program in
C, and use the local operating system to get the job done, is assumed. 
If you know how to work with {\tt Makefile}'s, even the better. 

After reviewing how the NEMO environment should be present,
Section~\ref{s:macro} describes some of the available macro packages. 
How an example program is written, compiled and used in NEMO is shown in
Section~\ref{s:example}.   In Section~\ref{s:cplusplus} we will show how
you can write programs in C++, and still use the NEMO libraries.
Finally, Section \ref{s:fortran} deals with those people who insist on
using FORTRAN, 

\section{The NEMO Programming Environment}

The modifications necessary to your UNIX environment in order to
access NEMO are 
extensively described in Appendix~\ref{a:setup}.
This not only applies
to a user, but also to the application programmer, although for the
latter the
{\bf static} environment is to be preferred here.
% You may even want to modify the {\tt NEMORC} \index{NEMORC, file}
% startup file and use your own private version,
% although this is not a recommended practice.
Relevant environment variables are described in Appendix~\ref{s-envvar}

In summary, the essential changes to your environment consist of three
additions to your local {\tt .cshrc}\index{.cshrc} startup file (or its
equivalent if you don't use the {\tt csh}-shell):

\begin{enumerate}

\item set the environment variable {\bf NEMO}\index{NEMO, environment}
 to the location of the root directory of NEMO, and
 the {\bf NEMOHOST}\index{NEMOHOST, environment}
 environment variable.

\item source the startup file {\tt \$NEMO/NEMORC}\index{NEMORC, file}.

\item add {\tt \$NEMOBIN}\index{NEMOBIN, environment}
to your search path in a convenient location (but before the system
directories where the C compiler is located).

\end{enumerate}

Although the environment variable {\bf NEMOHOST}\index{NEMOHOST,
environment} is used to allow the package to be shared across a number
of different architectures ({\it e.g.} Sun3's and Sun4's), it is 
recommended to use this environment
variable 
in a single-architecture environment 

Once the NEMO environment is merged into your UNIX environment,
most programs can be compiled as follows:

\begin{verbatim}
    cc -o try try.c -lnemo -lm
\end{verbatim}


\section{The NEMO Macro Packages}
\mylabel{s:macro}
We will describe a few of the most frequently used macro packages
available to the programmer.  They reside in the form of header include
files in a directory tree starting at {\tt \$NEMOINC}.  Your application
code would need references like:

\begin{verbatim}
    #include <stdinc.h>
    #include <snapshot/body.h>
    #include <snapshot/get_snap.c>
\end{verbatim}

Some of the macro packages are merely function prototypes, to
facilitate modern C compilers, and 
have associated object code in libraries in {\tt \$NEMOLIB} and programs
need to be linked with the appropriate ones.



\subsection{stdinc.h}

The macro package {\tt stdinc.h} \index{stdinc.h} provides all basic
definitions that ALL of NEMO's code must include as the first include
file.  It also replaces the often used {\tt stdio.h} include file in C
programs.  The {\tt stdinc.h} include file will provide us with a way to
standardize on future expansions, and make code more
machine/implementation independent ({\it e.g.} POSIX.1\index{POSIX.1}).  In
addition, it defines a more logical standard for C notation.  For
example, the normal C practice of using pointers to character for
pointer to byte, or integer for bool, tends to encourage a degree of
sloppy programming, which can be hard to understand at a later date. 

A few of the basic definitions in this package:

\begin{description}

\item[$\bullet$ {\tt NULL}:]
 macro for 0, used to distinguish null characters and null pointers.
This is often already defined by {\tt stdio.h}.\index{NULL}
There is potential trouble when NULL has been set to 
{\tt (void *)0}, \verb+'\0'+ is OK though. Example on IBM's AIX operating
system.

\item[$\bullet$ {\tt bool}:]
typedef for {\tt short int} or {\tt char},
used to specify boolean\index{bool} data. See also next item.
\footnote{The {\it curses} library also defines {\tt bool}, and this
made us change from {\tt short int} to {\tt char}}.

\item[$\bullet$ {\tt TRUE, FALSE}:] 
macros for 1 and 0, respectively, following normal C 
conventions.\index{TRUE}\index{FALSE}

\item[$\bullet$ {\tt byte}:]
typedef for {\tt unsigned char}, used to specify byte-sized data.
\index{byte}

\item[$\bullet$ {\tt string}:]
typedef for {\tt char *}, used to point to strings. Don't use
{\tt string} for pointers you increment, decrement or explicitly
follow (using *); such pointers are really {\tt char *}.
\index{string}

\item[$\bullet$ {\tt real, realptr}:]
typedef for {\tt float} or {\tt double} ({\tt float *} or {\tt double *},
respectively), depending on the use of the {\tt SINGLEPREC} flag. 
The default is {\tt double}. \index{real}\index{realptr}

\item[$\bullet$ {\tt proc, iproc, rproc}:] 
typedefs for pointers to procedures (void functions), integer-valued
functions and real-valued functions respectively. 
% This always confusing issue will become clear later on.
\index{iproc} \index{rproc}\index{proc}

\item[$\bullet$ {\tt local, permanent}:]
macros for {\tt static}. Use {\tt local} when declaring variables or
functions within a file to be local to that file. They will not appear
in the symbol table be usable as external symbols. Use 
{\tt permanent} within  a function, to retain their value upon
subsequent re-entries in that function.\index{local}\index{permanent}

\item[$\bullet$ {\tt PI, TWO\_PI, FOUR\_PI, HALF\_PI, FRTHRD\_PI}:]
macros for $\pi$, $2\pi$, $4\pi$, $\pi/2$ and $4\pi/3$, 
respectively.\index{PI}

\item[$\bullet$ {\tt ABS(x), SGN(x)}:]  
macros for absolute value and sign of {\tt x}, irrespective of the 
type of {\tt x}.\index{ABS}\index{SGN}. Beware of side effects.

\item[$\bullet$ {\tt MAX(x,y), MIN(x,y)}:]
macros for the maximum and minimum of {\tt x,y}, irrespective of the
type of {\tt x,y}.\index{MAX}\index{MIN}
Beware of side effects.

\item[$\bullet$ {\tt stream}:]
typedef for {\tt FILE *}. They are mostly used with the NEMO
functions {\tt stropen} and {\tt strclose}, which are functionally
similar to {\it fopen(3)} and {\it fclose(3)}, plus some added 
NEMO quirks. \index{stream}  (see section~\ref{ss-filestr} below)

\end{description}


\subsection{getparam.h}

The command line syntax described earlier in Chapter~\ref{c:iface}
is implemented by a small set of functions used by all conforming
NEMO programs.  A few function calls generally suffice to get the values
of the input parameters.  A number of more
complex parsing routines are also available, to be discussed in the
next subsection.

First of all, a NEMO program must define which 
{\bf program keywords}\index{program keywords}\index{keywords, program}
it will 
recognize. For this purpose it must define an array of {\tt string}s
with the names and the default values for the keywords, and optionally,
but STRONGLY recommended,
a one line help string for that keyword: \index{getparam.h}

\begin{verbatim}
    #include <stdinc.h>     /*   every NEMO module needs this  */
    #include <getparam.h>   /* needed when user interface used */

    string defv[] = {       /* definitions of the keywords */
        "in=???\n          Input file (a snapshot)",
        "n=10\n            Number of particles to view",
        "VERSION=1.1\n     14-jul-89 - 200th Bastille Day - PJT",
        NULL,
    };

    string usage = "example program";   /* def. of the usage line */

\end{verbatim}
\index{defv, string} \index{help, -string in keywords}
The ``{\it keyword=value}'' and ``{\it help}'' part of the string must be
separated by a newline symbol ({\tt $\backslash$n}). If no 
newline\index{newline, help} is present,
as was the case in earlier releases, no help 
string is available.\footnote{ZENO\index{ZENO} uses a different technique: ...}
\index{help, keywords}
The {\tt 'help=h'} command
line option displays the ``{\it help}'' part of string during
execution of the program for quick inline reference.\index{inline, help}
The ``{\it usage}'' part defines a string that is used as a one
line reminder what the program does. It's used by the various
invocations of the user interface.\index{usage, string}

The first thing a NEMO program does, is comparing the command
line arguments of the program (commonly called
{\tt string argv[]} in a C program) 
with this default vector of ``{\it keyword=value}''
strings ({\tt string defv[]}), and replace
appropriate reset values for later retrieval. This is done by calling
{\tt initparam}\footnote{It secretly assumes that {\tt argv[]} is
NULL terminated, which is not guaranteed on all UNIX 
implementations}\index{initparam} as the first step in your MAIN
program:

\begin{verbatim}
    main (argc, argv)
    int argc;
    string argv[];
    {
        initparam(argv,defv);
        . . .
\end{verbatim}

It also checks if keywords which do not have a default
value ({\it i.e.} were given ``{\tt ???}'') 
have really been given a proper 
value on the command line, if keywords are not specified twice, 
enters values of the system keywords etc.

There is a better alternative to define the main part of a NEMO program:
by renaming the main entry point {\tt main()} to 
{\tt nemo\_main()} \index{nemo\_main}, without any arguments, and
calling the array of strings with default '{\it key=val}'s
{\tt string defv[]}, the linker will automatically include
the proper startup code ({\tt initparam(argv,defv)}), the worker routine {\tt
nemo\_main()} itself, and the stop code ({\tt finiparam()}). 
The above section of code would then be replaced by a mere:

\begin{verbatim}
    nemo_main()
    {
        . . .
\end{verbatim}

This has the obvious advantage that various NEMO related 
administrative details
are now hidden from the application programmers, and occur automatically.
Remember that standard {\tt main()} already shields the application
programmer from a number of tedious setups (e.g. {\it stdio} etc.).
Within NEMO we have taken this one step further.
The example given later in Section~\ref{ss-example} also uses the
technique of calling the main entry point {\tt nemo\_main()}.

Once the user interface has been initialized, keyword values may be obtained at
any point during execution of the program by calling {\tt getparam()},
which returns a string\footnote{note that ANSI rules say
you can't write to this location in memory if they are direct references
to {\tt string defv[]}; this is something that may well be fixed in a future
release}:

\begin{verbatim}
    if (streq(getparam("n"),"0")
        printf(" You really mean zero or octal?\n");
\end{verbatim}

There is a whole family of {\tt getXparam()} \index{getXparam, parsing}
\index{parsing, getXparam} functions which 
parse\footnote{Depending on compiler switches at installation
the getXparam parsing includes full expressions} the
string in a value of one of the basic C types {\tt int, long, bool,}
and {\tt real}. It returns that value in that type: 

\begin{verbatim}
    #include <getparam.h>     /* defines ``int getiparam()'' */
    . . .
    int nbody;
    . . .
    if ( (nbody = getiparam("n")) <= 0) {
        printf("Cannot handle %d particles\n",nbody);
        exit(0);
    }
\end{verbatim}

Finally, there is a macro called {\tt getargv0()}, which returns
\index{getargv0, program name} \index{program name, getargv0}
the name of the calling program, mostly used for identification:

\begin{verbatim}
    if (getbparam("quit"))
        error("%s: early quit", getargv0());
\end{verbatim}

This is very useful in library routines, who normally would not be
able to know who called them. Actually, NEMO's {\tt error} function
already uses this technique, since it cannot know the name of the
program by whom it was called.
The {\tt error} function prints a message, and exits the 
program\index{error(3)}.
% unless {\tt recover} was called to catch error messages.

More detailed information can also be found in the appropriate manual
page: {\it getparam(3NEMO)} and {\it error(3NEMO)}.

\subsection{Advanced User Interface and String Parsing}
\mylabel{ss:advusr}

Here we describe {\it setparam} to add some interactive capabilities in
\index{setparam} a standard way to NEMO.  Values of keywords should only
be accessed and modified \index{keywords, modifying} this way.  Since
keywords are initialized/stored within the source code, most compilers
will store their values in a read-only part of data area in the
executable image.  Editing them may cause unpredictable behavior. 

If a keyword string contains an array of items of the same type, 
one can use either {\tt nemoinpX} or 
{\tt getXrange}, \index{parsing, nemoinpX}\index{parsing, getXrange}
\index{nemoinpX, parsing}\index{getXrange, parsing}
depending if you know how many items to expect in the string.
The {\tt getXrange} routines will allocate a new array which will contain
the items of the parsed string. If you do already have a
declared array, and know that all items will fit in there, the
{\tt nemoinpX} routines will suffice.

An example of usage:

\begin{verbatim}
    double *x = NULL;
    double y[NYMAX];
    int nxret, nyret;
    int nxmax=0;

    nyret = nemoinpd(getparam("y"), y, NYMAX);

    nxret = getdrange(getparam("x"), &x, &nxmax);
\end{verbatim}

In the first call the number of elements to be parsed
from an input keyword {\tt y=} is limited to {\tt NYMAX}, and is useful
when the number of elements is expected to be small or more
or less known. The actual
number of elements returned in the array {\tt y[]} is {\tt nyret}.

When the number of elements to be parsed is not known at all, or one needs
complete freedom, the dynamic allocation feature of {\tt getdrange} can
be used. The pointer {\tt x} is initialized to {\tt NULL}, as well as
the item counter {\tt nxmax}.
After
calling {\tt getdrange}, {\tt x} will point to an array of length
{\tt nxmax}, in which the first {\tt nxret} element contain the parsed
values of the input keyword {\tt x=}. Proper re-allocation will be done
when a larger space is need on subsequent calls.

Both routines return negative error return codes, see {\it nemoinp(3NEMO)}.

More complex parsing is also done by calling {\tt burststring} first
to break a string in pieces, followed by a variety of functions.

\subsection{Alternatives to nemo\_main}

It is not required for your program to define with {\tt nemo\_main()}. There
are cases where the user needs more control.  An example of this is the
{\tt falcON}\index{falcON} N-body code in {\tt \$NEMO/use/dehnen/falcON}.
A header file (see e.g. {\tt inc/main.h}) now defines main, instead of
through the NEMO library:

\begin{verbatim}
// in main.h:

extern string defv[];
extern string usage;

namespace nbdy {
  char version [200];                             // to hold version info
  extern void main();                             // to be defined by user
};

int main(int argc, char *argv[])                  // ::main()
{
  snprintf(nbdy::version,200,"VERSION=" /*..*/ ); // write version info
  initparam(argv,defv);                           // start NEMO
  nbdy::main();                                   // call nbdy::main()
  finiparam();                                    // finish NEMO
}
\end{verbatim}

and the application includes this header file, and defines the keyword list
in the usual way :

\begin{verbatim}
// in application.cc

#include <main.h>

string defv[] = { /*...*/, nbdy::version, NULL }; // use version info
string usage  = /*...*/ ;

void nbdy::main() { /*...*/ }                     // nbdy::main()

\end{verbatim}


\subsection{filestruct.h}

\mylabel{ss-filestr}
The {\it filestruct} package provides a direct and consistent way of passing
data between NEMO programs, much as {\it getparam} provides a way of passing
(command line) arguments to programs.  \index{filestruct.h}
For reasons of economy and
accuracy, much of the data manipulated by NEMO is stored on disk in
binary form.  Normally, data stored this way is completely
unintelligible, except to specialized programs which create and access
it.  Furthermore, this approach to data handling tends to be very
brittle: a trivial addition or alteration to the data stored in such a
file can force the tedious and error-prone revision of many programs. 
To get around these problems and provide an explicit, flexible, and
structured method of storing binary data,  we developed a collection of
general purpose routines to access binary data files. 

From the programmers point of view, a structured binary file 
\index{structured file} is a stream
of tagged data objects.  These objects come in two classes.  An {\it
item}\index{item}
is a single instance or a regular array of one of the following C
primitive types: {\tt char, short, int, long, float} or {\tt double}.  A
{\it set} \index{set}
is an unordered sequence of {\it items} and {\it sets}. This definition is
recursive, so fully hierarchical file structures are allowed, and indeed
encouraged.  Every set or item has a name {\it tag} \index{tag, item}
associated with it,
used to label the contents of a file and to retrieve objects from a set. 
Data items have a {\it type} \index{type, item}
and array {\it dimension} \index{dimension, item} attributed
associated with them as well. This of course means that there is a little
overhead, which may become too large if many small amounts of data
are to be handled. For example, a snapshot  with 128 bodies
(created by {\tt mkplummer}\index{mkplummer}) with double
precision masses and full 6 dimensional phase space coordinates
totals 7425 bytes, whereas a straight dump of only the essential
information would be 7168 bytes, a mere 3.5\% overhead. After an
integration, with 9 full snapshots stored and 65 snapshots with only
diagnostics output, the overhead is much larger: 98944 bytes of
data, of which only 64512 bytes are masses and phase space coordinates:
the overhead is 53\% (of which 29\% though are the diagnostics output,
such conservation of energy and angular momentum, cputime, center
of mass, etc.).


The filestruct package uses ordinary {\it stdio(3)} streams to access
input and output files; \index{stream}
hence the first step in using filestruct is to
open the file streams.  For this job we use the NEMO library routine
{\tt stropen()}, \index{stropen}
which itself is not part of filestruct.  {\tt
stropen({\it name,mode})} is much like {\tt fopen()} of {\it stdio}, but
slightly more clever; it will not open an existing file for output,
unless the {\it mode} string is {\tt "w!"}. 
% append ???
An additional oddity to
{\tt stropen} is that it treats the dash filename {\tt "-"},
\index{-,filename} as standard 
in/output,\footnote{Older versions of 
{\it filestruct} cannot handle binary files in pipes, since
filestruct uses fseek(3)}
and {\tt "s"} as a scratch file.
Since {\it stdio} normally flushes all buffers on exit, it is
often not necessary to explicitly close open streams, but if you do so,
use the matching routine {\tt strclose()}. This\index{strclose}
also frees up the table
entries on temporary memory used by the filestruct package. As in most
applications/operating systems a task can have a limited set of
open files associated with it. Scratch\index{scratch files} files 
are automatically deleted from disk
when they are closed. \index{scratch files}

Having opened the required streams, it is quite simple to use the basic
data I/O routines.  For example, suppose the following declarations have
been made:

\begin{verbatim}
    #include <stdinc.h>
    #include <filestruct.h>

    stream instr, outstr;
    int    nbody;
    string headline;

    #define MAXNBODY 100
    real    mass[MAXNBODY];
\end{verbatim}

(note the use of the {\tt stdinc.h} conventions).  And now suppose that, 
after some computation, results have been stored in the first {\tt nbody}
components of the {\tt mass} array, and a descriptive message has been
placed in {\tt headline}.  The following piece of code will write the
data to a structured file:

\begin{verbatim}
    outstr = stropen("mass.dat", "w");

    put_data(outstr, "Nbody", IntType, &nbody, 0);
    put_data(outstr, "Mass", RealType, mass, nbody, 0);
    put_string(outstr, "Headline", headline);

    strclose(outstr);    
\end{verbatim}

\index{put\_data()} \index{put\_string}
Data (the 4th argument in {\tt put\_data},
is always passed by address, even if one element is written. This
not only holds for reading, but also for writing, as is apparent from
the above example.
Note that no error checking is needed when the file is opened for
writing. If the file {\tt mass.dat}
would already have existed, {\tt error()} would \index{error(3)}
have been called inside {\tt stropen()} and aborted the
program. Names of tags are arbitrary, but we encourage you to use
descriptive names, although an arbitrary maximum of 64 is enforced
by chopping any incoming string.

The resulting contents of {\tt mass.dat} can be viewed with the
{\tt tsf} utility: \index{tsf}

\begin{verbatim}
    % tsf mass.dat
    int Nbody 010
    double Mass[8] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    char Headline[20] "All masses are equal"
\end{verbatim}

Note the octal representation {\tt 010=8} of {\tt Nbody}. **OLD**

% CANNOT DO THIS YET - CONCEPT NOT INTRODUCED HERE
% Also note
% that it is not {\bf required} to use sets, subsets, subsubsets etc.,
% but such encapsulation has certain advantages.

It is now trivial to read data from this file:

\begin{verbatim}
    instr = stropen("mass.dat", "r");

    get_data(instr, "Nbody", IntType, &nbody, 0);
    get_data(instr, "Mass", RealType, mass, nbody, 0);
    headline = get_string(instr, "Headline");

    strclose(instr);
\end{verbatim}

Note that we read the data in the same order as they were written.

\index{get\_data} \index{get\_string}
During input, the filestruct routines normally perform strict type-checking;
the tag, type and dimension supplied to {\tt get\_data()} must match
the attributes of the data item, written previously,
exactly. Such strict checking helps
prevent many common errors in using binary data. Alternatively, you can use
{\tt get\_data\_coerced()}, which is called exactly like {\tt get\_data()}, but
interconverts {\tt float} and {\tt double} 
values\footnote{The implementors of NEMO will not be held responsible for
any loss of precision resulting from the use of this feature.}.\index{real}

To provide more flexibility in programming I/O, a series of related items
may be hierarchically wrapped into a set:

\begin{verbatim}
    outstr = stropen("mass.dat", "w");

    put_set(outstr, "NotASnapShot");
       put_data(outstr, "Nbody", IntType, &nbody, 0);
       put_data(outstr, "Mass", RealType, mass, nbody, 0);
       put_string(outstr, "Headline", headline);
    put_tes(outstr, "NotASnapShot");

    strclose(outstr);    
\end{verbatim}

Note that each {\tt put\_set()}\index{put\_set}
must be matched by an equivalent  {\tt put\_tes()}\index{put\_tes}.
For input, corresponding routines {\tt get\_set()} and {\tt get\_tes()}
are used. These also introduce a significant additional functionality:
between a {\tt get\_set()}\index{get\_set}
and {\tt get\_tes()}\index{get\_tes}, the data items
of the set may 
be read in any order\footnote{tagnames must now be unique within an item,
as in a C {\tt struct}},
or not even read at all. For 
example, the following is also a legal way to access 
the {\tt NotASnapShot}\footnote{This is a toy model, shown for its
simplicity. The full {\tt SnapShot} format is discussed in 
Section \ref{ss:snapshot}}:

\begin{verbatim}
    instr = stropen("mass.dat", "r");

    if (!get_tag_ok(instr,"NotASnapShot"))
        error("File mass.dat is not a NotASnapShot\n");

    get_set(instr,"NotASnapShot");
    headline = get_string(instr, "Headline");
    get_tes(instr,"NotASnapShot");

    strclose(instr);    
\end{verbatim}

\index{get\_tag\_ok}
This method of ``filtering'' a data input stream clearly opens up many
ways of developing general-purpose programs. Also note that the 
{\tt bool} routine {\tt get\_tag\_ok()} can be used to control the
flow of the program, as {\tt get\_set()} would call {\tt error()} when
the wrong tag-name would be encountered, and abort the program.

The UNIX program {\tt cat} can also be used to catenate\index{cat}
multiple binary data-sets into one, 
\index{data, catenation}
{\it i.e.}

\begin{verbatim}
    %   cat mass1.dat mass2.dat mass3.dat > mass.dat
\end{verbatim}

The {\tt get\_tag\_ok} routine can be used to handle such multi-set data
files.  The following example shows how loop through such a combined
data-file. 

\begin{verbatim}
    instr = stropen("mass.dat", "r");

    while (get_tag_ok(instr, "NotASnapShot") {
       get_set(instr, "NotASnapShot");
       get_data(instr, "Nbody", IntType, &nbody, 0);
       if (nbody > MAXNBODY) {
            warning("Skipping data with too many (%d) items",nbody);
            get_tes(instr,"NotASnapShot");
            continue;
       }
       get_data(instr, "Mass", RealType, mass, nbody, 0);
       headline = get_string(instr, "Headline");
       get_tes(instr,"NotASnapShot");
       /*   process data   */
    }

    strclose(instr);    
\end{verbatim}

The loop is terminated at either end-of-file, or if the next object in
{\tt instr} is not a {\tt NotASnapShot}. 

It is easy to the skip for an item if you know if it is there:


\begin{verbatim}
    while(get_tag_ok(instr,"NotASnapShot"))     /*  ??????? */
        skip_item(instr,"NotASnapShot");
\end{verbatim}

The routine {\tt skip\_item()} is only effective, or for that matter
required, when doing input at the top level, {\it i.e.} not between
a {\tt get\_set()} and matching {\tt get\_tes()}, since I/O at deeper
levels is random w.r.t. items and sets. In other words,
at the top level I/O is sequential, at lower levels random.

A relative new feature in data access is the ability to
do random and blocked access to the data. Instead of using
a single call to {\tt get\_data} and {\tt put\_data}, the
access can be sequentially blocked using 
{\tt get\_data\_blocked} and {\tt put\_data\_blocked}, provided
it is wrapped using {\tt get\_data\_set} and {\tt get\_data\_tes},
for example:

\begin{verbatim}
   get_data_set    (instr, "Mass", RealType, nbody, 0);
   real *mass = (real *) allocate((nbody/2)*sizeof(real));
   get_data_blocked(instr, "Mass", mass, nbody/2);
   get_data_blocked(instr, "Mass", mass, nbody/2);
   get_data_tes    (instr, "Mass");
\end{verbatim}

would read in the {\tt Mass} data in two pieces into a smaller sized {\tt mass}
array. A similar mode exists to randomly access data with an item. A current
limitation of this mode is that such access is only allowed on one item at a
time. In this mode an item must be closed before the next one can be opened
in such a mode.\index{blocked I/O}\index{filestruct,blocked I/O}
\index{filestruct,random I/O}\index{random access}

\subsection{vectmath.h}

The {\tt vectmath.h} macro package \index{vectmath.h}
provides a set of macros to handle
some elementary operations on two, three or general $N$ dimensional
vectors and matrices.  The dimension $N$ can be picked by providing the
package with a value for the preprocessor macro {\bf NDIM}.  If this is
not supplied, the presence of macros {\bf TWODIM} and {\bf THREEDIM} will
be checked, in which case {\bf NDIM} is set to 2 or 3 respectively.
The default of {\bf NDIM}
when all of the above are absent,
is 3. Of course, the macro {\bf NDIM} must be provided before
{\tt vectmath.h} is included to have any effect.
Resetting the value of {\bf NDIM} after that,
if your compiler would allow it anyhow without an explicit {\tt \#undef},
may produce unpredictable results.

There are also a few of the macro's which can be used as a regular
C function, returning a real value, {\it e.g.} {\tt absv()} for the
length of a vector.

Operations such as {\bf SETV} (copying a vector) are properly defined for
every dimension, but {\bf CROSSVP} (a vector cross product) has a
different meaning in 2 and 3 dimensions, and is
absent in higher dimensions.

It should be noted that the matrices used here are true C
matrices, a pointer to an array of pointers (to 1D arrays), unlike
FORTRAN arrays, which only occupy a solid 2D block of memory. C
arrays take slightly more memory. For an example how to make C arrays and
FORTRAN arrays work closely together see e.g. {\it Numerical
Recipes in C} by {\it Press et al.} (MIT Press, 1988).
\index{Numerical Recipes} \index{Press W.}

In the following example a 4 dimensional vector is cleared:

\begin{verbatim}
    #define NDIM 4
    #include <vectmath.h>

    nemo_main()
    {
        vector a;           /* same as:   double a[4]  */

        CLRV(a);
    }
\end{verbatim}

{\it some more examples here - taken from some snap code}

\subsection{snapshots: get\_snap.c and put\_snap.c}
\mylabel{ss:snapshot}

These routines exemplify an attempt to provide truly generic I/O of
N-body data.  They read and write structured binary data files
conforming to the overall form seen in earlier sections.
Internally they
operate on {\tt Body} structures; \index{Body, structure}
A {\tt Body} has components accessed
by macros such as {\tt Mass} for the mass, {\tt Pos} and {\tt Vel} for
the position and velocity vectors, {\it etc.}.  Since {\tt get\_snap.c}
and {\tt put\_snap.c} use only these macros to declare and access
bodies, they can be used with any suitable {\tt Body} structure.  They
are thus provided as C source code to be included in the compilation of
a program. Definitions concerning Body's and snapshots are obtained 
by including the files {\tt snapshot/body.h} and {\tt snaphot/snapshot.h}.
\index{snapshot.h}\index{body.h}

A program which should handle a large number of particles, may decide to
include a more simple {\tt Body} structure, as is {\it e.g.} provided by
the {\tt snapshot/barebody.h} macro file.  This body only includes the
masses and phase space coordinates, which would only occupy 28 bytes 
per particle (in single precision), as opposed to the 100/104
bytes per particle
for a double precision {\tt Body} from the standard {\tt snapshot/body.h} 
macro file.  This last one contains {\tt Mass, PhaseSpace, Phi, Acc, Aux}
and {\tt Key}.

In the example listed under Table~\ref{src:snapfirst}
the first snapshot of an input file is
copied to an output file. 

\begin{table}[tb]
\caption[source: snapfirst.c]{\$NEMO/src/tutor/snap/snapfirst.c}
\mysrcfile{snapfirst.src}
\mylabel{src:snapfirst}
\small
\verbatimlisting{snapfirst.src}
\normalsize
\end{table}

Notice that the first argument of {\tt stropen()}, the filename, is
directly obtained from the user interface. The input file is
opened for reading, and the output file for writing.
Some history\footnote{See next section for more details on
history processing}
is obtained from the input file (would we not have done
this, and the input file would have contained history, 
a subsequent {\tt get\_snap()} call would have failed to find
the snapshot), and the first snapshot is read into an array
of bodies, pointed to by {\tt btab}. Then the output file
has the old history written to it (although any command line
arguments were added to that), followed by that first snapshot.
Both files are formally closed before the program then returns.

\subsection{history.h}

When performing high-level data I/O, as is offered by a package such as
{\tt get\_snap.c} and {\tt put\_snap.c}, there is an automated way to
keep track of data history. \index{history.h}

When a NEMO program is invoked, the program name and command line
arguments are saved by the {\tt initparam()} in a special history
database.  Most NEMO programs will write such history items to their
data-file(s) before the actual data.  Whenever a data-file is then
opened for reading, the programmer should first read these data-history
items. Conversely, when writing data, the history should
be written first. In case of the {\tt get/put\_snap} package:

\begin{verbatim}
    get_history(instr);
    get_snap(instr, &btab, &nbody, &time, &bits);
       /*     process data     */
    put_history(outstr);
    put_snap(outstr,&btab, &nbody, &time, &bits);
\end{verbatim}

\index{get\_history} \index{put\_history}
Private comments should be added with the {\tt app\_history()}
\index{app\_history}\index{add\_history}\index{readline, GNU}
\footnote{The old name, {\tt add\_history} was already used by the
GNU {\it readline} library}
When a series of snapshot is to be processed, it is recommended that
the program should only be output
the history once, before the first output of the snapshot, as in the
following example:

\begin{verbatim}
    get_history(instr);
    put_history(outstr);
    for(;;) {
        get_history(instr);     /* defensive but in-active */
        get_snap(instr, &btab, &nbody, &time, &bits);
            /*     process data  and decide when done */
        put_snap(outstr,&btab, &nbody, &time, &bits);
    }
\end{verbatim}

Note that the second call to {\tt get\_history()}, within the for-loop, 
is really in-active. If there happen to be history items sandwiched
between snapshots, they will be read and added to the history stack,
but not written to the output file, since {\tt put\_history()} was only 
called before the for-loop. It is only a defensive call:
{\tt get\_snap()} would fail since it expects only pure
{\tt SnapShot} sets (in effect, it calls {\tt get\_set(instr,"SnapShot")}
first, and would call {\tt error()} if no snapshot encountered).

\section{Building NEMO programs}
\mylabel{s:example}

Besides writing the actual code for a program, an application programmer
has to take care of a few more items before the software can be added 
and formally be accepted to NEMO.  This concerns 
writing\index{compiling, NEMO programs}\index{linking, NEMO programs}
the documentation and possibly a Makefile, the former one preferably 
in the form of standard UNIX manual pages ({\it man(5)}).
We have templates for both Makefile's and 
manual pages. Both these are discussed in detail in the next subsections. 

Because NEMO is a development package within which
a multitude of people are donating software and
libraries, linking a program can become cumbersome. In the most
simple case however (no graphics or mathematical libraries needed),
only the main NEMO library is needed, and the
following command should suffice to produce an executable:

\begin{verbatim}
    % cc -g -o snapprint snapprint $NEMOLIB/libnemo.a -lm
or:
    % cc -g -o snapprint snapprint -lnemo -lm
\end{verbatim}

The second form would only work if your {\it cc} compiler understands
the {\bf -L} switch, and the {\tt \$NEMOBIN/cc} has installed this
feature.  See Appendix~\ref{a:install} how to properly install this
script. 

For graphics programs a solution would be to use the {\bf YAPPLIB} 
environment variable. 
\index{YAPPLIB, environment}

An example of the compilation of a graphics program:
\begin{verbatim}
    % cc -g -o snapplot snapplot.c -lnemo $YAPPLIB -lm
\end{verbatim}

Each user is given a subdirectory in {\tt \$NEMO/usr}, under which
code may be donated which can be compiled into the running version
of NEMO. Stable code, which has been sufficiently tested and
verified, can be placed in one of the appropriate {\tt \$NEMO/src} 
directories. For proper inclusion of user contributed software
a few rules in the {\tt Makefile} have to be adhered to.

The {\tt bake}\index{bake, make script}
and {\tt mknemo}\index{mknemo, script}
script should handle compilation and
installation of most of the standard NEMO cases. Some programs, like 
the N-body integrators, are almost like complicated packages themselves,
and require their own Makefile or install script. For most
programs you can compile it by:

\begin{verbatim}
    % bake snapprint
\end{verbatim}

or to install:\footnote{this assumes you have some 
appropriate NEMO permissions}

\begin{verbatim}
    % mknemo snapprint
\end{verbatim}

\subsection{Manual pages}

\index{manual, online}
It is very important to keep a manual file (preferably in the UNIX man format)
online for every program. A program that does not have an accompanying manual
page is not complete. Of course there is always the inline help
({\tt help=}) that every NEMO program has.

To a lesser degree this also applies to the public
libraries. A template roff sample can be found in {\tt example.8}.
We encourage authors to have a MINIMUM set of sections in a man-page as
listed below. The ones with a '*' are considered somewhat less 
important:

\begin{list}{-}{\labelwidth 75pt \leftmargin 75pt}
\item[{\bf NAME}]
                  the name of the beast.
\item[{\bf SYNOPSIS}]
                  command line format or function prototype, include
                    files needed etc.
\item[{\bf DESCRIPTION}]
                  maybe a few lines of what it does, or not does.
\item[{\bf PARAMETERS}]
                  description of parameters, their meaning and default values.
        This usually applies to programs only. 
\item[{\bf EXAMPLES}]
                (*) in case non-trivial, but recommended anyhow
\item[{\bf DEBUG}]
                (*) at what {\tt debug} levels what output appears.
\item[{\bf SEE ALSO}]
                (*) references to similar functions, more info
\item[{\bf BUGS}]
                (*) one prefers not to have this of course
\item[{\bf TIMING}]
                (*) performance, dependence on parameters if non-trivial
\item[{\bf STORAGE}]
                (*) storage requirements - mostly of importance when programs
                  allocate memory dynamically, or when applicable for the
                  programmer.
\item[{\bf LIMITATIONS}]
                (*) does it have any obvious limitations?
\item[{\bf AUTHOR}]
                  who wrote it (a little credit is in its place) and/or who
                is responsible.
\item[{\bf FILES}]
           (*) in case non-trivial
\item[{\bf HISTORY}]
    date, version numbers, why updated, by whom (when created)
\end{list}


\subsection{Makefiles} 
\mylabel{ss:makefile}

Makefiles are scripts in which "the rules are defined to make targets", see
{\it make(1)} for many more details. In other words, the Makefile tells
\index{makefile, NEMO rules}
how to compile and link libraries and programs. NEMO uses Makefiles
extensively for installation, updates and various other system
utilities. Sometimes scripts are also available
to perform tasks that can be done by a Makefile.

There are basically three types of Makefiles in NEMO:

\begin{description}
\item{1.} The first (top) level Makefile. It lives in NEMO's root
directory (normally referred to as {\tt \$NEMO})
and can steer installation on any of a number of selected machines,
it includes some import and export facilities (tar/shar)
and various other system maintenance utilities.
At installation it makes sure all directories are present,
does some other initialization, and then calls Makefile's 
one level down to do the rest of the dirty work. The top level
Makefile is not of direct concern to an application programmer,
nor should it be modified without consent of the NEMO
system manager.

\item{2.} Second level Makefiles, currently in {\tt \$NEMO/src}
and {\tt \$NEMO/usr},
steer the building of libraries and programs by calling Makefiles in
subdirectories one more level down.  Both this 2nd level Makefile and
the one described earlier are solely the responsibility of NEMO system
manager.  You don't have to be concerned with them, except to know of 
their existence because your top level Makefile(s) must be callable
by one of the second level Makefiles. This interface will be described 
next.

\item{3.} Third level Makefiles live in source or user directories
{\tt \$NEMO/src/topic} and {\tt \$NEMO/usr/name} (and possibly below). 
They steer the installation of user specific 
programs and libraries, they may update NEMO libraries too.
The user writes his own Makefile, he usually splits up his directory 
in one or more subdirectories, where the real work is done by what
we could then call level 4 or even level 5 Makefiles. However, 
this is completely the freedom of a user.
The level 3 Makefiles normally have two kinds of entry points
(or 'targets'): the user 'install' targets are used by the user, 
and make sure this his sources, binaries, libraries, include files etc. are
copied to the proper places under {\tt \$NEMO}.
The second kind of entry point
are the 'nemo' targets and never called by you, the user;
they are only called by Makefiles
one directory level up from within {\tt \$NEMO} below during the
rebuilding process of NEMO, {\it i.e.} a user never calls 
a nemo target, NEMO will do this during its installation.
Currently we have NEMO install itself in two phases, resulting in
two 'nemo' targets: 'nemo\_lib' (phase 1) and 'nemo\_bin' (phase 2). 
A third 'nemo' target must be present to create a lookup
table of directories and targets for system maintenance. 
This target must be called 'nemo\_src', and must also call lower
level Makefiles if applicable.

This means
that user Makefiles {\bf MUST} have at least these three targets in order
to rebuild itself from scratch. In case a user
decides to split up his directories, the Makefiles must also
visit each of those directories and make calls through the same
entry points 'nemo\_lib' and 'nemo\_bin', 'nemo\_src'; 
a sort of hierarchical install process.
   
\end{description}

For more details see the template Makefiles in NEMO's sec subdirectories
and the example below in section~\ref{ss-example}.

We expect a more general install mechanism with a few more strict rules for
writing Makefiles, in some next release of NEMO. 


\newpage
\subsection{An example NEMO program}

\mylabel{ss-example}
Under Table~\ref{src:hello}
below you can find a listing of a very minimal NEMO program, 
``{\tt hello.c}'':\index{nemo\_main}\index{hello.c}

\begin{table}[tb]
\caption[source: hello.c]{\$NEMO/src/tutor/hello/hello.c}
\mylabel{src:hello}
\mysrcfile{hello.src}
\small
\verbatimlisting{hello.src}
\normalsize
\end{table}

and a corresponding example Makefile to install by {\em user} and 
{\em nemo} could look like the one shown under Table~\ref{src:makefile}

Note that for this simple 
example the {\tt Makefile} actually larger than
the source code, {\tt hello.c}, itself. Fortunately 
not every programs needs
their own Makefile,  in fact most programs can be compiled with
a default rule, via the {\tt bake}\index{bake, make script}
script. This generic makefile is used by the {\tt bake} command, 
and is normally installed in {\tt \$NEMOLIB/Makefile}, but check
out your {\tt bake} command or alias.

% \newpage
% \centerline{\bf Example Makefile}

\begin{table}[htb]
\caption[source: makefile]{Sample makefile - cf. \$NEMOLIB/Makefile}
\mylabel{src:makefile}
\mysrcfile{makefile.src}
\footnotesize
\verbatimlisting{makefile.src}
\normalsize
\end{table}

{\bf Warning:}\ The structure of this
so-called 'standard' NEMO Makefiles is still under 
debate, and will probably drastically 
change in some future release. Best is to check some
local Makefiles. A possible candidate is the GNU\index{make, GNU}
make\index{gmake, GNU} facility.

\section{Extending NEMO environment}

Let us now summarize the steps to follow to add and/or create new software to
NEMO.  The examples below are suggested steps taken from adding
Aarseth's {\tt nbody0} program to NEMO, and we assume him to have his
original stuff in a directory {\tt \~/nbody0}. 

\begin{description}

\item[1:] Create a new directory, {\tt "cd \$NEMO/usr ; mkdir aarseth"}
        and inform the system manager of NEMO that a new user should be
        added to the user list in {\tt \$NEMO/usr/Makefile}. You can also
        do it yourself if the file is writable by you.

\item[2:] Create working subdirectories in your new user directory,
        {\tt "cd aarseth ; mkdir nbody0"}.

\item[3:] Copy a third level Makefile from someone else, and substitute
        the subdirectory names to be installed for you, i.e. your new
        working subdirectories ({\tt 'nbody0'} in this case):
        {\tt "cp ../pjt/Makefile . ; emacs Makefile"}.

\item[4:] Go 'home' and install, {\tt "cd \~/nbody0 ; make install"}, assuming
        the Makefile there has the proper install targets. Check the
        target Makefile in the directory 
        {\tt \$NEMO/usr/aarseth/nbody0} what this last command
        must have done.
        
\end{description}

Actually, only step 1 is required. If a user cannot or does not want 
to confirm to the level 3/4 separation, he may do so, as long as the 
Makefile in level 3 (e.g. {\tt \$NEMO/usr/aarseth/Makefile}) contains the 
nemo\_lib, nemo\_bin and nemo\_src install targets. An example of adding
a foreign package that way is the {\tt GRAVSIM} package \index{GRAVSIM},
which has it's own internal structure. In the directory tree starting
at  {\tt \$NEMO/usr/mbellon/gravsim} an example of a different approach
is given. Sometimes public domain packages have been added to NEMO, and
its Makefiles have been adapted slightly to the NEMO install procedure.

\section{Programming in C++}
\mylabel{s:cplusplus}

Most relevant header files from the NEMO C libraries have been made
entrant for C++. This means that all routines should be available 
through:

\begin{verbatim}
    extern "C"  {    
        .....    
    }
\end{verbatim}

The only requirement is of course that the {\tt main()} be in C++.
For this you have to link with the NEMO++ library {\bf before}
the regular NEMO library. So, assuming your header (-I) 
and library (-L) include flags have been setup, you should be able to
compile your C++ programs as follows:

\begin{verbatim}
    % cppc -g -o test test.cc -lnemo++ -lnemo -lm
\end{verbatim}


\section{Programming in FORTRAN}
\mylabel{s:fortran}

Programming in FORTRAN \index{FORTRAN, programming in} can also be done,
but since NEMO is written in C and there is no 'standard' way to link
FORTRAN and C code, such a description is always bound to be system
dependent (large differences exist between UNIX, VMS, MSDOS, and
UNICOS is somewhat of a peculiar case).
Even within a UNIX environment there are a number of ways how the
industry has solved this problem (cf. Alliant). Most comments that
will follow, apply to the BSD convention of binding FORTRAN and C.

In whatever language you program, 
we do suggest that the startup of the program is done in C,
preferably through \index{nemo\_main} the {\tt nemo\_main()} function (see
Section~\ref{ss-example}).  
As long as file I/O is
avoided in the FORTRAN routines, character and boolean 
variables are avoided in arguments of C callable FORTRAN functions, 
all is relatively simple. Some care is also needed for multidimensional
arrays which are not fully utilized.
The only thing needed are C names of the FORTRAN routines to be called
from C.  This can be handled automatically by a macro package. 

Current examples can be found in the programs
{\tt nbody0} and {\tt nbody2}. In both cases data
file I/O is done in C in NEMO's  {\it snapshot(5NEMO)} format,
but the CPU is used in the FORTRAN code. 

Examples of proposals for other FORTRAN interfaces can be found
in the directory {\tt \$NEMOINC/fortran}.

% In case a truly full fledged FORTRAN-c interface is needed, look in
% the\index{GIPSY} directory {\tt \$NEMO/usr/gipsy/util/f2c}

Again this remark: the {\it potential(5NEMO)} assumes for now a BSD type
f2c interface, because character variables are passed.  This has not
been updated yet.  You would have to provide your own f2c interface to
use FORTRAN potential routines on other systems. 

Simple FORTRAN interface workers within the snapshot interface are
available in a routine {\tt snapwork(n,m,pos,vel,...)}.


{\it More description will follow}

\subsection{Calling NEMO C routines from FORTRAN}
\mylabel{ss:cfromf}

The NEMO user interface, with limited capabilities, is 
also available to FORTRAN\index{fortran, calling C}
programmers. First of all,
the keywords, their defaults and a help string must be made
available (see Section~\ref{ss:getparam}).
This can be done by supplying them as comments
in the FORTRAN source code, as is show in the following 
example listed under Table~\ref{src:testf2c}

\begin{table}[htb]
\caption[testf2c.f]{\$NEMO/src/kernel/fortran/test.f}
\mylabel{src:testf2c}
\mysrcfile{testf2c.src}
\small
\verbatimlisting{testf2c.src}
\normalsize
\end{table}

The documentation section between {\tt C+} and {\tt C-} can be 
extracted\index{ftoc}
with a NEMO utility, {\tt ftoc}, to the appropriate C module as follows:

\begin{verbatim}
    % ftoc test.f test_main.c
\end{verbatim}

after which the new {\tt test\_main.c} file merely has to be included on
the commandline during compilation. To avoid having to include
FORTRAN libraries explicitly on the commandline, easiest is
to use the {\tt f77} command, instead of {\tt cc}:

\begin{verbatim}
    % f77 -o test test.f test_main.c -I$NEMOINC -L$NEMOLIB -lnemo
\end{verbatim}

This only works if your operating supports mixing
C and FORTRAN source code on one commandline. Otherwise try:

\begin{verbatim}
    % cc -c test_main.c
    % f77 -o test test.f test_main.o -L$NEMOLIB -lnemo
\end{verbatim}
where the NEMO library is still needed to resolve the
user interface of course.

The other alternative would be:
\begin{verbatim}
    % f77 -c test.f
    % cc -o test test.o test_main.c -L$NEMOLIB -lnemo \
                -lF77 -lI77 -lU77 -lm
\end{verbatim}
with various possible complications with the new FORTRAN 1.3+ compiler
on SUN workstations. Browsing with  {\it nm(1)} through UNIX library
files to find undefined reference might be the only alternative
left to find out where the system has hidden them. See the 
{\tt FORLIBS}\index{FORLIBS, environment}
environment variable defined in {\tt NEMORC} startup file. 



\subsection{Calling FORTRAN routines from NEMO C}
\mylabel{ss:ffromc}

No official support 
is needed, although for portablility it would
be nice to include a header file that maps the symbol names
and such. 


\section{Debugging}
\mylabel{s:debugging}

Apart from the usual debugging methods that everybody knows about,  NEMO
programs usually have the following additional properties which can cut
down in debugging time. If not conclusive during runtime, you can either
decide to compile the program with debugging flags turned on, and run
the program through the debugger, or add more {\tt
dprintf}\index{dprintf(3)} or {\tt error}\index{error(3)} function
calls:

\begin{itemize}
\item
During runtime you can set the value for the {\tt debug=}\index{debug,
system keyword} (or use the equivalent {\tt DEBUG} environment variable)
system keyword to increase the amount of output. Note that only levels
0 (the default) through 9 are supported. 9 should produce a lot of output.

\item
During runtime you can set the value for the {\tt error=}\index{error,
system keyword} (or use the equivalent {\tt ERROR} environment variable)
system keyword to bypass a number of fatal error messages that you
know are not important. For example, to overwrite an existing file
you would need to increase {\tt error} by 1.

\end{itemize}

