%% This file is to be included by latex in nemo.tex
%
% chapter: Shell Scripts
\myfile{shell.inc}
\mylabel{a:shell}

Shell scripts can make life a lot easier.  If you have never written
them, and have to do a lot of repetitive work, this appendix is for you. 
\index{shell script} It will give you a small tour through creating and
using them, but for a more complete coverage proper literature is
recommended.\footnote{recommended: Anderson and Anderson, {\it The
Unix C Shell Field Guide}.  (Prentice Hall,1986)}
In the second part of this appendix we will introduce you to an
advanced usage of shell script, to make them look more like
real NEMO programs with a ``keyword=value'' type interface.

\section{C Shell Scripts}

A shell script is a plain ASCII text file and can hence be created with
any text editor.  It is a list of commands which you would normally have
issued yourself on the commandline (quite similar in concept to command
files under VMS or batch files under MSDOS).  In addition, most shells
have the capability of command flow logic ({\tt goto}, {\tt
if/then/else/endif} etc.) retrieving the command line parameters,
defining and using (shell) variables etc. 

Under the Unix environment one can also choose which shell to use,
although we shall only give examples in the most commonly
used shell, the C-shell ({\tt csh}). As an example,
we will show a shell script which copies all files from one directory to
a new one, also creating that new directory. The new directory must not
exist yet, otherwise the script will fail with an error message. It's
usage would look like:

\begin{verbatim}
    % csh -f scriptname dir1 dir2 
\end{verbatim}
or shorter:
\begin{verbatim}
    % scriptname dir1 dir2
\end{verbatim}
The second form, in which the script is not told through which 
shell it should process its commands, is the recommended practice.
The first line of this script file must then contain the line
% This technique is called hash-bang, and may not be available
% on all Unix implementations. ???
\begin{verbatim}
    #! /bin/csh -f
\end{verbatim}
to denote the script is to be run via the C-shell, which on standard
Unix systems is located in {\tt /bin/csh}. In this case the script 
needs to be made executable, {\it i.e.}
\begin{verbatim}
    % chmod +x scriptname
\end{verbatim}

In effect your operating
system will then issue the first form of the command, The second form
has the advantage that you don't have to remember which shell to use,
and in the end saves a few keystrokes, always considered a big issue in 
Unix!

Now lets look at the full text of the script first:

\footnotesize\begin{verbatim}
    #! /bin/csh -f
    #
    # Example of a shell script to copy all files from one directory 
    # to another. The input directory must not contain any subdirectories, 
    # and it will not copy any so-called (hidden) dot-files.
    #
    ##                 	 check if called properly
    if ($#argv != 2) then
        echo "Usage: $0 dir1 dir2"
        echo "copies all files from one directory to another"
        goto done
    endif
    ##                  save command line args in variables
    set dir1=$1
    set dir2=$2
    ##                  check if dir1 indeed is an existing dir
    if (! -d $dir1) then
        echo "$dir1 is not a directory" ; exit 1
    endif
    ##                  check if dir2 does not exist
    if (-e $dir2) then
        echo "$dir2 already exists" ; exit 1
    endif
    ##                  create new dir2
    mkdir $dir2
    if ($status != 0) goto error
    ##                  loop through all files in dir1
    foreach file ($dir1/*)
        if (-d $file) then
            echo "Skipping $file (is a directory)"
        else
            echo "Copying $file"
            cp $file $dir2
        endif
    end
    ##           Labels to jump to exit OK (done) or not OK (error)
    done:
     exit 0
    error:
     exit 1
\end{verbatim}
\normalsize

A few things can be noted:

\begin{itemize}

\item Comments are lines that start with a {\tt \#}, but the first
 line of the script
 must contain this strange construction 
 {\tt "\#! /bin/csh"}\footnote{Some older versions of UNIX do not 
 understand this technique}
 to tell it how to execute itself. By default, if the first
 line would not contain such a directive, the script 
 would be executed by the (more primitive) Bourne shell ({\tt /bin/sh}).
 Any options that would normally be supplied to the shell must then
 be given in the first-line shell directive, as with the {\tt -f} option
 in the example above.

\item Shell variables can be created using the
{\tt set name=value} construct; they are henceforth
\index{set,csh}
referenced by prepending the shell variable name with a 
dollar: {\tt echo \$name} would simply display the {\tt value} of 
the {\tt name} shell variable ({\tt echo} is a UNIX command).

\item Command line arguments are in special 
shell variables {\tt 0, 1, 2} etc.
{\tt \$0} is the name of the script, {\tt \$1} the first argument
(if present), etc. All command line arguments 
{\tt \$1, \$2,.. } are also pre-defined in an shell variable
{\tt argv}, which is actually a list (like an array). It can be
referenced as {\tt \$argv[1], \$argv[2], ...} etc. This form has
the advantage that the construct \verb+$#argv+ (the same
as \verb+$*+) can be used to find the
number of elements in the list. In our example we want exactly two,
hence the first check. A shell variable list can be initialized using
brackets, {\it e.g.}: {\tt set name=(val1 val2 val3)}. In this case
\verb+$#name+ is 3, \verb+$name[2]+ has the value {\tt val2} etc.

\item The {\tt if (test) then/else/endif} is used
 to test a condition and control command flow. The {\tt !} is
used to negate a condition, and multiple conditions can be
used together in the {\tt or} ({\tt !!}) and {\tt and}
({\tt \&\&}) boolean operators. The {\tt -e} tests
if the next argument is a file or directory that exists, 
and the {\tt -d} if the next argument is a directory.

\item The {\tt foreach name (list)} construct takes all elements from
   the {\tt list} that follows, and executes all commands until 
    \index{foreach, csh}
    a matching {\tt end} is encountered. Inside the loop the current
    shell variable is defined in {\tt name}.
   
\item Labels are defined by a name, followed by a colon, 
    {\tt label:} jumping to a label is be done by using {\tt goto label}.
    \index{label, csh} \index{goto, csh}

\item Scripts should by terminated properly by an {\tt exit} command,
 optionally  returning an error code (0=no error) to the caller.
 This gives the caller the opportunity to catch faulty behavior.
 After each command the {\tt status} shell variable contains the
 {\tt exit} value ({\tt \$status}) of that command.

\item After each shell command the special shell variable {\tt status}
contains the exit status of the previous command. We have used it
once to check of the {\tt mkdir} command completed successfully,
if not, a direct exit by jumping to the label {\tt error} was
provided for.

\end{itemize}

\ifnemo
More examples of shell scripts in NEMO can be found in
the directory {\tt \$NEMO/csh} (public scripts) and for the
more adventurous in {\tt \$NEMO/src/scripts} (maintenance scripts).
\fi

\section{Parseargs C-shell scripts}
\mylabel{s:parseargs}

Writing user friendly shell scripts often requires a fair amount
of administrative work, especially if
the number of variables to the script is variable.

If the public domain program {\tt parseargs} \index{parseargs}
is available, that approach may used. 
Shell scripts created with this option have the
\index{plummer, model}
syntax (using the example {\tt plummer.csh} C-shell script below):
\begin{verbatim}
    % plummer.csh [-n <nbody>] [-s <seed>] <file>
\end{verbatim}
in the short notation\footnote{like the UNIX command-line syntax}
using single character options, or in the longer form:
\begin{verbatim}
    % plummer.csh [+nbody=<nbody>] [+seed=seed>] <file>
\end{verbatim}
In this longer form the difference between
the normal NEMO user interface is that keywords
have to be preceded with a {\tt +} symbol, plus it is possible to
create scripts with ``parameters'' without an associated keyword.
Perhaps this last option should not be used (this example uses
it in the {\tt <file>} filename parameter) in order to keep 
in pace with core NEMO programs.


%               Example script:
\newpage
Example script {\tt plummer.csh}:

\small\begin{verbatim}
#!/bin/csh -f
#
#  Example C-shell script, using parseargs, to create an Nbody Plummer model
#
set NAME="`basename $0`"

setenv ARGUMENTS "\
  '?', ARGHIDDEN, argUsage, NULL,    'Help {print usage and exit}', \
  'n', ARGOPT,    argInt,   nbody,   'Nbody {# bodies in disk}', \
  's', ARGOPT,    argInt,   seed,    'Seed {Random number seed}', \
  ' ', ARGREQ,    argStr,   file,    'File {Output Filename}', \
  ' ', ARGLIST,   argStr,   argv,    'argv {any remaining arguments}', \
  ENDOFARGS \
"
## set defaults ##
set nbody='1000'    ## default number of bodies in disk
set seed='0'        ## default inityial seed
set other=()        ## remaining parameter (ignored in here)

## parse command-line ##
parseargs -s csh -e ARGUMENTS -u -- "$NAME" $argv:q >/tmp/tmp$$
if ( $status != 0 ) then  ## improper syntax (or just wanted usage)
        rm -f /tmp/tmp$$
        exit 2
endif

## evaluate output from parseargs & remove temporary file
source /tmp/tmp$$
rm -f /tmp/tmp$$

## echo arguments ##
echo "ARGUMENTS:"
echo "=========="
echo file=$file:q
echo nbody=$nbody:q
echo seed=$seed:q
if ($#argv > 0) echo Additional Positional Parameters=$argv:q

## get down to the work
#
echo " ### Making initial conditions"
mkplummer out=$file nbody=$nbody seed=$seed
\end{verbatim}\normalsize
