%% This file is to be included by latex in ug.tex and nemo.tex
%
% Chapter:  Miriad shell
%	Some \if to resolve miriad vs. nemo references.
%
%       With some acknowledgements to Alex Rudolph  (march 1991)
%
\myfile{miriad.inc}
\mylabel{s:miriad}

The {\tt miriad} front-end works on any simple terminal,
\index{miriad, front-end}
and is not restricted to a particular operating system or computer type.
To invoke\footnote[1]{Optional command line switches are summarized in
Table~\ref{t:miriad2}}
{\tt miriad}, type:

\begin{verbatim}
        > miriad
\end{verbatim}
The usual system prompt will be replaced by the {\tt Miriad\%} prompt:
\begin{verbatim}
        Miriad%
\end{verbatim}
and {\tt miriad} will read a file, {\tt lastexit}, with 
\index{lastexit} the values of all
parameters saved when you last exited {\tt miriad} (see EXIT below).
This file will be created/read in/from your current working directory.

We shall now describe all {\tt miriad} commands in more detail. A summary
of the commands, and their syntax, is given in Table~\ref{t:miriad1} at
the end of this section. Commands not known to {\tt miriad} are simply
passed to your host operating system. This means on VMS the command
{\tt DIR}, and on UNIX the command {\tt ls}, would still be usable
and valid commands\footnote[2]{UNIX aliases are not supported
in {\tt miriad}}.

\subsection*{INP, GO, and TASK} 

To inspect the inputs of a task, as well as to select the task, {\it e.g.}
{\tt histo} type
\begin{verbatim}
        Miriad% inp histo
\end{verbatim}
{\tt miriad} will 
show the parameters of the task along with the values, if any, 
previously set.  For example, if, the first time you run {\tt miriad}, you type:
\begin{verbatim}
        Miriad% inp histo
\end{verbatim}
{\tt miriad} will reply by writing:
\begin{verbatim}
        Task:   histo
        in       =
        region   = 
        range    =  
        nbin     = 
\end{verbatim}
and will replace the {\tt Miriad\%} prompt with a task prompt
\begin{verbatim}
        histo%
\end{verbatim}
indicating that you have chosen the task {\tt histo}.  You can also choose the
task {\tt histo} without using {\tt inp} by typing:
\begin{verbatim}
        Miriad% task histo
\end{verbatim}
{\tt miriad} will then replace the {\tt Miriad\%} prompt with the 
{\tt histo\%} prompt, but
the inputs will not be printed out.  In either way
a parameter can be set by typing any of:
\begin{verbatim}
        histo% in=gauss
        histo% set in gauss
\end{verbatim}
(the first form being preferred). In either case retyping:
\begin{verbatim}
        histo% inp
\end{verbatim}
will result in {\tt miriad} replying:
\begin{verbatim}
        Task:   histo
        in       = gauss
        region   =  
        range    =  
        nbin     =      
\end{verbatim}
Tasks are run either by typing {\tt go} {\it taskname} at 
the {\tt Miriad\%} prompt
(advisable only if you know you like the inputs) or by typing {\tt go} at the
{\it taskname}{\tt \%} prompt.  Thus, in the above example, typing:
\begin{verbatim}
        histo% go
\end{verbatim}
would result in the task {\tt histo} running with {\tt in=gauss} as the only
assigned paramater (all the other parameters were set to their default
values).  

Any task can be run regardless of the chosen task, merely by typing:
\begin{verbatim}
        Miriad% go itemize
\end{verbatim}
{\tt miriad} then executes the task {\tt itemize}
(using whatever input parameters it finds from the {\tt lastexit} file
or a previous run) and changes the default 
task and prompt to:
\begin{verbatim}
        itemize%
\end{verbatim}

\subsection*{GOB}

The {\tt gob} command does exactly the same as the {\tt go} command,
except it runs the task in the \index{background}
background (UNIX), or spawns the task\index{spawn} 
(VMS). In either way, {\tt miriad} is immediately ready for new commands,
although the output of the task just started may still return
output to the terminal.

\subsection*{HELP and ?}

When typed at the {\tt Miriad\%} prompt, the {\tt help}
command gives a list of all
the available tasks.  This feature is useful when you use {\tt miriad} for the
first time.  When typed at the {\it taskname}{\tt \%} prompt, the 
help command gives
information on the chosen task.  As with other commands, typing
\begin{verbatim}
        histo% help itemize
\end{verbatim}
displays the help file for {\tt itemize} but does not change the default 
taskname ({\it i.e.} {\tt histo}). \index{help, miriad}
To display the list of available tasks from a {\it taskname} {\it \%}
prompt, type {\tt help ?}.  

The {\tt ?} command lists general help about the {\tt miriad} 
interface, including
most of the information contained in this section.

\subsection*{EXIT, END, QUIT}

{\tt exit} and {\tt end} both exit {\tt miriad}, and, 
if any parameter values have been
changed from the previous time you ran {\tt miriad}, all parameter values are
saved in a file {\tt lastexit}.  
The parameter values from this file are read
in automatically when you next run {\tt miriad}.  
{\tt quit} exits {\tt miriad} without
saving the present parameter values in {\tt lastexit}.  This command is useful
if you have changed parameter values you do not wish to save.

\subsection*{UNSET and RESET}

It is often convenient to assign a parameter value to the default.  This
assignment is accomplished with the command {\tt unset}.  For example, if the
task {\tt histo} had been run with inputs:
\begin{verbatim}
        in       = gauss
        region   =  
        range    = 0,1
        nbin     = 10
\end{verbatim}
but you wanted to run histo with the default for nbin, you would type:
\begin{verbatim}
        histo% unset nbin
\end{verbatim}
Multiple parameters can be unset on the same line.  Thus, to assign both
{\tt range} and {\tt nbin} to their default values, you would type:
\begin{verbatim}
        histo% unset range nbin
\end{verbatim}
Then, typing:
\begin{verbatim}
        histo% inp
\end{verbatim}
will result in {\tt miriad} replying
\begin{verbatim}
        Task:   histo
        in       = gauss
        region   =  
        range    =  
        nbin     =
\end{verbatim}
The command reset can be used to unset all parameter values (not just
those of the present task selected); use with caution.

\subsection*{SAVE and LOAD}

As noted above, when you exit {\tt miriad}, all your inputs will be saved in
the file {\tt lastexit}.  However, many parameters are used by more than one
program and, at present, only one value of each parameter is allowed
throughout {\tt miriad} (global parameters).  The command {\tt save}
writes the current parameter values to a default file ({\it taskname}{\tt .def} 
if the task
{\it taskname} has been selected) or to a user specified file.  This default
file (or a user specified file) can be read in using the command {\t load}. 
Note that {\tt save} writes out {\bf all} parameter values, not just those for
the specific task chosen.
Example:
\begin{verbatim}
        histo% save
        histo% save histo1.pars
\end{verbatim}
In the first case the parameters are saved in a file {\tt histo.def}, whereas
the second case saves them in a user specified file {\tt histo1.pars}.


\subsection*{VIEW}

Rather than type in parameter assignments one at a time, it is possible
to invoke a text editor (EDT in VMS, vi or the editor given by the
{\bf EDITOR} \index{EDITOR, environment}
environment variable in UNIX) to edit the {\tt taskname.def} file.  As
usual, if no taskname is given, view will edit the {\tt taskname.def} file
indicated by the {\tt taskname\%} prompt.  If a taskname is specified, 
{\it e.g.}:
\begin{verbatim}
        histo% view itemize
\end{verbatim}
the existing 
file {\tt itemize.def} is edited, and the selected task is changed to
{\tt itemize}.  Note that if the {\tt itemize.def} file did not exists, 
one is created.

\subsection*{CD}

The {\tt cd} (change directory) command\footnote[1]{Note that the 
{\tt miriad} {\tt cd} command does not understand special UNIX symbols}
allows you to change your
working directory from within {\tt miriad}.  Note that the file {\tt
lastexit} is written to your current working directory, but that you will
exit {\tt miriad} to the directory you started it from. You can always
find out your current working directory with a local host command,
{\it i.e.} {\tt pwd} on UNIX and {\tt SHO DEF} on VMS.
{\it perhaps cd with no args should say current working dir}

\subsection*{VERSION}

The {\tt version} command displays the which version of {\tt miriad} you are
running, and how it was compiled\footnote[2]{This latter information is
useful only for the advanced user, so if you don't understand it,
don't worry.}  If you have a problem with {\tt miriad}, especially one you
have never seen before, check if the version of {\tt miriad} has
changed.  If so, you may want to tell you local MIRIAD manager about the
problem, if you suspect that it is related to the new version of {\tt miriad}.


%%%

%------------------------------------------------------------------------
%NEMO
% \subsection*{TCL}
% An experimental command processor may been installed in your version
% of Miriad. If so, the {\tt VERSION} command should tell you so. See
% \index{TCL, command language}
% {\tt "TCL: An Embeddable Command Language"} by {\it John K. Ousterhout}
% for more details how to use this.
%
%========================================================================

\subsection{Command and command line switches overview}

The two tables summarize all of {\tt miriad}'s commands
(Table~\ref{t:miriad1}) and the optional command line switches
(Table~\ref{t:miriad2}). As an example if you want
to keep all default files located globally accross the 
filesystem, independant of the project you're working on, you
could create a subdirectory in you home directory, and use that
when you startup miriad. A UNIX example:
\begin{verbatim}
    % mkdir ~/def
    % miriad -d ~/def
\end{verbatim}

\begin{table}
\caption[Miriad shell command overview]{Miriad shell command overview}
\mylabel{t:miriad1}
\centering
\begin{tabular}{|l|l|l|}\hline
Command   &  Syntax             & Comments \\ \hline
&& \\
=       &  par = value          & assignment (see also deprected SET) \\
!       &  !cmd par1 par2 ...   & run a cmd in parent shell\\
?       &                       & help \\
\^~C    &                       & attempt to interrupt current running task \\
alias   &  alias [name] [value] & set or show aliases \\
cd      &  cd directory         & change directory \\
exit    &  exit                 & exit miriad, and save parameters in lastexit \\
go      &  go [taskname]        & start up task \\
gob     &  gob [taskname]       & start up task in background \\
help    &  help [taskname]      & online help on task\\
inp     &  inp [taskname]       & overview inputs \\
input   &  input cmdfile        & read commands from a command file\\
load    &  load [keyfile]       & load task parameters from keyfile\\
quit    &  quit                 & quit miriad, and does not save lastexit \\
reset   &  reset                & reset all variables \\
save    &  save [keyfile]       & save task parameters to keyfile\\
set     &  set par value        & assigment ({\it deprecated}) \\
task    &  task taskname        & set new default taskname \\
unalias &  unalias al1 al2 .... & unset aliases \\
unset   &  unset var1 var2 .... & unset variables (blank them) \\
version &  version              & current version of miriad and capabilities\\
view    &  view [taskname]      & edit task parameters in a keyfile\\ \hline
\end{tabular}
\end{table}

\begin{table}
\caption{Miriad shell command line switches}
\mylabel{t:miriad2}
\centering
\begin{tabular}{|l|l|l|}\hline
Switch   &  Environment Variable & Comments \\ \hline
&& \\  
-f lastexit\_file & -         & global keyword file [{\tt lastexit}] \\
\index{lastexit}
-b mirbin\_dir    & MIRBIN    & executables \\
\index{MIRBIN, environment}
-d mirdef\_dir    & MIRDEF    & default files \\
\index{MIRDEF, environment}
-p mirpdoc\_dir   & MIRPDOC   & documentation files \\
\index{MIRPDOC, environment}
-g               &           & turn debugging on \\ 
                 & PAGER     & pager for reading documents [{\tt more}] \\
\index{PAGER, environment}
                 & MIRPAGER  & document pre-formatter [{\tt doc}] \\
\index{MIRPAGER, environment}
?                &           & help \\ \hline
\end{tabular}
\end{table}
%%


