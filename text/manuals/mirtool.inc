%% This file is to be included by latex in ug.tex
%
% mirtool - A SunView Front End
%
%       19-feb-91   updated for new usersguide from older..nemotool  PJT
%        7-apr-01   now deprecated in NEMO V3
%
\myfile{mirtool.inc}
\mylabel{s:mirtool}

For Sun workstations \index{Sun} running {\tt suntools} or {\tt sunview}
the program {\tt mirtool} aids \index{mirtool, front-end}
users in running 
\ifnemo 
NEMO or
\fi 
MIRIAD programs in a window environment. 

{\tt mirtool} itself is batch oriented; {\it i.e.}, it sets up, runs,
and controls the tasks, but does not interact directly with them.
% {\it e.g.} the {\tt AIPS} or {\tt GIPSY} \index{AIPS} \index{GIPSY}
% interfaces.
Tasks which require user interaction must open a window
themselves. Since this is generally not done, except for some graphics
programs, programs that expect input from the keyboard may hence act 
up\footnote[1]{In technical terms: {\it stdin} is closed}.

To invoke\footnote[2]{Optional command line switches are summarized in 
Table~\ref{t:mirtool1}} {\tt mirtool}, type:
\begin{verbatim}
    % mirtool &
\end{verbatim}
%%% BEGIN NEMO only
\ifnemo
for MIRIAD users, and
\begin{verbatim}
    % mirtool -b $NEMOBIN -p $NEMODOC &
\end{verbatim}
for NEMO users. An alias {\tt nemotool} \index{nemotool} is available for
these purposes.
\fi
%%% BEGIN NEMO only
in which the trailing \verb+&+ is used to put the program in the
background. 
After this a window appears containing a few subwindows with
confusing stew of subwindows with
buttons, labels, inputs, and messages.  
At the top of the window is the main control panel where
all of the important bits of information hang out. A few buttons
are present: they are activated (or toggled) with the LEFT BUTTON
on the mouse.
The input items appear as a name printed in boldface, ending
in a colon after which a value can be filled in by hand or by 
pulling down a menu using the RIGHT BUTTON (examples are
{\bf name} and {\bf logfile}).
Clicking the LEFT BUTTON
on these name items also activates the cursor in the panel, such
that a value can be typed in. The MIDDLE BUTTON will attempt to
give help (be careful with the {\bf Help} button though - see below). 

% \ifnemo
% \begin{figure}
%  \vspace{6.5in}
% \caption[Example of NEMO'`s user interfaces on a Sun workstation]
% {Example of NEMO's user interfaces on a Sun workstation. The lower 
% left window uses the {\tt miriad(nemo)} shell, the two windows on the 
% right are created by {\tt mirtool}.}
% \end{figure}
% \fi

The {\bf Name} blank is the most important input in the window, it is
where you type in the name of the task to be run.  If you forgot
the list of available tasks, or merely dislike typing, then
clicking the RIGHT MOUSE button will produce a menu displaying all the
currently available tasks.  Releasing the mouse button with any taskname 
highlighted
will result in that taskname being loaded in just as if it were 
typed in {\bf name}. Subsequently the names of the
parameters and their values from a previous run (if present) are
loaded into the middle left window, just below the {\bf name} input item.

The {\bf Quit} button quits {\tt mirtool}
(no surprise).  Jobs which are currently running
will continue to do so.  These jobs will also continue to send their
output to the specified log file without mirtool's intervention, since
there is no connection anymore between mirtool and the spawned
task.

The {\bf Save} button will save the current parameter values for
the current task in the appropriate keyword \index{keyword, file}
file.  These values are
automatically saved when a task is run, but not when the user quits or
switches to another task. 

The {\bf Run} button starts execution of the task specified in {\bf name}.
In the right window it will then display the job number, status and name.

The {\bf Help} button loads the general mirtool help file into
the help window. This only works when the help window has been
activated using the MIDDLE BUTTON on the current task. This is more
a bug than a feature. To get help on a particular task, it must
have been selected via the {\bf Name} blank, after which the MIDDLE
BUTTON must be clicked on the name or any of the keyword blanks in
the entry subwindow below. \index{help, mirtool}

The {\bf Local} button is actually a toggle button. When clicking on it,
{\bf Remote} will appear, and jobs will then
be run on the {\bf hostname} (see below) machine. This feature
has been disabled on most implementations right now.

The {\bf hostname} blank is the name of a remote accessible computer where
the job should be run.  The user must be able to execute the "{\tt rsh}"
command on this remote host.  For now, NFS must be used to share disks,
as the {\tt rsh} command needs to change directory to the one listed in the
{\bf directory} blank (see below). \index{remote execution}
A list of available hosts can be obtained by clicking on the
RIGHT MOUSE button on the {\bf hostname} blank.
This will produce a menu, like the task-{\bf name}
menu, of names which can be loaded into the hostname blank by releasing
the mouse button on a particular name. Again, this option has been disabled
on most implementations right now.

The {\bf directory} blank is the default directory where the tasks run. 
This is here so that the user can specify files and data sets as if they
were local; {\it i.e.}, instead of using {\tt /usr/nemo/Dataset}, the
{\bf directory} can be set to {\tt /usr/nemo}, and the data set,
usually as one of keywords of a task, can be
specified by {\tt Dataset}.  At startup it defaults to the directory
where mirtool was started up. 

The {\bf logfile} blank is the name of the current log file where output
from executed tasks should go.  Once the logfile is set, only tasks
which are run thereafter will use this logfile -- past tasks will
\index{logfile, mirtool}
continue to use the previous log file, or throw away the output, if
there was none.  Clicking the RIGHT BUTTON on the {\bf logfile}
blank will cause a logfile to open up (if a filename has been
entered);
({\it i.e.} a new window appears somewhere at the bottom).
The logfile can be hidden, erased, printed, and stop/started with 
the pulldown menu which appears thereafter using the RIGHT BUTTON.
The logfile should preferably be specified from the root directory,
{\it i.e.}, starting with a ``/'', since the interaction between the 
logfile location and the current directory can be strange. 

As mentioned before, the right subwindow in mirtool lists the current 
running (R) and stopped/paused (S) jobs.  With the RIGHT MOUSE button and a
pulldown menu these jobs can be temporarily stopped (STOP) started again
(CONT), or killed (KILL).  If the menu program is exited and re-entered,
these jobs here will not re-appear.  (A bug is present in the sense that
the status may not be correct when jobs were restarted from pausing)

In order to run mirtool, the following environment variables should be
set (see Appendix~\ref{a:setup}).
{\bf MIRBIN} \index{MIRBIN, environment}
is the directory where the binaries are found. 
{\bf MIRPDOC} \index{MIRPDOC, environment} is the
directory where the doc files "*.doc" are found.  
{\bf MIRDEF} \index{MIRDEF, environment} is the
directory where the users saved default (keyword) files go 
(the user must have local write permission).  
A handy default is the current directory {\tt "."}.  {\bf
MIRLOG} \index{MIRLOG, environment}
is the path of the users log file, if one is desired.
It can be changed interactively though by using the {\bf Logfile:} 
blank. 

You can also put something like the following line 
into your {\tt .rootmenu} file\footnote[1]{This file is used by
sunview to determine what can be selected from the root pulldown
menu.}:
\index{.rootmenu}
\begin{verbatim}
  "mirtool" /usr/local/bin/mirtool >& /dev/console
\end{verbatim}
Of course you will have to tailor the exact location of the {\tt mirtool}
program, since environment variables cannot be used in the {\tt .rootmenu}
file. Consult the Sun manuals for more help on operating {\tt sunview}.

\begin{table}
\caption{Mirtool shell command line switches (SunOS sunview)}
\label{t:mirtool1}
\centering
\begin{tabular}{|l|l|l|}\hline
Switch   &  Environment Variable & Comments \\ \hline
&& \\
-b mirbin\_dir    & MIRBIN    & executables \\
\index{MIRBIN, environment}
-d mirdef\_dir    & MIRDEF    & default files \\
\index{MIRDEF, environment}
-p mirpdoc\_dir   & MIRPDOC   & documentation files \\
\index{MIRPDOC, environment}
-l mirlog\_file   & MIRLOG    & filename of log \\
\index{MIRLOG, environment}
-g               &           & turn debugging on \\ 
                 & MIRPAGER  & document pre-formatter [doc] \\ 
\index{MIRPAGER, environment}
?                &           & help \\ \hline
\end{tabular}
\end{table}

