%% This file is to be included by latex in nemo.tex
%
%  26-mar-94  for V2.2
%   7-apr-01  for V3.0
% Appendix: setup
\myfile{setup.inc}
\mylabel{a:setup}

This Appendix describe how you have to modify your UNIX environment in
order to use NEMO.  It requires nothing more than a few additions to
your standard {\tt .cshrc} startup file. 
\index{.cshrc} \index{.login} We do however
differentiate here between a {\bf static} and {\bf dynamic} setup.

Overall installation of the NEMO package will be discussed separately,
in Appendix~\ref{a:install} and requires some knowledge of the UNIX
operating system. 
 
\section{Static Setup}
\mylabel{s:staticins}

A {\bf static} setup provides NEMO automatically every time you log in
and may consist of up to four modifications to the {\tt .cshrc} file:
    
(1) The environment variable 
{\bf NEMO}\footnote{The environment variable {\bf NEMOPATH} 
is\index{NEMOPATH, environment} from V1 is now invalid}
should\index{NEMO, environment}
be defined to tell NEMO where its root directory is.
This environment variable 
is used to derive many subsequent references to locations
of data files, documentation, program binaries, libraries etc.
Add the following line to your
{\tt .cshrc} file before any other references to NEMO are made,
{\it e.g.}:

\begin{verbatim}
    setenv NEMO /usr/nemo
\end{verbatim}

(or whatever your local root directory for NEMO is).

(2) The optional environment variable {\bf NEMOHOST}\index{NEMOHOST,
environment} should be set
to refer to the HOSTTYPE\index{HOSTTYPE, environment} you want to run.
If not present, NEMO will attempt to resolve it. Some UNIX shells 
define an environment variable HOSTTYPE, which could also be used.

(3) Add the statement\index{NEMORC, file}\index{startup file, NEMORC}

\begin{verbatim}
    source $NEMO/NEMORC
\end{verbatim}

to your {\tt .cshrc} file. 

(3) Your search path should include 
{\tt \$NEMOBIN}, preferably in\index{PATH, environment}
the beginning of the path definition.
The reason for this specific location is that we often use
a slightly modified {\it cc(1)} compiler (script),
{\it e.g.}, you may then have something like \index{NEMOBIN, environment}

\begin{verbatim}
    set path=(. $NEMOBIN /bin /usr/bin /usr/local/bin )
\end{verbatim}

in your {\tt .cshrc} file. 

If you do this the very first time, make these modifications permanent for
your current terminal session, {\it e.g.}:
\begin{verbatim}
    % source .cshrc ; rehash ; echo $PATH
\end{verbatim}

\section{Dynamic Setup}
\mylabel{s:dynains}

A {\bf dynamic} setup provides NEMO only after a startup
load command, usually dubbed {\tt nemo}, issued during an
interactive terminal session. A script 
{\tt \$NEMO/nemo.rc}\index{nemo.rc, file}\index{startup file, nemo.rc}
is available for this purpose. The following two modification are then
necessary to your {\tt .cshrc} file, instead of the above described
procedure:

\begin{verbatim}
    setenv NEMO /usr/nemo
    alias nemo 'source $NEMO/nemo.rc'
\end{verbatim}

Note the single quotes, to allow so-called late evaluation.
Whenever the command ``{\tt nemo}'' is issued, the NEMO environment
is loaded with whatever the current value of \$NEMO is.
After this has been done, the {\tt nemo} alias is 
replaced by another one to prevent re-entry. The alias
{\tt omen}\index{omen} will unload NEMO from the environment.

\section{Tailoring}
\mylabel{s:tailoring}

If you create a file {\tt .nemorc}\index{.nemorc}  in your
home directory, it will be read after the standard  {\tt \$NEMO/NEMORC}
(or {\tt nemo.rc})
file, and this is where you would want to add your private NEMO additions,
or even override the things you don't want.
