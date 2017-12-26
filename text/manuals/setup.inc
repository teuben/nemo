%    -*- LaTeX -*-     for inclusion in nemo.tex
%
%  26-mar-94    for V2.2
%   7-apr-01    for V3.0
%  26-dec-2017  for V4.0, greatly simplified
%
% Appendix: setup
\myfile{setup.inc}
\mylabel{a:setup}



This Appendix describes how you have to modify your UNIX environment in
order to use NEMO.

Overall installation of the NEMO package will be discussed separately
in Appendix~\ref{a:install}, and requires some knowledge of your Unix
operating system.  Here we assume somebody has already installed NEMO
somewhere on the system for you. 
 
\section{Setup}
\mylabel{s:staticins}

The simplest would be something like

\begin{verbatim}
    source /local/pkg/nemo/nemo_start.csh
\end{verbatim}

if you are using a {\tt csh} type shell. This example is in the case when NEMO had
been installed in {\tt /local/pkg/nemo/}.  In case you were using the {\tt bash}
shell, this example would read

\begin{verbatim}
    source /local/pkg/nemo/nemo_start.sh
\end{verbatim}

If you don't know which shell you are using, the command

\begin{verbatim}
    echo $SHELL
\end{verbatim}

should tell you.

A word of caution: NEMO modifies a number of environment variables
(e.g. {\tt \$PATH}) in order to work. It can happen that this will
interfere with other packages which do the same.

If you want to automatically make NEMO part of your shell, you will need to add
the {\tt source} command to your {\tt .cshrc} or {\tt .bashrc} file, depending if
you use {\tt csh} or {\tt bash} resp.\footnote{Mac users need to use .{\tt bash\_login}}


\section{Tailoring}
\mylabel{s:tailoring}

If you create a file {\tt .nemo\_start.csh}
(or {\tt .nemo\_start.sh}) in your
home directory, it will be read after the standard setup,
and this is where you would want to add your private NEMO additions,
or even override the things you don't want.

\section{Interactive vs. Login shell}

A reminder on Unix shells might be in place here. The two shells
that most Unix shells are derived from are the orginal Bourne shell
({\tt sh}) and C-shell (\tt csh}).  Their modern equivalents are {\tt bash}
and {\tt tcsh}.

Shells are either a login shell or an interactive shell. Depending on
this nature, the shell will read one or more startup files. So placing
the appropriate {\tt nemo\_start.(c)sh} file can be daunting under
different Unix implementations (Linux and Mac seems to have taken a
different approach).



  
% interactive bash:  /etc/profile, ~/.bash_profile, ~/.bash_login, ~/.profile

% BASH: /etc/profile  ~/.profile  ~/.bash_login 1) ~/bash_profile
% CSH:  /etc/csh.login  /etc/csh.cshrc  ~/.cshrc ~/.login
