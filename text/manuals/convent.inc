%% This file is to be included by latex in nemo.tex
% 26-mar-94  for V2.2
%  Opening statements: Conventions used in this manual

\myfile{convent.inc}

\noindent
The following typographical conventions are used in this
manual:

Text in {\it italic}, such as {\it image(3NEMO)}, mean
a reference to a UNIX manual page. In this case the command

\begin{verbatim}
    % man 3 image
\end{verbatim}

\noindent
would bring up the manual page {\tt \$NEMO/man/man3/image.3}. See also
{\it man(1)}.

Text in {\tt verbatim} are used to display the contents of source files,
or sample interactive sessions. The latter come in a few categories,
where the system prompt denotes the system you're on. 
Currently you may see samples from a UNIX shell (csh), where the command
is preceded by the percent ({\tt \%}) symbol, and VMS DCL, where
the dollar ({\tt \$}) is used:

\begin{verbatim}
    % ls
    
    $ DIR
\end{verbatim}

Text in {\bf boldface} is used to denote UNIX  environment variables.


