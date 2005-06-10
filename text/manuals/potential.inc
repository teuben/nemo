% latex include; file goes into nemo.tex or so.
% From $NEMO/src/orbit/potential/data one ought to run
%		make ctex
% and use the resulting file 'potctex.inc' in the appendix.
% That will be the formal reference material for potentials,
% this file will be merely for inclusion in the main document,
% plus a figure.
%
\myfile{potential.inc}
\mylabel{ss:potential}

Here we list some of the standard potentials available in {\tt NEMO},
in a variety of units, so not always $G=1$!

Recall that most {\tt NEMO} program use the keywords {\tt potname=}
for the identifying name, {\tt potpars=} for an optional
list of parameters 
and {\tt potfile=} for an optional text string,for example
for potentials that need some kind of text file.
The parameters listed in {\tt potpars=} will always
have as first parameter the pattern speed in cases
where rotating potentials are used. 
A Plummer potential with mass 10 and
core radius 5 would be hence be supplied 
as: {\tt potname=plummer potpars=0,10,5}.  The plummer potential
ignored the {\tt potfile} keyword.

\bigskip
{\bf plummer:} Plummer potential (BT, pp.42, eq. 2.47)
$$
    \Phi = -  {  M  \over 
                {   {(r_c^2 + r^2)}^{1/2} }  }
$$
\begin{enumerate}
\item[$\Omega_p$] Pattern Speed
\item[$M$]      Total mass
\item[$r_c$]    Core radius
\end{enumerate}

%-----------------------------------------------------------------
%\begin{figure}[htp]
% \ifdoplot
%    \begin{center}
%        \leavevmode
%	\epsfxsize=12cm
%        \epsffile{plummer.ps}
%    \end{center}
% \else
%    \vspace{1.0in}
%    \begin{center}
%      {\it To generate this plot run ....}
%    \end{center}
%    \vspace{1.0in}
% \fi
% \caption[Plummer potential]{Plummer potential properties}
% \mylabel{f:plummer}
%\end{figure}
%-----------------------------------------------------------------
