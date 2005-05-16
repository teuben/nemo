%% This file is to be included by latex in nemo.tex
%
% Appendix: GSL
\myfile{gsl.inc}
\mylabel{a:gsl}

This Appendix describes some of the enhanced features if
NEMO was compiled with the\index{GSL}
GNU Scientific Library 
(GSL)\footnote{{\tt http://sources.redhat.com/gsl/}}.

\section{Installation}

During NEMO installation GSL should have been automatically
detected if present in one of a few common location
(e.g. {\tt /usr/local/include/gsl}), though the
configure flag
{\tt --with-gsl-prefix=} can also be used to force a different
location. Equally so, using 
{\tt --disable-gsl} can be used to avoid setting NEMO up to 
use GSL.
If not, and you install GSL in a specific location, you will
have to edit {\tt \$NEMOLIB/makedefs} and possible
{\tt \$NEMOLIB/config.h}, or re-install NEMO with the correct flags
or hope auto-detection will work.
\smallskip
To remind, installation of GSL can be done in two ways:

\begin{enumerate}
\item {\bf end user}:
tar
\begin{verbatim}
  tar zxf gsl-1.4.tar.gz
  cd gsl-1.4
  ./configure --with-prefix=$NEMO/opt
  make 
  make install
\end{verbatim}

\item {\bf developer} (or even maintainer):
\index{cvs, GSL}

\begin{verbatim}
  cvs -d :pserver:anoncvs@sources.redhat.com:/cvs/gsl checkout gsl
  ./autogen.sh
  ./configure --enable-maintainer-mode
  make
\end{verbatim}

\end{enumerate}


\section{Features}

The following components in NEMO make use of GSL:

\begin{enumerate}

\item xrandom: the random number generator can use any of the many
generator types that GSL supports. In classic NEMO we use the portable
(but slightly modified) {\it ran3} routine from Numerical Recipes, using the
standard {\bf seed=}\index{seed=, keyword}
but with GSL enabled, this keyword can now be used with an optional
random number generator type, viz. {\bf seed={\it seed-number},{\it generator-type}},
e.g. {\bf seed=0,ran3}. The default for GSL is their {\tt mt19937} generator,
though good old {\tt ran3} is also available. The NEMO program {\tt xrandom}
will help with testing and selecting a good generator.
\smallskip
The xrandom program will also display two new keywords:
{\tt gsl=} and {\tt pars=} with which different distributions can be
generated. This duplicates some of the {\tt gsl-randist.c} program in the
GSL distribution itself.

\item spline:

\item fitting:
This is under development

\item histogram:
There is {\tt gsl-histogram.c}, but this is not used yet in NEMO.

\end{enumerate}

\section{Future}




