%  This file is a LaTeX file, to be included by nemo.tex
%
%  Glossary:
%
\myfile{glossary.inc}
\mylabel{a:glossary}

\begin{itemize}

\item[{\bf bodytrans}] Dataformat that is used to perform\index{bodytrans}
arbitrary operations on expression variables used in snapshot's.

\item[{\bf ccd}] Synonymous for image; most programs in NEMO which handle
images start or end with ``{\tt ccd}''.

\item[{\bf fie}] Most expressions
that you give to program 
keywords\index{keywords, program}\index{program keywords}
are parsed by {\it nemofie}
and eventually {\it fie}. (Nomenclature borrowed from\index{GIPSY} GIPSY)

\item[{\bf fits}] ``Flexible Interchange Transport System'', a standard
dataformat used to interchange data between machines. Commonly used for
images.\index{FITS}

\item[{\bf history}] Each NEMO dataset normally contains a
data history in the form of of history items at the 
beginnging of the dataset.\index{hisf}\index{history, data}
They are normally kept track of by the data processing programs,
and can be displayed with the program {\tt hisf}.

\item[{\bf image}] Dataformat in NEMO, used to represent 2- and 3-D
images/data cubes. See also {\bf ccd}.

\item[{\bf miriad}] Another astronomical data reduction package, 
from which we have borrowed\index{miriad}
some user interfaces ({\tt miriad} and {\tt mirtool})
which are plug-compatible with our command-line syntax.

\item[{\bf orbit}] Dataformat in NEMO used to represent a stellar
orbit; most programs in NEMO which handle orbits start or end
with ``{\tt orb}''.

\item[{\bf potential}] Dataformat in NEMO used to represent a potential;
most programs in NEMO which handle potentials start or end
with ``{\tt pot}''.

\item[{\bf program keyword}] Keywords that are defined by the program
only. They can be seen by using the {\tt help=} keyword
(in itself being a system keyword).

\item[{\bf review}] A small user interface that pops up when a 
program is interrupted. Type {\tt quit} to exit it, or 
{\tt ?} for help. This feature of the user interface may
not be installed in your version.

\item[{\bf set}] Compound hierarchical data-structure of a
structured file. They are the equivalent of a C structure.

\item[{\bf snapshot}] Dataformat used in NEMO to represent 
an N-body system. Most programs that handle {\it snapshot}'s
in NEMO start or end with ``{\tt snap}''.

\item[{\bf structured file}] The binary data NEMO writes is in a
hierarchical structured format. Programs like {\tt rsf}, {\tt rsf}
and {\tt csf} perform general and basic I/O functions on such
files. They are hierarchical structured sets.

\item[{\bf system keyword}] Global keyword that every NEMO program
knows about, and are not listed in the (program) keywords that
can be seen by issuing e.g. {\tt help=} 
(in itself being a system keyword).

\item[{\bf yapp}] ``Yet Another Plotting Package'', the library definition
that is used by all programs that produce graphics output. It is kept
very simple. The {\tt yapp=} system keyword controls the
graphics device definitions/capabilities.

\end{itemize}

