%% This file is to be included by latex in nemo.tex
%
% Appendix: Potentials
\myfile{potname.inc}
\mylabel{a:potential}

This Appendix\footnote{Automatically generated from CTEX comments
in the {\tt \$NEMO/src/orbit/potential/data} source code}
lists a number of potentials that are available in NEMO.
Most NEMO programs that deal with potentials have three 
program keywords associated with 
potentials: {\bf potname=, potpars=} and {\bf potfile=},
describing their name, optional parameters and optional
associates filename (or other textual information).
Each section below details a potential, and 
explains the usage of the {\bf potpars=} and {\bf potfile=}
keywords. Mostly $G=1$, unless otherwise mentioned.



\input potctex.inc
