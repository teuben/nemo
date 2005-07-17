%% This file is to be included by latex in nemo.tex
%
% Appendix: Potentials
\myfile{potname.inc}
\mylabel{a:potential}

This Appendix
lists a number of potentials 
\footnote{Automatically generated from CTEX comments 
in the {\tt \$NEMO/src/orbit/potential/data} source code}
that are currently available in NEMO.
Most NEMO programs that deal with potentials have three 
program keywords associated with 
potentials: {\bf potname=, potpars=} and {\bf potfile=},
describing their name, optional parameters and optional
associates filename (or other textual information).
Each section below details a potential
and  explains the usage of the {\bf potpars=} and {\bf potfile=}
keywords. 
The section title is the actual {\bf potname=} to be used for
this potential.
Mostly $G=1$, unless otherwise mentioned.

A word of caution: {\it potentials} are being re-designed into 
a more general  concept of {\it accelerations}, although with the
intent to keep them backwards compatible.



\input potctex.inc
