%% This file is to be included by latex in nemo.tex
%
% chapter: ACS
\myfile{acs.inc}
\mylabel{a:acs}

In this section we shall discuss the interaction between
NEMO and the used code in {\it The Art of Computational Science}
(hereafter ACS)\index{ACS}, an
(online) series of books on how to build a computational lab.



\section{The data format}

The dataformat introduced in ACS is a very simple ascii-table
based one, in which a header denotes the number of particles
($N$) and the time ($T$), followed by $N$ rows of 7
columns in which the mass, positions and velocities in cartesian
coordinates are stores. The NEMO utility {\it tabtos} can read
these files, and with {\it snapprint} they can be created out
of NEMO's binary {\it snapshot} format:

\footnotesize\begin{verbatim}
  % mkplummer - nbody=4 | snapprint - m,x,y,z,vx,vy,vz header=t 
4 0
0.25 -0.380733 -0.639301 -0.123983 0.200564 -0.10937 -0.0599709
0.25 -0.463769 1.16927 -0.686732 0.085217 0.449426 0.153703
0.25 1.1554 -1.18444 0.810501 -0.037432 -0.263928 0.0507987
0.25 -0.310902 0.654469 0.000214173 -0.248349 -0.0761281 -0.14453
\end{verbatim}\normalsize

note that in this example the default outputformat 
was used in snapprint, resulting in low precision output. Use
{\tt format=\%.16g} to get the ideal format to preserve 
double precision accuracy.

\footnotesize\begin{verbatim}
  % sphere1 -n 4 | tabtos - - nbody,time mass,pos,vel | tsf -
seed = 1074484824
[reading 4 bodies at time 0.000000]
char History[48] "tabtos - - nbody,time mass,pos,vel VERSION=1.3j"
set SnapShot
  set Parameters
    int Nobj 4
    double Time 0.00000
  tes
  set Particles
    int CoordSystem 66306
    double Mass[4] 0.250000 0.250000 0.250000 0.250000
    double PhaseSpace[4][2][3] 0.470784 -0.750808 -0.0257877 0.00000
      0.00000 0.00000 -0.0537154 0.348303 -0.306566 0.00000 0.00000
      0.00000 0.423827 0.629383 -0.424370 0.00000 0.00000 0.00000
      -0.637115 0.0135874 -0.267421 0.00000 0.00000 0.00000
  tes
tes

\end{verbatim}\normalsize

\section{A first N-body integration}

Chapter 9 in ACS/Vol 1 (pp173) discusses how 
to setup a star cluster, in particular the program
{\it sphere1}.

Chapter 8 introduced an N-body code called {\it nbody\_sh}.

With the conversion programs it is now trivial to
convert data back and forth, including the usage of
pipes along the way. After this any analysis can be done
using NEMO programs:

\footnotesize\begin{verbatim}
% sphere1 -n 10 |\                                  ACS program creates a cold sphere or stars
      nbody_sh1 -o 0.25 -t 10 |\                    ACS program to integrate 
      tabtos - run1.snap nbody,time mass,pos,vel    NEMO program to convert to binary snapshots

% snapplot run1.snap nxy=4,4                        multiple panel XY projection
% snapxyz run1.snap - | xyzview -                   animation
% snapcenter run1.snap . report=t                   notice the center of mass

\end{verbatim}\normalsize