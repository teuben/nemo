%    -*- LaTeX -*-     for inclusion in nemo.tex
%
% \section{Orbits}


\myfile{orbit.inc}
\mylabel{s:orbit}


In this section we will describe how to integrate individual
stellar orbits, display and analyze them. Be aware that
although 3D orbits can be computed the number of utilities to
analyze them is rather limited.

Orbits are normally stored in datafile (see also
{\it orbit(5NEMO)}), and a close conceptual relationship exists
between a (single-particle type) {\bf snapshot} and an {\bf orbit}:
an orbit is an ordered series of phase-space coordinates
whereas a snapshot is a series of particles with no particular
order, but all at the same time.

Since orbits will be computed in an analytical potential, we assume for
the remainder of this section that you have familiarized yourself with
how to supply potentials to orbit integrator programs. They all share
the same triple ``{\tt potname=, potpars=, potfile=}'' keyword
interface, as described in Section \ref{s:potential}. Many
examples of the tricky {\tt potpars=} keyword are given in Appendix
\ref{a:potential}.

\subsection{Initializing}

There are a few programs with which orbits can be initialized:

\begin{itemize}

\item
{\bf mkorbit} is the most straightforward program. You can give simply 
give it all 6 phase space coordinates, and an orbit file consisting of 
this one point is generated. It is also possible to give the potential
in which the particle is to move, and 5 phase space coordinates plus
the energy, or even 4 phase space coordinates and the energy
plus the total angular momentum or angular momentum along the Z axis
(for axisymmetric systems).

Let's start with an example of creating a simple orbit by itself
with no associated potential.

\begin{verbatim}
% mkorbit out=orb1 x=1 y=0 z=0 vx=0 vy=0.2 vz=0
### Warning [mkorbit]: Potential potname= not used; set etot=0.0
pos: 1.000000 0.000000 0.000000
vel: 0.000000 0.200000 0.000000
etot: 0.000000
lz=0.200000                              

% tsf orb1
char History[59] "mkorbit out=orb1 x=1 y=0 z=0 vx=0 vy=0.2 vz=0 VERSION
  =3.2b"
set Orbit
  set Parameters
    int Ndim 03
    double Mass 1.00000
    double IOM[3] 0.00000 0.200000 0.00000
    int Nsteps 01
  tes
  set Potential
  tes
  set Path
    double TimePath[1] 0.00000
    double PhasePath[1][2][3] 1.00000 0.00000 0.00000 0.00000 0.200000
      0.00000
  tes
tes                                          
\end{verbatim}

\item
{\bf perorb} is a program that for given initial conditions
(similar to the ones described in {\tt mkorbit} above) attempts
to calculate periodic orbits in that potential. The output file
will be a file with one (or more) orbits. This is a bit of an
advanced program, and will not be covered here.

\item
{\bf stoo} is a program that can take a particle position from
a snapshot, and turn it into an orbit. For example, sampling some
initial conditions from the positions of stars in a Plummer sphere,
we could use the following small C-shell code to find some
statistical properties of this selected set of 
orbits\footnote{For the careful reader:
{\tt mkplummer} and {\tt potname=plummer} actually
have different units, and as such this experiment is not 
properly set up.}

\begin{verbatim}
    mkplummer out=p100 nbody=p100
    foreach i (`nemoinp 0:100:10`)    
        stoo in=p100 out=orb$i ibody=$i
        orbint orb$i orb$i.out 10000 0.01 10000 potname=plummer
        orbstat orb$i.out
    end
\end{verbatim}

The reverse program, {\tt otos} turns an orbit into a snapshot, and 
may come in handy since the snapshot package has far more advanced
analysis programs.

\end{itemize}

\subsection{Integration}
\begin{itemize}
\item
{\bf orbint} integrates orbits from given initial conditions. If the
input orbit has more than 1 step, the last step is taken as the
initial conditions. Although the {\tt potname=, potpars=, potfile=}
keywords can be given, if the input orbit contains...

\begin{figure}[htb]
    \begin{center}
    \leavevmode
    \epsfysize=8cm
    %TODO \epsfbox[90 90 510 510]{orb1.ps}
    \end{center}
\caption{Sample orbit 1 ({\tt orb1.out})}
\mylabel{f:orbit1}
\end{figure}


\item
{\bf perorb} finds periodic orbits, and stores a full period which should
close the orbit. This program finds periodic orbits in the XY plane
(i.e. currently it will only find 2D orbits) by searching for the
centers of invariant curves in the surface of section.

\item
{\bf henyey} also finds periodic orbits, but uses Henyey's 
method\footnote{see also van Albada \& Sanders, (1982, MNRAS, 201, 303)}.
This program has however not been released to the public version of
NEMO.

\end{itemize}

\subsection{Display}

\begin{itemize}
\item
{\bf orbplot} is the only orbit plotting program we currently have.
For more sophisticated display {\tt tabplot} and/or
{\tt snapplot} would have to be used after transforming the data.
Also {\tt snapplot} uses the powerful {\it bodytrans} expression
parser to plot arbitrary
body related expressions, although {\tt orbplot} can
handle both {\tt x, y, z} and {\tt vx, vy, vz} for the
{\tt xvar=} and {\tt yvar=} keywords. An example of the output of
{\tt orbplot} is given in Figure \ref{f:orbit1}.

\end{itemize}

\subsection{Analysis}

\begin{itemize}

\item
{\bf orbstat} is an example of a simple program that reads orbits,
and displays statistics of it's 2D (x-y-) coordinates: 
maximum extent, as well as statistics of the angylar momentum.
This program is not suited for 3D orbits yet.

\begin{verbatim}
% orbint orb1 orb1.long 
% orbstat orb1.out
# T     E       x_max   y_max   u_max   v_max   j_mean  j_sigma
1000 -0.687107 1 0.999958 0.746764 0.746611 0.2 3.83111e-09
\end{verbatim}

\item
{\bf orbfour} performs a variety of fourier analysis on the

\begin{verbatim}
% orbint orb1 orb1.long 100000 0.01 10000 10 plummer 
INIPOTENTIAL Plummer: [3d version]
Pattern speed=0
0.000000 0.020000 -0.707107     -0.6871067811865
100.000000 0.277794 -0.964901     -0.6871067811856
200.010000 0.020912 -0.708019     -0.6871067812165
300.020000 0.271222 -0.958329     -0.6871067812194
400.030000 0.023376 -0.710483     -0.6871067812465
500.040000 0.259253 -0.946360     -0.6871067812551
600.050000 0.027415 -0.714522     -0.6871067812765
700.060000 0.242979 -0.930086     -0.6871067812904
800.070000 0.033056 -0.720163     -0.6871067813065
900.080000 0.223694 -0.910801     -0.6871067813241
Energy conservation: 2.00138e-10                               
% orbfour orb1.long amode=t
<R> N A0 A1 A2 A3 A4 B1 B2 B3 B4
1 10001 0.000360461 0.334714     0.000150399 -0.000472581 -0.000158864
                   -0.000667155  0.000228086 -0.000725406  0.000103029

% orbfour orb1.long amode=f
<R> N C0 C1 P1 C2 P2 C3 P3 C4 P4
1 10001 0.000360461 
        0.334715      -0.114202 
        0.000273209   56.5992 
        0.000865763 -123.083
        0.000189349  147.035
\end{verbatim}

\item
{\bf orbsos} computes surface of section coordinates. Since this program
does not plot, but produces a simple ascii table, you can pipe the output
into {\tt tabplot}:

\begin{verbatim}
% orbsos orb1.long y | tabplot - 3 4  xlab=Y ylab=VY
% orbsos orb1.long x | tabplot - 3 4  xlab=X ylab=VX
\end{verbatim}

will plot either a Y-VY or X-VX surface of section.

\begin{figure}[htb]
    \begin{center}
    \leavevmode
    \epsfysize=8cm
    %TODO \epsfbox[90 90 510 510]{orb1sos.ps}
    \end{center}
\caption{Surface of Section for sample orbit 1 ({\tt orb1.long})}
\mylabel{f:orbit2}
\end{figure}


\item
{\bf orbdim}
computes the dimensionality of an orbit, i.e.  how
many integrals of motions it has.  Although it requires very long
integration times to accurately compute this, it is completely
automatic, and does not require an analysis like that for a surface of
section (which is also graphic).  It is based on an interesting paper
by Carnevali \& 
Santangelo\footnote{Carnevali, P. \& Santangelo, P., 1984. ApJ 281 473-476}. 


\item
{\bf otos} transforms an orbit back into a snapshot, thereby giving you
the much richer set of analysis tools that are available for
{\it snapshot}'s.

\end{itemize}



