%% This file is to be included by latex in nemo.tex
%
% Appendix: Units and Coordinate Systems
\myfile{coordsys.inc}
\mylabel{a:coordsys}

\section{Coordinate Systems}

Astronomy is well known for its confusing coordinate systems: nature and math don't
always look at things through the same mirror. For example, the so-common (mathematical)
right handed coordinate system that we call X-Y-Z does not neatly
fit in with our own Galaxy, which
rotates counter-clockwise (meaning the angular moment vector points to the galactic
south pole). Sky coordinates (e.g. Right-ascension Declination) are is pinned on the
sky, looking up, instead of on a sphere, looking down, and become a left-handed coordinate
system where the "X" coordinate increases to the left. 

Here are some examples, and their respective NEMO programs that deal with this. See also
their {\it manual} pages for more detailed information.


\begin{itemize}

\item
{\tt mkgalorbit}: 
orbits in our galaxy have to deal with the UVW space velocities
in the galactic coordinate system. There is no formal definition of a spatial XYZ 
system, other than Z=0 being the galactic plane. So, where to put the sun if
the galactic center is (0,0,0) is purely by convention. It so happens that 
(-R0,0,0) is convenient since galactic longitude and latitude can be
easily expressed as $atan2(y,x)$ and $atan2(z,\sqrt{x^2+y^2})$ resp. 

\item
{\tt mkspiral, mkdisk}:
These programs use 
the {\tt sign=} keyword to set the sign of the angular moment vector. Positive
means thus counter-clockwise rotation in this convention. Of course with
tools like {\tt snapscale} and {\tt snaprotate} snapshots can always be
re-arranged to fit any schema.

\item
{\tt snaprotate}:
in order for a known object (e.g. a disk) to be viewed as an ellipse with given
position angle and ellipticity this program uses a series of Eulerian angle
rotations. Apart from the astronomical convention to using position angle
as a counter clockwise angle measured from north, these numbers do no trivially
convert to the angles we use on the projected sky.
There are some examples in the manual page, which we briefly highlight here. 

\begin{enumerate}

% Olling also uses -1 for clockwise and +1 for counter clock wise rotation

\item 
The first example describes the projection of disks of spiral galaxies. 
If kinematic information is also known, 
the convention is to use the position angle of the receding
side of the galaxy. The inclination still has an ambiguity, which could
be used to differentiate based on which is the near and far side, but this would
result in values outside the commonly used $0..90$ range. Thus the sense
of rotation (sign of the angular momentum)
is a more natural way, with
again {\tt sign=-1} for counter-clockwise rotating (as seen for us projected on
the sky)

\begin{figure}[htb]
    \begin{center}
    \leavevmode
    \epsfysize=5cm
    \epsfbox{galaxyrot.eps}
    \end{center}
\caption[Example galaxy disks]
{Example galaxy disks:
clockwise (M33, left) and counter-clockwise (M51, right), assuming trailing
spiral arms}
\mylabel{f:galaxyrot}
\end{figure}

Here are NEMO commands to create an example velocity fields of these two galaxies with
the right orientation and velocity field (with an arbitrary rotation curve of course):

\footnotesize
\begin{verbatim}

mkdisk - 1000 sign=+1 mass=1 |\                    # clock wise rotating
 snaprotate - - theta=-30,-160 order=yz |\
 snapgrid - vel-m33.ccd moment=-1

mkdisk - 1000 sign=-1 mass=1 |\                    # counter clock wise rotating
 snaprotate - - theta=+22,170  order=yz |\
 snapgrid - vel-m51.ccd moment=-1

\end{verbatim}
\normalsize

\item 
The second example is that of a barred galaxy, or two nested disks if you wish. Here
an addition angle, the difference between the major axis and that of the bar), is
a parameter.

\end{enumerate}


\item
{\tt snapgalview}

\item
{\tt ccdfits}


\end{itemize}

\section{Units}



