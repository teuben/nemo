%    -*- LaTeX -*-     for inclusion in nemo.tex
% It also includes potential.inc
%
%   Chapter: examples, the start of the cookbook...
\myfile{examples.inc}
\mylabel{c:examples}

Now that we have a reasonable idea how NEMO is structured and used, we
should be ready to go through some real examples  Some of the examples
below are short versions of shell scripts\footnote{where applicable, the
examples in this chapter are written in the C-shell language}
available online in one of the directories
(check {\tt \$NEMO/csh} and perhaps {\tt \$NEMOBIN}).
The manual pages
{\it programs(8NEMO)} and {\it intro(1NEMO)} \index{programs(8)}
\index{intro(1)} are useful to find (and cross-reference) programs
if you're a bit lost. Each program manual should also 
have some references to closely related programs.


\section{N-body experiments}

In this section we will describe how to set up an N-body experiment, run,
display and analyze it.  In the first example, we shall set up a head-on
collision between two spherical "galaxies" and do some simple analysis.

\subsection{Setting it up}

In Chapter~\ref{c:filestr} we already used {\tt mkplummer} to create 
a Plummer model;
\index{Plummer, model}
here we shall use the program {\tt mkommod} ("MaKe an Osipkov-Merritt
MODel") \index{Osipkov-Merritt, models} \index{Merritt - see Osipkov}
to make two random N-body realizations of a King model \index{King, models}
with dimensionless central potential $W_c = 7$ and 100 particles each. 
The small number of particles is solely for the purpose of getting
results within a reasonable time. Adjust it to whatever you can afford
on your CPU and test your patience and integrator
(see Appendix~\ref{a:bench} benchmarks).

\small\begin{verbatim}
    1% mkommod in=$NEMODAT/k7isot.dat out=tmp1 nbody=100 seed=280158
              
\end{verbatim}\normalsize

These models are produced in so-called RMS-units \index{units, rms}
in which the
gravitational constant G=1, the total mass M=1, and binding energy E=--1/2.
In case you would like virial units
\footnote{Virial units are the preferred units, see also:
Heggie\index{Heggie D.} \& Mathieu\index{Mathieu R.}, E=--1/4,
in: {\it The use of supercomputers in stellar
dynamics} ed. Hut\index{Hut P} \& McMillan\index{McMillan S}, 
Springer 1987, pp.233}\index{units, virial}
the models have\index{virial, units}\index{rms, units}
to be rescaled using {\tt snapscale}:\index{snapscale(1)}

\small\begin{verbatim}
    2% snapscale in=tmp1 out=tmp1s rscale=2 "vscale=1/sqrt(2.0)"
\end{verbatim}\normalsize

In the case that your user interface was not compiled with the 
{\bf NEMOINP}\footnote{This can be found out by 
using the program nemoinp(1NEMO) or {\tt help=?}.}
directive, the {\tt vscale} expression has to be calculated by you,
{\it i.e.} {\tt vscale=0.707107}. Also note the use of the quotes in
the expression, to prevent the shell to give special meaning to
the parenthesis, which are shell {\bf meta} characters.
\index{meta, shell characters}

The second galaxy is made in a similar way\index{mkommod(1)}, with
a different seed of course:

\small\begin{verbatim}
    3% mkommod in=$NEMODAT/k7isot.dat out=tmp2 nbody=100 seed=130159
\end{verbatim}\normalsize

This second galaxy needs to be rescaled too, if you want virial units:


\small\begin{verbatim}
    4% snapscale in=tmp2 out=tmp2s rscale=2 "vscale=1/sqrt(2.0)"
\end{verbatim}\normalsize

We then set up the collision by stacking the two snapshots, albeit with
a relative displacement in phase space.  The program {\tt snapstack} was exactly
written for this purpose:\index{snapstack(1)}


\small\begin{verbatim}
    5% snapstack in1=tmp1s in2=tmp2s out=i001.dat \
                    deltar=4,0,0 deltav=-1,0,0
\end{verbatim}\normalsize

The galaxies are initially separated by 4 unit length and approaching
each other with a velocity consistent with infall from infinity
(parabolic encounter). The particles assembled in the data file
{\tt i001.dat} are now ready to be integrated.

To look at the initials conditions we could use:

\small\begin{verbatim}
    6% snapplot i001.dat xrange=-5:5 yrange=-5:5
\end{verbatim}\normalsize
which is displayed in Figure~\ref{f:encounter1}.

\begin{figure}[t]
 \ifdoplot
   % \PSinsert{encounter1.ps}{7.0}{6.0}{-1.0}{-1.0}{1.0}
   % \PSinsert{encounter1.ps}{7.5}{7.5}{0.5}{0.5}{1}{0}
    \begin{center}
        \leavevmode
        \epsfxsize=10cm
        \epsfysize=10cm
        %TODO \epsffile{encounter1.ps}
    \end{center}
 \else
    \vspace{1.0in}
    \begin{center}
      {\it To generate this plot run ....}
    \end{center}
    \vspace{1.0in}
 \fi
 \caption[Encounter: initial conditions]{Initial conditions for 
the encounter as set up in this section}
 \mylabel{f:encounter1}
\end{figure}

\subsection{Integration using hackcode1}

We then run the collision for 20 time units, with the standard
N-body integrator based on the Barnes 
\index{Barnes J.} \& Hut \index{Hut P.} "hierarchical tree" 
algorithm\footnote{see also their paper in: Nature, Vol. 324, pp 446 (1986).}:
\index{hackcode1(1)}


\small\begin{verbatim}
    7% hackcode1 in=i001.dat out=r001.dat tstop=20 freqout=2 \
       freq=40 eps=0.05 tol=0.7 options=mass,phase,phi > r001.log
\end{verbatim}\normalsize

The integration frequency relates to the integration timestep as {\tt
freq} = $1/\Delta t$, the softening length {\tt eps} = $\epsilon$, and
opening angle or tolerance {\tt tol} = $\theta$.  A major output of
masses, positions and potentials of all particles is done every {\tt
1/freqout} = 0.5 time units, which corresponds to about 1/5 of a
crossing time.  The standard output of the calculation is diverted to a
file {\tt r001.log} for convenience.  This is an (ASCII) listing,
containing useful statistics of the run, such as the efficiency of the
force calculation, conserved quantities etc.  Some of this information
is also stored in diagnostic sets in the structured binary 
output file {\tt r001.dat}. 

As an exercise, compare the output of the following two commands:
\index{tsf(1)}


\small\begin{verbatim}
    8% more r001.log
    9% tsf r001.dat | more
\end{verbatim}\normalsize


\subsection{Display and Initial Analysis}

\mylabel{s-dispanal}
As described in the previous subsection, {\tt hackcode1} writes various
diagnostics \index{diagnostics, N-body}
in the output file. A summary of conservation of energy and center-of-mass
motion can be graphically displayed using {\tt snapdiagplot}:
\index{snapdiagplot(1)}


\small\begin{verbatim}
    10% snapdiagplot in=r001.dat
\end{verbatim}\normalsize

The program {\tt snapplot} displays the evolution of the particle
distribution, in projection; \index{snapplot(1)}


\small\begin{verbatim}
    11% snapplot in=r001.dat
\end{verbatim}\normalsize

Depending on the actual graphics ({\it yapp}) \index{yapp}
interface of snapplot,
you may have to hit the RETURN key, push a MOUSE BUTTON or just
WAIT to advance from one to the next frame.

The {\tt snapplot} program has a very powerful tool built into it
which makes it possible to display any ``projection'' the user wants.

As an example consider:

\small\begin{verbatim}
    12% snapplot in=r001.dat xvar=r yvar="x*vy-y*vx" xrange=0:10 \
               yrange=-2:2 "visib=-0.2<z&&z<0.2&&i%2==0"
\end{verbatim}\normalsize

plots the angular momentum of the particles along the z axis,
$J_z = x*v_y - y*v_x$,
against their radius, $r$, but only for the even numbered particles,
(\verb+i%2==0+) within
a distance of 0.2 of the X-Y plane ($-0.2<z \&\& z<0.2$).
Again note that some of the expressions are within quotes, to prevent
the shell of giving them a special meaning. 

The {\tt xvar}, {\tt yvar} and {\tt visib} expressions are fed to the
C compiler (during runtime!) and the resulting object file is then
dynamically loaded \index{loadobj}
into the program for 
execution\footnote{loadobj, the dynamic object loader, 
does not works on all UNIX implementations}.
The expressions must contain legal C expressions and depending
on their nature must return a value in the context of the
program. {\it E.g.} {\tt xvar} and {\tt yvar} must return a 
real value, whereas {\tt visib} must return a boolean (false/true or
0/non-0) value. This should be explained in the manual page of the
corresponding programs.

In the context of snapshots, the expression can contain 
basic body variables which
are understood to the {\it bodytrans(3NEMO)} routine. \index{bodytrans(3)}
The real
variables {\tt x, y, z, vx, vy, vz} are the cartesian phase-space
coordinates, {\tt t} the time,
{\tt m} the mass, {\tt phi} the potential,
{\tt ax,ay,az} the cartesian acceleration and {\tt aux} 
some auxiliary information.
The integer variables are
{\tt i}, the index of the particle in the snapshot (0 being the
first one in the usual C tradition) and {\tt key}, another
spare slot.

For convenience a number of expressions have already been pre-compiled
(see also Table \ref{t:bodytrans}),
{\it e.g.} the radius {\tt r=$\sqrt{x^2+y^2+z^2}$=sqrt(x*x+y*y+z*z)},
and velocity {\tt v=$\sqrt{v_x^2+v_y^2+v_z^2}$=sqrt(vx*vx+vy*vy+vz*vz)}. Note that {\tt r} and
{\tt v} themselves cannot be used in expressions, only the basic
body variables listed above can be used in an expression. 

When you need a complex expression that has be used over and
over again, it is handy to be able to store these expression under
an alias for later retrieval. 
With the program {\tt bodytrans} \index{bodytrans(1)}
it is possible to save such compiled expressions object files under
a new name. 

%---------------------------------------------------------------------------
\newcommand{\rvec}{{\bf \={r}}}
\newcommand{\avec}{{\bf \={a}}}
\newcommand{\vvec}{{\bf \={v}}}
\small\begin{table}
\caption[Precompiled bodytrans expressions]{Some precompiled bodytrans expressions}
\mylabel{t:bodytrans}
\centering
\medskip
\begin{tabular}{|l|l|l|} \\ \hline
name        	&  type & expression \\ 
&& \\ \hline
%%
0	& int	& 0 \\
1	& int   & 1 \\
i 	& int 	& i \\
key	& int   & key (see also {\it real} version below)\\
0	& real  & 0.0 \\
1	& real  & 1.0 \\
ar	& real  & (x*ax+y*ay+z*az)/sqrt(x*x+y*y+z*z) \\
        &       & \ \  or: (\rvec$\cdot$\avec)/$|$\rvec$|$ \\
aux	& real  & aux \\
ax	& real	& ax \\
ay	& real	& ay \\
az	& real	& az \\
etot	& real 	& phi+0.5*(vx*vx+vy*vy+vz*vz) \\
        &       & \ \   or: $\phi$ + \vvec$^2$/2 \\
i	& real	& i \\
jtot	& real  & sqrt(sqr(x*vy-y*vx)+sqr(y*vz-z*vy)+sqr(z*vx-x*vz)) \\
        &       & \ \   or: $|$\rvec$\times$\vvec$|$ \\
key	& real	& key  (see also {\it int} version above)\\
m	& real	& m \\
phi	& real	& phi \\
r	& real  & sqrt(x*x+y*y+z*z) \\
        &       & \ \   or: $|$\rvec$|$ \\
t	& real	& t \\
v	& real	& sqrt(vx*vx+vy*vy+vz*vz) 
                  or: $|$\vvec$|$ \\
vr	& real  & (x*vx+y*vy+z*vz)/sqrt(x*x+y*y+z*z) \\
        &       & \ \    or: (\rvec$\cdot$\vvec)/$|$\rvec$|$ \\
vt 	& real  & sqrt((vx*vx+vy*vy+vz*vz)-sqr(x*vx+y*vy+z*vz)/(x*x+y*y+z*z)) \\
        &       & \ \   or: $\sqrt{}$(\vvec$^2$-(\rvec$\cdot$\vvec)$^2$/$|$\rvec$|^2$)\\
vx	& real	& vx \\
vy	& real	& vy \\
vz	& real	& vz \\
x	& real	& x \\
y	& real	& y \\
z	& real	& z \\ 
glon    & real  & $l$, atan2(y,x)*180/PI              [-180,180] \\
glat    & real  & $b$, atan2(z,sqrt(x*x+y*y))*180/PI  [-90,90] \\
mul     & real  & $(-vx\sin{l} + vx\cos{l})/r$ \\
mub     & real  & $(-vx\cos{l}\sin{b} - vy\sin{l}\sin{b}+vz\cos{b})/r$ \\
xait    & real  & Aitoff projection x [-2,2] T.B.A. \\
yait    & real  & Aitoff projection y [-1,1] T.B.A. \\
\hline
\end{tabular}
\end{table}\normalsize
%---------------------------------------------------------------------------

As usual an example: \index{bodytrans(1)}


\small\begin{verbatim}
    13% bodytrans expr="x*vy-y*vz" type=real file=jz
\end{verbatim}\normalsize

saves the expression for the angular momentum in a real
\index{angular momentum}
valued bodytrans expression file, {\tt btr\_jz.o} which can in future programs
be referenced as {\tt expr=jz} (whenever a real-valued bodytrans
expression is required), {\it e.g.}

\small\begin{verbatim}
    14% snapplot i001.dat xvar=r yvar=jz xrange=0:5
\end{verbatim}\normalsize

Alternatively, one can handcode a {\it bodytrans} function, compile it,
and reference it locally. This is useful when you have truly complicated
expressions that do not easily write themselves down on the commandline.
The $(x,y)$ AITOFF\index{aitoff projection} projection are an example of
this. For example, consider the following code in a (local working directory)
file {\tt btr\_r2.c}:


\small\begin{verbatim}
    #include <bodytrans.h>

    real btr_r2(b,t,i)
    Body *b;
    real t;
    int  i;
    {
        return sqrt(x*x + y*y);
    }
\end{verbatim}\normalsize

By compiling this:


\small\begin{verbatim}
    15% cc -c btr_r2.c
\end{verbatim}\normalsize

an object file {\tt btr\_r2.o} is created in the local directory,
which could be used in any real-valued bodytrans expression:


\small\begin{verbatim}
    16% snapplot i001.dat xvar=r2 yvar=jz xrange=0:5
\end{verbatim}\normalsize

For this your environment variable {\bf BTRPATH} \index{BTRPATH,
environment} must have been set to include the local working directory,
designated by a dot.  Normally your NEMO system manager will have set
the search path such that the local working directory is searched before
the system one (in {\tt \$NEMOOBJ/bodytrans}).

\subsection{Movies}

In the previous subsection we showed that the program {\tt snapplot}
\index{snapplot,movie}
displays \index{movie}
snapshots of the system at selected times. When
displaying such snapshots in rapid succession, the illusion of
a movie can be obtained.
Normally however, the time to search and read the data from disk, calculate
the required {\it bodytrans} projection and display it is much too
long to give this impression, except possibly for very small number of 
particles. Still one lacks basic frame manipulations.

Two solutions are offered: 

\subsubsection{snapplot, movie}

First, the snapplot program has a {\tt frame=} keyword, whose
value is the filename (more properly the base 
of the filename) of the saved bitmap\index{bitmap}
of the current image on the screen.  The format of this bitmap is system
and {\it yapp} interface dependent, {\it i.e.} it depends on which
version of snapplot is used.  In other words: make sure you use
the right compiled version of snapplot, and check your local
documentation.  It's name may even be as obscure as {\tt sp}, {\tt
snapplot\_sv} or so, but it may also be hidden under the official {\tt
snapplot} name itself. 

There should then exist a program which manipulates and plays the {\tt
frame} (-raster) files \index{rasterfile}
back at a high enough rate to call
it a movie.  This is normally
achieved by putting the frames into memory, with preferably a memory
bitmapped display device or a fast screen-loader.  On a SUN workstation
\index{SUN} the local format used is {\it rasterfile(5)},
and public domain
programs such as {\it movie(1NEMO)} and {\it movietool(1NEMO)}
can be used to re-display and
manipulate frame files.  A disadvantage with {\tt movie} is that one has
to be outside of suntools, though this may well be worthwile, because of
the large number of
options in the menu of {\tt movie}.  With the in-house written clone of
movie, called {\tt movie\_sv}, it is possible to display {\tt frame}
files from within suntools.  The menu of {\tt movie\_sv} is not as
sophisticated, and the number of frames which fit in memory is not as
large, but the advantage of doing it from within suntools is sometimes
more important. By choosing a smaller size of the rasterfile this
last disadvantage can even be circumvented. {\tt yapp=256} is the
default. The normal procedure\footnote{This assumes during installation
the {\tt snapplot} program was compiled with {\bf \$YAPP\_SV} and
renamed as {\tt sp}}
 to create frames and display them within
suntools is:

\begin{verbatim}
    17% sp snap.out frame=movie1
    18% movie_sv frame=movie1
\end{verbatim}

Combining the {\it bodytrans(3NEMO)} capabilities of {\tt snapplot}
with {\tt movie}, it is also possible to look at the 3D structure
of an N-body system. The program {\tt makepath} \index{makepath} and 
shell script utility {\tt 3dmovie}\index{3dmovie}\index{movie,3D}
generate various 'flyby's or
'fly-around's to get an idea of the 3D structure of a system, though this
is mainly intended for batch use, is cumbersome 
in interactive mode and better methods are available though (see below).

The current {\tt movie} programs for SUN workstations can only 
handle black and white rasterfiles. {\it check if this is a bug
or a feature. ** also some bug when running within openwindows **}

Another making movies: one can also use a fancy {\tt snapplot\_cg}
version, and record image by image on a high quality recording device,
or interface with a 
{\it genlock}\footnote{interface device between the 
video output of a image display
device and a VCR recorder which allows frame-by-frame recording} device,
or make slides.


\subsubsection{snap3dv, snapxyz, xyzview, xgobi}

Another approach is to store the computed 2D or 3D coordinates (that you
would normally view with {\tt snapplot} in\index{3D viewing}
a special format, and use an interactive 3D viewer for this.

The first step would be to create that intermediate format. Depending on the
viewer, and different format (or program) needs to be selected to do this.
We have currently two programs available: 
{\tt snap3dv}\index{snap3dv}: writes a number of formats, all intended
for exported outside of the NEMO environment (typically they are all ascii
files), and {\tt snapxyz}\index{snapxyz}: writes an {\bf xyzc} 
file (a regular NEMO structured file which you can view with {\tt tsf}).

As for the different ascii format, see the comments in the manual
page of {\tt snap3dv} which programs can be used.

As for the {\bf xyzc} format: the program {\tt xyzview}\index{xyzview}
can be used to interactively view any 3D coordinates of a snapshot
at different times. In addition, one can also store a 3D vector
information to the 3D coordinates (could be velocity, or force field or
magnetic field direction) of a point, and display these. {\tt xyzview} can
keep a limited number of point data in memory, and display them in
rapid succession, in the view/zoom/velocity mode that can be changed
dynamically. The viewing conditions can be stored in a file, which could
be used to create offline movies in a mode described above.

{\tt xgobi} is\index{xgobi}
is versatile multi-variate data browsing and analysis tool
which suits itself quite well to some interactive 3D viewing of 
N-body data.
\begin{verbatim}
    snapprint in=snap.out times=10 | xgobi
\end{verbatim}

You will need to obtain {\tt xgobi} via independant means, it does not come with
NEMO. A new version of xgobi is under development, see {\tt http://www.ggobi.org/}. 

\subsection{Advanced Analysis}

Determine position of the proper "dynamical center of mass" of a system,
{\it e.g.} using the 50\% most bound particles.  For this, it's best to
save the potentials of the particles during the integration.  Most
N-body integrators will have an option to let you save the potentials of
the particles.  In {\tt hackcode1} one has to use the {\tt
options=mass,phase,phi} keyword to save this relevant information.

Here is an example of analyzing the resulting output dataset of the
merger we ran earlier. Using {\tt snapcenter}\index{snapcenter}
the snapshot is centered weighing each particle
by the third power of the potential, 
followed by a computation of the lagrangian radii using
{\tt snapmradii}\index{snapmradii}.
The resulting table is piped into {\tt tabmath}\index{tabmath}
and the logarithm of all radii is computed (the first columns
contains the time) after which the converted table is piped into
tabplot.


\begin{verbatim}
    set rad=0.01,0.05:0.95:0.05,0.99
    set nrad=(`nemoinp $rad | wc -w`)
    @ nrad++

    snapcenter r001.dat - '-phi*phi*phi' report=f |\
        snapmradii - $rad |\
        tabmath - - %1,`nemoinp 2:$nrad format="log(%%%.f)" separ=,` all |\
        tabplot - 1 2:$nrad line=1,1 yapp=2/xs
\end{verbatim}



\subsection{Generating models}

Besides a variety of programs of the kind of {\tt mkplummer}, {\tt mkommod},
{\tt mkexpdisk} etc., 
models can also be generated by calculating appropriate tables
(containing a run of density, potential and radius)
and feeding them into programs which translate such tables into a 
snapshot. An example of a program in the making is {\tt anisot}.

An alternative way is to use a package such as {\tt Mathematica} to
integrate the differential equations. With the following example
\index{mathematica, tables} we \index{Plummer, mathematica table}
leave this as an exercise to the reader:  %% PIET did this

\small\begin{verbatim}
  << NumericalMath/RungeKutta.m

  rpsipsiprime = RungeKutta[{psiprime, -(2/r)psiprime + Exp[-psi]},
    {r, psi, psiprime}, {10^-5, 10^-10/6, 10^-5/3}, 50, 10^-8,
    MaximumStepSize->0.5];

  rpsi = rpsipsiprime[[Table[n,{n,Length[rpsipsiprime]}],{1,2}]];

  r = rpsi[[Table[n,{n,Length[rpsi]}],1]];

  psi = rpsi[[Table[n,{n,Length[rpsi]}],2]];

  rho = Exp[-psi];

  rrhotranspose={r,rho};
 
  rrho = Transpose[rrhotranspose];

  logrlogrho = Map[Log, rrho, {2}];

  c = Log[10.];

  log10rlog10rho = logrlogrho / c;

  shortrrho = Take[rrho,30];

  PlotODESolution[shortrrho, 1, 2, 
        AxesLabel-> {"r", "rho"}, 
        PlotLabel -> "Isothermal Sphere"]

  PlotODESolution[log10rlog10rho, 1, 2, 
        AxesLabel-> {"log r", "log rho"}]
\end{verbatim}\normalsize

\subsection{Handling large datasets}

On of NEMOs weaknesses is also it's strong point: programs must
generally be able to fit all their data in (virtual) memory.
Although programs usually free memory associated with data
that is not needed anymore, there is a very clear maximum
to the number of particles it can handle in a snapshot. By 
default\footnote{one can recompile NEMO in single precision and define
{\tt body.h} with less wastefull members}
a particle takes up about 100 bytes, which limits the size of a 
snapshots on workstations somewhat.

It may happen that your data was generated on a machine which had
a lot more memory\index{memory, large}\index{data, large}
then the machine you want to analyze your data on.
As long as you have the diskspace, and as long as you don't need
programs that cannot operate on data in serial mode, there is
a solution to this problem. Instead of keeping all 
particles in one snapshot, they are stored in several snapshots
of (equal number of) bodies, and as long as all snapshots have the same
time and are stored back to back, most programs that can
operate serially, will do this properly and know about it.
Of course it's best to split the snapshots on the machine with
more memory\index{snapsplit}:

\begin{verbatim}
    % snapsplit in=run1.out out=run1s.out nbody=10000
\end{verbatim}


If it is just one particular program (e.g. snapgrid\index{snapgrid})
that needs a lot of extra memory, the following may work:


\begin{verbatim}
    % snapsplit in=run1.out out=- nbody=1000 times=3.5 |\
        snapgrid in=- out=run1.ccd nx=1000 ny=1000 stack=t
\end{verbatim}


\section{Images}

In this section we will some examples of NEMO's image format.
The manual pages {\it image(5NEMO)}) describes the data format,
whereas {\it image(3NEMO)}) introduces image I/O library routines.

We will give some examples on how to create an image,
create a contour\index{contour diagram} diagram, and
export the image as a FITS file and use it within another package. 
Specific examples are also given how to read in that FITS file
in AIPS and IRAF.

Images in NEMO are stored (in memory as well as
disk) as double precision floating point 
numbers\footnote{one can re-install NEMO to work in single precision},
which limits 
programs to how large an image can be dealt with.
There are also the {\it xyio(3NEMO)} routines which allow 
row oriented access to the image data, but there are not
many programs who use these routines.

\input image.inc



\subsection{Making an image from a snapshot}


The simplest way to grid a snapshot into an image 
is {\tt snapccd}, but this program has 
been superseded by the much more powerful tool {\tt snapgrid}, which
is based on {\tt snapplot}. An example: \index{gridding,snapshots}
\index{snapccd(1)} \index{snapgrid(1)}


\small\begin{verbatim}
    19% snapgrid snap_in image_out xrange=-10:10 yrange=-10:10 \
        nx=128 ny=128 xvar=x yvar=z evar="m*m"
\end{verbatim}\normalsize

grids the $X$ and $Z$ coordinates of the snapshot {\tt snap\_in} to an
$128 \ast 128$ image {\tt image\_out}.  The range in gridded coordinates
is from -10 to 10 in both $X$ and $Z$, with pixel coordinates
defined in the center of a cell. Note that the\index{emissivity}
emissivity ({\tt
evar=m*m}) is given as the square of the mass, which could be applicable
for ionized hydrogen gas when the mass in the snapshot would have been
indicative of 
the gas density.  In the case of stars or neutral hydrogen gas {\tt
evar=m} would have been more appropriate (which is actually the default).

Often it is then desirable to smooth the image to improve the signal to noise
ratio; although the degree of smoothing\index{smoothing, image}
depends on
the average number of 'objects'
per pixel.  Example: \index{ccdsmooth(1)}


\small\begin{verbatim}
    20% ccdsmooth image_in image_out gauss=0.3
\end{verbatim}\normalsize

smooths the image to a 
circular beam with FWHM (Full Width Half Maximum) of 0.3 physical
units. In the 
above gridding example  this amounts to about 2 pixels.

{\tt snapgrid} can use any {\it bodytrans} expression for most of its
gridding and emission related variables. It
essentially computes the following two dimensional image:
$$
	I_p(x,y) = {{ \sum_{i=1}^{N} e_i z_i^p} \over {\Delta{x}\Delta{y}}}
$$
where a sum is taken over all $N$ particles along the line of sight,
in a pixel of size $\Delta{x}$ and  $\Delta{y}$
and a {\tt moment=}$p$ is computed in $z$ over
the emission $e$.
This is for {\tt mean=f}, typical in the case where
each particle represents a mass, and an image such as this would
conserve mass (a surface density is actually computed).
However, in the case where a particle
represents a density, it is more appropriate to ignore the cell size,
and compute a mean per cell
$$
	I_p(x,y) = {{ \sum_{i=1}^{N} e_i z_i^p} \over N }
$$
instead (the case {\tt mean=t}). For a galactic velocity field
one normally chooses $z=-v_z$, i.e. {\tt zvar=-vz}.

\subsection{Galactic and Extragalactic objects}

{\tt snapgrid} has a specific choice of defaults which would make
observers of extragalactic objects, 
{\it i.e.} external observers from the
positive Z axis, happy ({\it i.e.} {\tt xvar=x yvar=y zvar=-vz}).
However with the help
of a few other tools in NEMO one can also make galactic observers
happy.

First the extragalactic objects: \index{extragalactic, observing}
A total intensity map is generated with the defaults
arguments of snapgrid. Channel maps at specific 
velocities can be generated using {\tt snapgrid zrange=zmin:zmax} or
{\tt zrange=zmean,zsig}, depending on the required velocity beam.
A velocity map is also easy to generate:
The raw zeroth and first order moment maps are 
saved ({\tt snapgrid moment=0,1})
and smoothed ({\tt ccdsmooth}) after which they can be divided
({\tt ccdmath}) \index{ccdmath(1)}
resulting in a velocity map. \index{position-velocity maps}
Shortcuts are available for mean velocity and dispersion
using {\tt moment=-1} and  {\tt moment=-2} resp., though these modes
allocate extra memory for the additional images needed to perform
the operations inside of snapgrid.
Position-velocity diagrams
can be directly generated using {\tt snapplot/grid xvar=x yvar=-vz}. If 
position-velocity maps need to be smoothed, remember that it may have 
to be
done in two steps (independently in 'position' and 'velocity'), because
the current version of ccdsmooth can only do circular 
gridded beams if the beam
is two dimensional. A more detailed example is given in the next subsection.

{\tt snaprotate, snapshift, snapscale} are helpful tools to project a
galaxy before all these operations are performed. 

For galactic objects: \index{galactic, observing}
one must choose an internal viewpoint (x,y,z) and (vx,vy,vz) somewhere
``inside'' the object, and make this the new origin \index{snapspin(1)}
using {\tt snapspin} or {\tt snapshift}.
%% PJT CHECK: they both seemingly do the same thing, yet snapshift should
%% be the right one


A rotation ({\tt snaprotate}) \index{snaprotate(1)}
may also be necessary.
Then use {\tt snapplot/grid xvar=vr yvar='atan2(y,x)' evar='m/(x*x+y*y+z*z)'}, 
where the Y axis of the plot will be a longitude between $-\pi$ and $\pi$, 
and the X axis the radial velocity, to create the familiar position-velocity
diagram. The {\tt evar=} keyword is only needed in gridding the snapshot.

Pretty pictures can be obtained using {\tt ccdplot}, which can
combine a contour\index{contour} map, overlayed on a greyscale image. 

\begin{figure}[t]
 \ifdoplot
   % \PSinsert{encounter1.ps}{7.0}{6.0}{-1.0}{-1.0}{1.0}
   % \PSinsert{encounter1.ps}{7.5}{7.5}{0.5}{0.5}{1}{0}
 \else
    \vspace{1.0in}
    \begin{center}
      {\it To generate this plot run ....}
    \end{center}
    \vspace{1.0in}
 \fi
 \caption[Position-Velocity diagram of a galactic disk]{Position-Velocity 
diagram of a galactic disk as seen by an internal observer}
 \mylabel{f:posvel1}
\end{figure}

\subsection{Extragalactic velocity field}

In this example we shall make a velocity field of a a particle
representation of a disk galaxy, with stars on circular orbits in
centrifugal balance with a fixed background potential.  

For disk-stars on circular orbits the program {\tt mkdisk} is 
useful\footnote{Only relevant keywords are shown, remaining take
their default value, use {\tt help=} to see them}:

\small\begin{verbatim}
   21% mkdisk out=disk1 potname=expdisk potpars=0,1,0.5 rmax=2 mass=1
\end{verbatim}\normalsize

In this the mass of the disk was set to non-zero,
in order to assign a finite emission to each star later on. If you 
would plot it's configuration with a program like {\tt snapplot},
you would see
not only a more-or-less constant surface density but also
that the disk is infinitely thin ({\tt snapplot yvar=z}). Also the particles
are indeed on circular orbits ({\tt snapplot xvar=r yvar=vt}), 
and there is no velocity dispersion ({\tt snapplot xvar=r yvar=vr}).

Viewed from the positive Z axis (the default 
with {\tt snapplot xvar=x yvar=y}) 
we would see no radial velocities in the disk;
in order to get a realistic looking velocity
field, we would have to rotate the model around a line of nodes (say the X
axis) using {\tt snaprotate}:

\small\begin{verbatim}
    22% snaprotate in=disk1 out=disk1.r theta=60 order=x
\end{verbatim}\normalsize

To plot the radial velocity the program {\tt snapplot} is used by
assigning a different symbol ({\tt psize=}) to different radial
velocities:

\small\begin{verbatim}
    23% snapplot in=disk1.r psize=0.1*vz
\end{verbatim}\normalsize

To simulate a true observation we shall use the program {\tt snapgrid} 
to grid\index{gridding, snapshots}
the discrete snapshot data ({\tt x,y,z,vx,vy,vz}) from the file {\tt disk1.r} 
onto a CCD-like device: a rectangular pixel array (matrix), with a value
(brightness, velocity etc.) associated with each pixel.
Since we are interested in the radial velocity field the
zero-th and first order moment maps need to be obtained,
and divided to get a radial velocity field:
$$
    \left< v \right> = - { \int{I(z) V_z(z) dz } \over {\int{I(z) dz}} } \eqno(1)
$$
Here $I(z)$ and $V_z(z)$ are the intensity and radial velocity
along the line of sight. Note the extra $-$ sign, 
to conform to the astronomical convention that positive
velocity means negative {\tt vz} if viewed from the positive Z axis.
In NEMO the denominator and numerator in
eq. (1) are evaluated as follows:

\small\begin{verbatim}
    24% snapgrid in=disk1.r out=mom0 moment=0 zvar=-vz evar=m
    25% snapgrid in=disk1.r out=mom1 moment=1 zvar=-vz evar=m
\end{verbatim}\normalsize

Since the data will be noisy, it is best to smooth the data a bit. 
Smoothing must however be done before the maps are divided ({\it why?}). 
Since the default pixel size is 4/64=0.0625 a 
Gaussian beam with a FWHM of 0.15 is used
to convolve the data. We would
use the programs {\tt ccdsmooth} and {\tt ccdmath} in the following order:

\small\begin{verbatim}
   26% ccdsmooth in=mom0 out=mom0s gauss=0.15
   27% ccdsmooth in=mom1 out=mom1s gauss=0.15
   28% ccdmath in=mom0s,mom1s out=disk1.vel fie=%2/%1
\end{verbatim}\normalsize
 
The final output file, {\tt disk1.vel}, now contains the radial velocity 
field map at an inclination of $60^o$. It
can be displayed with programs like {\tt ccdplot} and {\tt ds}.
{\tt ccdplot} is a NEMO program, capable of plotting contours as well
as greyscale (if given the right graphics device driver). {\tt ds} is
general purpose image display program
and displays a map in color on a sun workstation (for this, {\tt ds} must
have been installed to understand the NEMO file format).


\begin{figure}[t]
 \ifdoplot
   % \PSinsert{encounter1.ps}{7.0}{6.0}{-1.0}{-1.0}{1.0}
   % \PSinsert{encounter1.ps}{7.5}{7.5}{0.5}{0.5}{1}{0}
 \else
    \vspace{1.0in}
    \begin{center}
      {\it To generate this plot run ....}
    \end{center}
    \vspace{1.0in}
 \fi
 \caption[Velocity field of a galactic disk]{Velocity field of a galactic disk}
 \mylabel{f:velfie1}
\end{figure}

You can also convert the NEMO image file to a FITS file. A FITS file is a true
astronomical standard, which can be read
into any other image processing package ({\tt ds} can also read FITS files)
({\tt AIPS, IRAF, MIDAS, Miriad}). Creating it can be done as follows:

\small\begin{verbatim}
    29% ccdfits in=velfie out=fits1
\end{verbatim}\normalsize

\subsection{Integrated Color Maps}

A true color map can be created from a snapshot by assigning a color
to the particles in the snapshot. This can be done by 
gridding the snapshots twice, each one assigning the particles with a 
different emissivity. In the example below a snapshot is assigned
an artificial radial color gradient. The particles are given a linearly
increasing emissivity, according to their ranking in radius. First,
we must make sure the snapshot is sorted in radius properly:

\small\begin{verbatim}
    30% snapsort in=snap.dat out=tmp1 rank=r
\end{verbatim}\normalsize

and next the snapshot is gridded twice:

\small\begin{verbatim}
    31% snapgrid in=tmp1 out=tmp1_1 evar=m
    32% snapgrid in=tmp1 out=tmp1_2 evar='m*(i+1)'
\end{verbatim}\normalsize

This means the color at the center would be $-2.5log(I_1/I_2) = 0$,
whereas at the edge the color would be $2.5log(nbody)$.
The images are best viewed when smoothed, and then divided
and taken the {\tt log} of. The factor 2.5 is left out here,
because the scaling is arbitrary:

\small\begin{verbatim}
    33% ccdsmooth in=tmp1_1 out=tmp1_1s gauss=0.3
    34% ccdsmooth in=tmp1_2 out=tmp1_2s gauss=0.3
    35% ccdmath in=tmp1_1s,tmp1_2s out=color_ccd.dat fie='log(%2/%1)'
\end{verbatim}\normalsize

........... {\it more to come}

\subsection{Extracting Rotation Curves from Galactic Velocity Fields}

As an example of image analysis we consider the extraction of a
rotation curve from an (axisymmetric)
disk galaxy. We shall assume the velocities
have been extracted already, but not consider the various 
tricky methods that exist to do this.

\smallskip
Two programs exist with which most scenarios can be played out to extract
rotation curves: {\tt rotcur} applies the tilted ring method, 
where a fixed rotation speed is assumed in a set of rings, of which
all geometric parameters (center, systemic velocity, position angle and
inclincation) can be either fitted or kept fixed at a given value. For
well behaving galaxies better signal to noise in the fitted parameters can
be achieved by the second program, {\tt rotcurshape}, which fits a shape
function to a disk. In both cases fitting occurs in a fully non-linear
sense, so initial values for all parameters need to be supplied. The
programs have not been written as to make reasonable estimates.

\smallskip
It is also worth noting that the output units in this program are
arcsec for radii, and km/s for velocities, since we are often using
real observations\footnote{note FITS files use degrees and m/s for 
distances and velocities}
For simulations you will need to use an appropriate
mnemonic or actual number for scaling in both distance and velocity
to get the output format with the appropriate precision. As an example,
if you use (Nbody based) virial units, you will most likely want to
use {\tt units=arcmin,100} or even {\tt units=1000,1000}, depending
on your taste.

\subsubsection{rotcur}

{\tt rotcur} applies the tilted ring method. This particular version has
been derived from the original sheltran version within GIPSY, although
many things still work the same compared to the 
current GIPSY version.

\small\begin{verbatim}
  1% rotcur help=h
in               : Input image velocity field [???]
radii            : Radii of rings (arcsec) []
vrot             : Rotation velocity []
pa               : Position angle (degrees) []
inc              : Inclination (degrees) []
vsys             : Systemic velocity []
center           : Rotation center (grids w.r.t. 0,0) [center of map] []
frang            : Free angle around minor axis (degrees) [20]
side             : Side to fit: receding, approaching or [both] []
weight           : Weighting function: {uniform,[cosine],cos-squared} []
fixed            : Parameters to be kept fixed {vsys,vrot,pa,inc,xpos,ypos} []
ellips           : Parameters for which to plot error ellips []
beam             : Beam (arcsec) for beam correction [no correction] []
dens             : Image containing containing density map []
tab              : If specified, this output table is used in append mode []
resid            : Output of residuals in a complicated plot []
tol              : Tolerance for convergence of nllsqfit [0.001]
lab              : Mixing parameter for nllsqfit [0.001]
itmax            : Maximum number of allowed nllsqfit iterations [50]
units            : Units of input {deg, arcmin, arcsec, rad, #},{#} for length and velocity [deg,1]
blank            : Value of the blank pixel to be ignored [0.0]
inherit          : Inherit initial conditions from previous ring [t]
reuse            : Reuse points from previous rings if used before? [t]
fitmode          : Basic Fitmode: cos(n*theta) or sin(n*theta) [cos,1]
nsigma           : Iterate once by rejecting points more than nsigma resid [-1]
imagemode        : Input image mode? (false means ascii table) [t]
wwb73            : Use simpler WWB73 linear method of fitting [f]
VERSION          : 2-jun-04 PJT [2.12]
\end{verbatim}\normalsize

As an example we shall consider the galaxy NGC 6503, for which VLA
data are publicly available on \index{ADIL}\index{NGC, 6503}
ADIL\footnote{NCSA :: ADIL code number 95.DW.01.01}.

\small\begin{verbatim}
   1% fits in=NGC6503.MOM1 - | ccdmath - vel1.ccd %1/1000
\end{verbatim}\normalsize

notice we're already converting the m/s in FITS to km/s in the {\tt ccd} file.
The angular distances are in degrees. In the printout below the output of
{\tt tsf} and relevant sections of {\tt fitshead} have been merged to
show their correspondence:

\small\begin{verbatim}
  set Parameters
    int Nx 512
    int Ny 512
    int Nz 1
    double Xmin 267.533                   // CRVAL1  =    2.67250000000E+02 /
    double Ymin 69.8322                   // CRVAL2  =    7.01166666667E+01 /
    double Zmin 1.00000
    double Dx -0.00111111                 // CDELT1  =     -1.111111138E-03
    double Dy 0.00111111                  // CDELT2  =      1.111111138E-03 
    double Dz 1.00000
    double Xrefpix 0.00000                // CRPIX1  =      2.560000000E+02
    double Yrefpix 0.00000                // CRPIX2  =      2.570000000E+02
    double Zrefpix 0.00000
    double MapMin -102.844                // DATAMIN =     -1.028444453E+05
    double MapMax 144.537                 // DATAMAX =      1.445368906E+05
    int BeamType 0
    double Beamx 0.00000
    double Beamy 0.00000
    double Beamz 0.00000
    char Namex[9] "RA---SIN"
    char Namey[9] "DEC--SIN"
    double Time 0.00000
    char Storage[5] "CDef"
    int Axis 0
  tes
\end{verbatim}\normalsize



\small\begin{verbatim}
\end{verbatim}\normalsize

\subsubsection{rotcurshape}

\small\begin{verbatim}
  1% rotcurshape help=h
in               : Input image velocity field [???]
radii            : Radii of rings (arcsec) []
pa               : Position angle (degrees) []
inc              : Inclination (degrees) []
vsys             : Systemic velocity []
center           : Rotation center (grids w.r.t. 0,0) [center of map] []
frang            : Free angle around minor axis (degrees) [0]
side             : Side to fit: receding, approaching or [both] []
weight           : Weighting function: {[uniform],cosine,cos-squared} [u]
fixed            : Geometric parameters to be kept fixed {vsys,xpos,ypos,pa,inc} []
ellips           : ** Parameters for which to plot error ellips []
beam             : ** Beam (arcsec) for beam correction [no correction] []
dens             : Image containing containing density map to be used as weight []
tab              : If specified, this output table is used in append mode []
resid            : Output of residual field []
fit              : Output the fit? or the residuals [f]
tol              : Tolerance for convergence of nllsqfit [0.001]
lab              : Mixing parameter for nllsqfit [0.001]
itmax            : Maximum number of allowed nllsqfit iterations [50]
units            : Units of input {deg, arcmin, arcsec, rad, #},{#} for length and velocity [deg,1]
blank            : Value of the blank (pixel) value to be ignored [0.0]
nsigma           : Iterate once by rejecting points more than nsigma resid [-1]
imagemode        : Input image mode? (false means ascii table) [t]
rotcurmode       : Full velocity field, or rotcur (r,v) fit only [f]
load             : dynamically loadobject file with rotcur_<NAME> []
rotcur1          : Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values []
rotcur2          : Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values []
rotcur3          : Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values []
rotcur4          : Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values []
rotcur5          : Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values []
VERSION          : 13-jan-05 PJT [1.3b]
\end{verbatim}\normalsize


\section{Tables}

NEMO has a few programs that manipulate table files,
\index{table, files}
although much more can be done 
with standard UNIX programs such as {\it awk(1)}. 
A few examples are given how they can be used together.  Programs such
as {\it mongo(1L)} \index{mongo} can be used to display results.  
If you're in the possession of {\tt sm}\index{sm}, most of the 
functionality of the table programs can be reproduced with {\tt sm}.
A quick and dirty plot can also be made with {\it tabplot(1NEMO)}. 
\index{tabplot(1)}

The following example shows how {\it tabmath(1NEMO)} 
and {\it awk(1)}\index{awk}\index{tabmath(1),awk}
can do the same thing:

\small\begin{verbatim}
    1% tabmath tab_in tab_out %1+%2
    2% awk '{ print $0,$1+$2 }' tab_in > tab_out
\end{verbatim}\normalsize

One can also use the table programs in UNIX pipes, and use NEMO's
feature of denoting a '-' (dash)  as filename for standard
input/output files:


\small\begin{verbatim}
    3% awk '{print $1,$3,$5}' tab_in | tabmath - tab_out \
                                "ifgt(%1,%2,sin(%3),cos(%3))" 
\end{verbatim}\normalsize

\index{awk}
Don't ask why such a complicated {\tt ifgt} construct, it's just
an example. See {\it nemofie(3NEMO)} for the syntax options
of the third ({\tt newcol=}) keyword of {\tt tabmath}. 

\subsection{Making an image from a table}

Any scatterdiagram \index{scatterdiagram, image}
can now easily be turned into an image by using the snapshot interface!
This would be temporary solution until the need for this would turn
into a program {\it tabccd(1NEMO)}.\index{gridding,tables}\index{tabgrid(1)}
The example below also demonstrates how
existing tools can be effectively combined to create a new
tool!

\small\begin{verbatim}
#! /bin/csh -f
#         - transform table into image - 
#       DEMO version: no bells and whistles
set in=$1               # infile (table)
set out=$2              # outfile (image)
set xcol=$3             # columns from table to use
set ycol=$4
set xrange=$5           # gridding area
set yrange=$6
set nx=$7               # number of pixels to use
set ny=$8
set sx=$9               # some smoothing 
set sy=$10

set tmp=tmp$$           # temp name for intermediate results

#     convert table to ASCII "205" snapshot (see atos(1NEMO))
awk  "END {print NR}"                $in     > $tmp.1
echo "3"                                    >> $tmp.1
echo "0.0"                                  >> $tmp.1
awk '{print 1.0}'                    $in    >> $tmp.1
awk '{print $'$xcol',$'$ycol',0.0}'  $in    >> $tmp.1
awk '{print 0.0,0.0,0.0}'            $in    >> $tmp.1
#                       convert to snapshot
atos $tmp.1 $tmp.2
#                       convert to image
snapgrid $tmp.2 $tmp.3 xrange=$xrange yrange=$yrange nx=$nx ny=$ny zvar=vz
#                       smooth image a bit
ccdsmooth $tmp.3 $tmp.4 gauss=$sx dir=x
ccdsmooth $tmp.4 $out  gauss=$sy dir=y
#                       write a FITS file
ccdfits $out $out.fits
#                       clean up mess
rm -f $tmp.*
\end{verbatim}\normalsize

\section{Potential}
\mylabel{s:potential}

Programs which need an external potential ({\it e.g.} orbit integrators)
can obtain these via the so-called potential descriptors.  They are
implemented in NEMO as loadable object files\footnote{Not all operating
systems allow the programmer to use this feature - see your local
loadobj \index{loadobj} implementation}.  To the user interface
this commonly
appears as a set of three program keywords {\tt potname=}, {\tt
potpars=} and {\tt potfile=}; they signify the identifying name of the
potential, its parameters and an associated filename.  The last two are
optional, since the potential may not need parameters or an associated
file(s). 

For example, the program {\tt potlist} lists the value of the potential
and forces at selected gridpoints:

\small\begin{verbatim}
    % potlist potname=harmonic potpars=0,3,2 x=0:3:1 y=0:3:1 z=0:6:2 dr=0.001
    x y z ax ay az phi phixx phiyy phizz rho dr time
    0 0 0 -0 -0 -0 0 9 4 1 1.11408 0.001 0 
    1 1 2 -9 -4 -2 8.5 9 4 1 1.11408 0.001 0 
    2 2 4 -18 -8 -4 34 9 4 1 1.11408 0.001 0 
    3 3 6 -27 -12 -6 76.5 9 4 1 1.11408 0.001 0 

\end{verbatim}\normalsize

The usage of the colon separated implied do-loop 
in the {\tt x=, y=} and {\tt z=} keywords assumes that
{\it herinp(3NEMO)} has been implemented.
{\tt potlist} will also take first order derivatives of the
force, to test Poissons equation.  This specific potential, {\tt
harmonic}, happens to have 4 parameters, although only 3 were given
in the example. The fourth one will take some default present in
the descriptor. The first parameter of all
potentials should be the pattern speed\footnote{Although if you supply
your own potential you could cheat and bypass
this}. The second through fourth parameters are the $\omega_X$, 
$\omega_Y$ and $\omega_Z$ harmonic coefficients resp., 
where the potential is given as:

$$
    \Phi(x,y,z) =  {1\over 2} \omega_X^2 x^2
                 + {1\over 2} \omega_Y^2 y^2
                 + {1\over 2} \omega_Z^2 z^2
$$

Although NEMO comes supplied with a small number of standard potential
descriptors, it is relatively easy to make your own
ones.  In Section \ref{buildpot}
we will describe how to add your own
potential descriptors. Next we shall present
a few examples from the standard list of available
potentials, the full listing can be found in Appendix~\ref{a:potential}.

\subsection{A few potentials}
\input potential.inc

\subsection{How to build your own potential descriptors}
\mylabel{ss:buildpot}

Although this subject really is one that should
be deferred to Chapter \ref{c:progr}, we will now present
a simple prototype ``definition'' for a potential\index{potential}
descriptor in C and\index{fortran}
Fortran\footnote{FORTRAN is not supported on all architectures}

\small\begin{verbatim}
    void inipotential (int *npar, double *par, char *name);
    void potential (int *ndim, double *pos, double *acc, double *pot, double *time);

    SUBROUTINE INIPOTENTIAL(NPAR, PAR, NAME)
    SUBROUTINE POTENTIAL(NDIM, POS, ACC, POT, TIME)
\end{verbatim}\normalsize

As you can and will see more of, a potential descriptor is in origin
really a C or FORTRAN source code file, that needs two (FORTRAN) subroutines
or (C) functions with the callable names {\tt
inipotential}\index{inipotential}
and {\tt potential}.\index{potential}  Their arguments must conform to
the specification given above.
Because we do want to allow Fortran source code
as well, all arguments are called by reference in C.

Programs which need a potential descriptor will automatically compile
your source code (if needed) and
load the object code into the program for usage.  The 
repository of standard NEMO potential
descriptors (as object files) lives in {\tt \$NEMOOBJ/potential}, and is
automatically searched when the environment variable {\bf POTPATH} is
appropriately
set.  Note that the two subroutines themselves are not called directly by
the user, but by a workhorse routine from the standard NEMO library. 
This hides much of the interface for the programmer. More details on this
technique can be found in Chapter \ref{c:progr} ({\it still to come}).

Below is a fully commented listing of the {\tt harmonic} potential,
as an example of such a potential descriptor
given in the C language. 
It, and other potentials, can be found in source code form
in the directory  {\tt \$NEMO/src/orbit/potential/data}

\footnotesize\begin{verbatim}
/*
 * harmonic.c: procedures for initializing and calculating 
 *             the forces and potential of a harmonic potential
 */
#include <stdinc.h>                 /* formal NEMO include file */
local double omega = 0.0;       /* defined but not used in here */
local double h[3] = {1.0,1.0,1.0};        /* default parameters */
/*----------------------------------------------------------------
 * INIPOTENTIAL: initializes the potential.
 *      input: npar, the number of parameters
 *             par[] an array of npar parameters
 *      If npar=0 defaults are taken (remember to initialize them 
 *      as static (local) variables in this file)
 *----------------------------------------------------------------
 */
void inipotential (int *npar, double *par, string name)
{
    int i;

    if (*npar>0)
        omega = par[0];
    for (i=1; i<(*npar); i++)
       h[i-1] = sqr(par[i]);
    if (*npar > 4) 
        warning("Only 4 parameters used in Harmonic Potential");

    dprintf (1,"INI_POTENTIAL Harmonic Potential\n");
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  wx^2,wy^2,wz^2= %f %f %f\n",h[0],h[1],h[2]);
}
/*----------------------------------------------------------------
 *  POTENTIAL: the worker routine. Determines at any given point 
 *      (x,y,z) the forces and potential. 
 *      Note that this routine is good for 1, 2 as well as 3D
 *----------------------------------------------------------------
 */
void potential (int ndim,double *pos,double *acc,double *pot,double *time)
{
    int    i;
        
    *pot = 0.0;
    for (i=0; i<*ndim; i++) {
        (*pot) += h[i]*sqr(pos[i]);
        acc[i] = -h[i]*pos[i];
    }
    *pot *= 0.5;
}
\end{verbatim}\normalsize

\section{Orbits}
\mylabel{s:orbits}

\input orbit.inc

\section{Exchanging data}
\mylabel{s:exch-data}

The exchange of (binary) data between machines of different architecture
is often a painful process.  For NEMO binary structured files, we have
devised a general portable way to port files between machines, even if
both have different low level file formats ({\it e.g.} SUN OS and the
Cray UNICOS OS).  \index{Convex} \index{SUN} A different solution has
been used by the MIRIAD package, which writes it's data always in the
same (IEEE) format.  \index{IEEE} This means a different layer of
translation routines is needed for certain architectures, notably VMS
and Unicos.  \index{VMS} \index{Unicos} A similar mechanism is expected
to be used in NEMO in some future release. 

In case the other machine has a \index{data, conversion}
\index{conversion, data} totally different file format, it's handy to
have the data in simple ASCII table format.  NEMO also allows import and
export of N-body data through an ASCII format described in {\tt atos}
\index{atos(1)} It can be used directly for multiple-snapshot data, but
example shell script are available to transport data. 

In the case of N-body data there is no standard format to store the
particle information, and we are subject to someone's favorite format. 
We will encounter a few, and show examples how to convert them under
NEMO.  In the case of images, there happens to be an
astronomical\index{FITS, data conversion} standard: FITS\footnote{See:
Wells et al.  (1981), A\&A Suppl.  44, 363.}, we will discuss a few
applications here too.  Tables can also be transferred in an extended
form of the\index{FITS, table files} FITS format\footnote{See: Harten et
al.  (1988), A\&A Suppl.  73, 365.}, although here the ASCII format may
do equally well.  Even N-body snapshots can be written in FITS format,
for an example see the toy program {\tt snapfits} which uses the now
deprecated FITS Random Group Format.

\subsection{NEMO data files in general}
\mylabel{ss:data}

Here is a neat trick to exchange NEMO data files between systems of
different binary file format.  On machine 1 the data is saved in (UNIX)
compressed ASCII format: \index{tsf(1)} \index{compress}


\small\begin{verbatim}
    m1% tsf r.dat maxprec=t allline=t | compress > r.data.Z
\end{verbatim}\normalsize

The data can then be transported to machine 2 (in binary mode
of course if data was compressed), 
and saved in the local binary structured file format:
\index{zcat} \index{rsf(1)}


\small\begin{verbatim}
    m2% zcat r.data.Z | rsf in=- out=r.dat
\end{verbatim}\normalsize

It turns out that for most data files the compressed ASCII file in full
precision is about as large as the original binary file.  The example
above also shows that, by using pipes, machine 1 and 2 never need to
store the full ASCII version of the file, which will in general be about
4 times as large as the binary file(s).  Note again that a dash filename
\index{-,filename}
is interpreted as standard input/output in the NEMO environment (see
also {\it stropen(3NEMO)}), but one should be warned here that some
older versions of structured files could not be used in 
pipes.\index{pipes}\index{structured files, pipes}

For machines which support I/O redirection in the ftp program, an even
more efficient solution is possible by redirecting
the (compressed) data from the
other machine into local binary structured format:
\index{ftp}


\small\begin{verbatim}
    m2% ftp m1
    ftp> binary
    ftp> get r.data.Z "| zcat | rsf - r.dat"
\end{verbatim}\normalsize

The compressed ASCII data never needs to be stored on the local
disk directly. The data is uncompressed and passed to {\tt rsf}
through a pipe.

\subsection{Snapshot Data}

To import a snapshot into NEMO format one can
use {\tt atos} or write the data in this ASCII (also
referred to as the "205") format.\index{atos(1)}\index{205 format, atos}
In particular a snapshot which is already in
table format with masses, positions and velocities in columns 
1,2-4,5-7, can be converted to snapshot format using a simple
C-shell. For example,

\small\begin{verbatim}
    % table_to_snapshot tab_file snap_file
\end{verbatim}\normalsize

with the following simplified version of the {\tt table\_to\_snapshot}
C-shell script (without any bells and whistles) using {\tt
awk}\index{awk} ({\tt tabmath}\index{tabmath} could have been used
similarly): \index{table, to snapshot}


\small\begin{verbatim}
#! /bin/csh -f
#       table_to_snapshot: demo version
set infile=$1           # input table (m,x,y,z,vx,vy,vz)
set outfile=$2          # output snapshot
set tmpfile=tmp$$       # a temporary scratch name

awk "END {print NR}"   $infile   > $tmpfile
echo "3"                        >> $tmpfile
echo "0.0"                      >> $tmpfile
awk '{print $1}'       $infile  >> $tmpfile
awk '{print $2,$3,$4}' $infile  >> $tmpfile
awk '{print $5,$6,$7}' $infile  >> $tmpfile

atos $tmpfile $outfile
rm $tmpfile
\end{verbatim}\normalsize

The full version of this script can be found in in {\tt \$NEMO/csh}.

\subsection{Image Data}

For images the situation is a little better because there exists a
standard in the astronomical community: the FITS format (see also {\it
fits(5NEMO)}).  {\tt ccdfits} convert a NEMO image to standard FITS
diskfile, which can be read and manipulated by various image processing
packages.  The reverse program, {\tt fitsccd}, is also available, and
can convert most FITS images into NEMO's {\it image(5NEMO)} format. 
When importing FITS images into NEMO, always be concerned with the
units, since NEMO insists that the origin be at (0,0,0). 

\subsubsection{AIPS}
In {\bf AIPS} the following can be done:
\index{AIPS}

Suppose your fits file is stored in a directory, which we will
call {\bf \$dir} ({\it e.g.} {\tt set dir=/usr/nemo/fits}) and the filename
is FITS.DAT (most filenames MUST be in upper case in AIPS), then the AIPS
task {\tt IMLOD} can be used to read the fits file (or from tape, see below):
\index{AIPS,IMLOD} \index{IMLOD,AIPS}

\small\begin{verbatim}
    % setenv XX $dir        # make sure this is set for AIPS
    % aips                  # login/start up AIPS
    ...                     # (some more login stuff here)
    > TASK 'IMLOD'          # set task name
    > INFILE 'XX:FITS.DAT'  # set up input FITS filename
    > INNAME 'TEST'         # set some name for output AIPS file
    > GO                    # run the task
    > MCAT                  # check if file TEST there 
    > GETN ...              # get map number for file TEST
    > TVLOD                 # load it on tv
    > TVFIDDLE              # change contrast
    > EXIT                  # quit AIPS
    % 
\end{verbatim}\normalsize

Using the tape-interface is a bit more cumbersome: dump the 
FITS.DAT file to a tape,
using {\it dd(1)} \index{dd, reading FITS files}
with a block\index{FITS,AIPS} size of 2880 
bytes\footnote{Newer versions of AIPS now allow you to use a blocking
factor which writes blocks in multiples of 2880 bytes; e.g. a blocking factor
of 10 needs block size 28800 bytes.}, and have
{\tt IMLOD} read the data from tape. This in case the disk interface
will not work. Even on DEC-VMS machines the
direct disk fits file may be used (this has done been successfully in
the GIPSY package - see {\it ccdfits(1NEMO))}. \index{FITS,GIPSY}
\index{FITS, ccdfits} \index{ccdfits(1)} \index{GIPSY, FITS}

The reverse process can also be used to write AIPS files to disk in FITS
format using the task {\tt FITTP}, as show in the following example:
\index{AIPS, FITTP} \index{FITTP, AIPS}

\small\begin{verbatim}
    % setenv XX $dir        # make sure this is set for AIPS
    % aips                  # login/start up AIPS
    ...                     # (some more stuff here)
    > MCAT                  # check directory
    > TASK 'FITTP'          # set task name
    > INNAME 'TEST'         # set name of input AIPS file
    > OUTFILE 'XX:FITS.DAT' # filename for output FITS file
    > GO                    # run the task
    > EXIT                  # quit AIPS
    % 
\end{verbatim}\normalsize

The file will then be in {\tt \$dir/FITS.DAT}, make sure that
the program {\tt aips} has write permission in that directory.
Again, if the disk interface does not work, the file has to be dumped
to tape, and read to disk using {\it dd(1)}. 
For an example see {\it ccdfits(1NEMO)}.

\subsubsection{IRAF}

IRAF \index{IRAF, startup}
is normally started up by issuing the {\tt cl} command
(you may need a {\tt login.cl} startup file in your current or home directory).
Converting an existing file in IRAF format into a FITS is very simple,
as is illustrated in the following example:\index{wfits, IRAF}

\small\begin{verbatim}
    % cl                            # startup IRAF
    cl> dataio                      # go into the dataio area
    da> wfits iraf_file fits_file   # and convert it
    da> logout                      # quit IRAF
    %
\end{verbatim}\normalsize

The parameter {\tt bitpix} may have to be set to 16 or 32, if you don't
like it's default.  The complementary IRAF program {\tt rfits} converts
a FITS file into IRAF format. 

\subsubsection{MIRIAD}

MIRIAD can currently not read or write FITS images with {\tt BITPIX=8}.

\subsubsection{IDL}
\subsubsection{MIDAS}



