% bodies.inc include file for latex - nemo.tex eventually
%% at some point	- Peter Teuben - 14-oct-90
%	24-mar-91  for latest


\section{Comments on a proposed new body structure}

\myfile{bodies.inc}
\mylabel{s:bodies}

\subsection{Filestruct: Basic Data I/O}
All data I/O is to be handled by {\tt filestruct}. Data is tagged
by name, type and dimension, and can be hierarchically structured
in {\tt set}'s (examples follow).

\begin{itemize}
\item {\tt copy\_item(i, o, tag): } recursively copy an item from input
	to output.
\item {\tt get\_data(i,tag,typ,dat,dimN,...,dim1,0): } read data
\item {\tt get\_data\_coerced(i,tag,typ,dat,dimN,...,dim1,0): } read data, 
	but convert to local floating point (double/float)
\item {\tt get\_data\_sub: } worker for above ({\it could be secret})
\item {\tt get\_dims(i,tag): } enquire dimensions of tag
\item {\tt get\_dlen(i,tag): } enquire data length of tag
\item {\tt get\_type(i,tag): } enquire type of tag
\item {\tt get\_set(i,tag): } acces named set in input
\item {\tt get\_string(i,tag): } read a string
% \item {\tt get\_tag(i): } get next tag, if there is one (***)
\item {\tt get\_tag\_ok(i,tag): } check if subsequent set/data has this tag
\item {\tt get\_tes(i, tag): } end acces named set in input
\item {\tt list\_tags(i): } list all tags or which get\_tag\_ok is TRUE.
\item {\tt put\_data(o,tag,typ,dat,dimN,...,dim1,0): } write data
\item {\tt put\_data\_sub: } worker for above ({\it could be secret})
\item {\tt put\_set(o, tag): }	begin a named set in output
\item {\tt put\_string(o, tag, dat): } write string 
\item {\tt put\_tes(o, tag): } begin a named set in output
\item {\tt skip\_item(i): } flush pending item
\item {\tt strclose(i/o): } remove stream from strtable
\end{itemize}

Including {\tt strclose} here may seem a bit strange, since the
{\tt stropen} function is defined by itself. Since a (secret)
table is maintained within {\tt filestruct}, {\tt strclose} needs
to have access to this table in order to remove allocated items
when  the file is closed.

\subsection{Old format!}

The old (internal) 
format of a body is defined through the {\tt <snapshot/body.h>}
include file, or one of it's friends. \index{snapshot, old}
Here's a printout of the full structure definition, together with the
usuall access macros in comments behind their members.

\small\begin{verbatim}
typedef struct {
   real   bodymass;      /* Mass(b)                               */
   vector bodyphase[2];  /* Phase(b) Pos(b) Vel(b)     barebody.h */
   /*-------------------------------------------------------------*/
   real   bodyphi;       /* Phi(b)                                */
   vector bodyacc;       /* Acc(b)                                */
   real   bodyaux;       /* Aux(b)                                */
   int    bodykey;       /* Key(b)                     body.h     */
   /*-------------------------------------------------------------*/
   char  *bodystory;     /* Story(b)                   mybody.h   */
   ---------------------------------------------------------------*/
} body, *bodyptr;
\end{verbatim}\normalsize

In full double precision this body takes up 104 bytes on a 
4-byte word machine
which doesn't have alignment restrictions; before
the {\tt *bodystory} 
\footnote[1]{The bodystory is an experimental member since sep-90}
was included, it took up 104 bytes on the Sparc (because of alignment)
and 100 bytes on the SUN3's. In single precision 
({\bf -DSINGLEPREC -fsingle} on a sun) it takes up 56 bytes. By taking
less relevant components, such as defined by
{\tt barebody.h}, we can go down to a mere 28 bytes per particle.
In fact, if your name is Jim Peebles, you probably only want single
precision positions: 4x3 = 12 bytes/particle, such that 12Mb of
memory gives you just about the 1 Million particles that make him
happy.

The {\tt *bodystory}, or Story's for short, were introduced
to experiment with something similar to the story-mechanism in Starlab.
It points to a NULL terminated character string, who's interpretation
is completeley open and must be defined by the application programmer.
It is a convenient way (but potentially destructive) to pass information
accross levels, and can be acted upon. If no information is present for
a body, the pointer is  NULL. The convention is that a
{\tt Body} can have multiple stories, separated by newlines. 
{\it This causes a problem in the symmetry between {\tt tsf/rsf}, which
has to be fixed}.

The biggest complaint about the current body is the absence of any
decent handling of unknown data. Missing data can be handled by
the application programmer, but any data components not known
to the body, are discarded. We'll come back to that later.

Any future snapshot I/O should still be able to read old snapshot
data.

\newpage
\subsection{N-body snapshot}
An example of a {\tt tsf}-printout of a snapshot with body could look as
follows: \index{story}

\small\begin{verbatim}
set SnapShot
  set Parameters
    int Nobj 012 
    double Time 0 
  tes
  set Particles
    int CoordSystem 0201402 
    double Mass[10] .100000 .100000 .100000 .100000 .100000 .100000 
      .100000 .100000 .100000 .100000 
    double PhaseSpace[10][2][3] .123714 .263689 -.853972 -.138099 
      0.0132169 -.299212 0.0674330 .637106 -.565022 .188248 .340942 
      -.117123 .569567 -.995788 1.01115 -.485795 -.443185 .138675 
      -1.57804 -.869878 -1.29730 -.185420 -.360114 0.0304421 -.709491 
      . . .
  tes
  set Story
    int Nstories 02
    int IStories[1] 02 05
    char story_2[14] "SN at t=1.821"
    char story_5[16] "Just a Red Giant"
  tes
tes
\end{verbatim}\normalsize

\newpage
\subsection{SPH-body snapshot}
A recently introduced SPH-body has the following structure.
\index{sph}

\small\begin{verbatim}
set SnapShot
  set Parameters
    double Time 2.25000 
    int Nobj 0140000 
    int Ngas 020000 
  tes
  set Particles
    int CoordSystem 0201402 
    double Mass[49152] 3.43320e-06 3.43320e-06 3.43320e-06 3.43320e-06 
      3.43320e-06 3.43320e-06 3.43320e-06 3.43320e-06 3.43320e-06 
      3.43320e-06 3.43320e-06 3.43320e-06 3.43320e-06 3.43320e-06 
      . . .
    double PhaseSpace[49152][2][3] -.248320 -.152840 .385600 .328440 
      -0.0414680 .269840 -.255530 -.152860 .385020 .339990 -.382760 
      .244530 -.251950 -.151300 .385330 .264970 -.196550 .270550 
      . . .
    double Density[8192] 311.690 277.780 302.880 357.140 292.900 
      243.670 298.990 173.450 310.530 238.610 92.2610 165.470 351.560 
      88.6400 340.250 271.420 210.110 316.260 281.510 165.540 290.510 
      . . .
    double Temperature[8192] 8728.70 8755.90 8745.50 8755.80 8753.40 
      8722.60 8723.20 8753.90 8758.70 8792.30 8961.70 8813.40 8756.10 
      8961.60 8746.40 8743.90 8784.70 8752.00 8759.40 8838.70 8749.40 
      . . .
    double SmoothLength[8192] 0.00264160 0.00288380 0.00282510 
      0.00251820 0.00263010 0.00255190 0.00250500 0.00316350 
      0.00276170 0.00284600 0.00387870 0.00306980 0.00259260 
      . . .
    double Potential[49152] -9.54700 -9.49110 -9.54720 -9.58510 
      -9.52920 -9.57280 -9.54690 -9.44330 -9.46030 -9.51760 -9.31220 
      -9.38130 -9.57730 -9.24140 -9.57240 -9.53320 -9.46750 -9.54520 
      . . .
  tes
tes
\end{verbatim}\normalsize
In here the first {\tt Ngas} bodies of the {\tt Nobj} total are
special SPH particles with some obvious extra information per
SPH particle in this case. See {\it snapmstat(1NEMO)} for a useful 
application.

This structure shows that the arrays in the {\tt Particles} set do not
all necessarely have to be of the same size.  Header variables, plus a
convention that the SPH particles come first, make that structure
unique. 

We shall use this concept in the next upgrade. {\tt Nobj} remains
the variable to denote the {\bf total} number of objects,
whereas a new variable,{\tt Nbody}, shall be used to denote the first 
(and true) bodies (again of which the first {\tt Ngas} could be SPH particles),
and the remaining number {\tt Nobj-Nbody}  are quasi particles,
or {\tt node}'s possibly in a tree structure.



\newpage
\subsection{Passing unknown information}

Normally a bitmap is maintained of the information
present in the {\tt Particles} item, and all information is
loaded. The {\tt get\_snap/put\_snap} would be more efficient
if the routines are loaded with a structure containing global
information. I believe this is also an easier way to make the
step up to the new format, to be described in the next section.
I'll contrast two structures, next to each other, one to be
used in linear body I/O (the default), but more efficient I/O
can be done by keeping the items in homogeneous arrays.
The I/O call would look a lot simpler, and is easier to
adapt for upgrades. E.g. we could have:

\small\begin{verbatim}
    sb snapbody;
    ...
    SB_THIS(snapbody) = that;
    get_newsnap(instr,&snapbody);
    time = SB_TIME(snapbody)
    ...
\end{verbatim}\normalsize
with a corresponding structure
\small\begin{verbatim}
<get_snap.h>                    <get_snapshot.h>
______________________________________________________________
typedef struct{                 typedef struct {
    int Nobj;                     int Nobj;
    Body *b;                      real *mass;
                                  real *phase;
                                  real *phi;
                                  real *acc;
                                  real *aux;
                                  int  *key;
                                  char **story;
    real time;                    real time;
    int bits;                     int bits;
    int  NIunknown;               int  NIunknown;   /* ?? */
    int  NRunknown;               int  NRunknown;   /* ?? */
    int  *Iunk;                   int  *Iunk;	    /* ?? */
    real *Runk;                   real *Runk;       /* ?? */
    some  *treestruct;            some  *treestruct;
} sb, *sbptr;                   } ss, *ssptr;
\end{verbatim}\normalsize
Probably the structure also needs information how liberal to do
the I/O, which have to be set before the {\tt get\_newsnap} would
be called.

This structure, as listed above assumes too much about the other 
dimensions of the data, and it hence not of general interest.

A better solution perhaps is to maintain
a structure, which reflects the actual hierarchical data structure 
on disk. The {\tt sb} or {\tt ss} structure contains a pointer
to this structure if there is such data. In {\tt filesecret.h} this
is defined as an {\tt item}:
\small\begin{verbatim}
typedef struct {		% filestruct item
  string itemtyp;     /* type string, listed in filestruct.h */
  int    itemlen;     /* length associated with above type */
  string itemtag;     /* name given this item by application */
  int   *itemdim;     /* int-string of dimensions, or NULL */
  byte  *itemdat;     /* the real goodies, if any, or NULL */
  long   itempos;     /* where the item began in input stream */
} item, *itemptr;
\end{verbatim}\normalsize
This is the most basic data structure known to Nemo, and is
normally not visible to the programmer.

A common {\tt select=} keyword and corresponding routines could be
formalized to describe the snap I/O. The arguments to the keyword
then contain directive what to read/write etc.

\newpage
\subsection{New format?}

With the new ideas that we have just discussed, we can
propose a new format for the {\tt Body} structure. It would still
exist in more or less its original form. 
Access routines, similar to {\tt get\_snap}
and {\tt put\_snap} can still be used to read in a snapshot in
the old way: a simple linear array of bodies.

However, an additional function has to be called subsequently
to impose a tree structure on a newly to be created
array of pointers to the bodies (which are linearly arranged in 
memory at startup time). 

The new structure could be as follows (still a proposal).

\small\begin{verbatim}
typedef struct {
    real   bodymass;
    vector bodyphase[2];
    char  *bodystory;
    char  *extbody;         /* all other stuff */
} body, *bodyptr;
\end{verbatim}\normalsize

The essential difference here is the pointers to the {\tt real} array
{\tt bodyraux} and {\tt int} array {\tt bodyiaux} to
contain an arbitrary amount of real or integer data 
per body. If they are {\tt NULL}, no information is
present. The {\tt Parameters} set (item) must contain the second
dimension of these new data items.
In addition there might be a pointer to all the other unknown stuff.
We could keep the Body more compatible with the olders and, and
keep the {\tt Aux} and {\tt Key} field, whereas the rest is caught 
under the {\tt unkown} umbrella.


When 'new' extended format is read by old programs, the tree-extension
({\tt Nobj-Nbody}) would be lost!
\newpage
For example, we could now have:
\small\begin{verbatim}
set SnapShot
  set Parameters
    int Nobj 012
    int NRaux 01    
    int NIaux 02
    double Time 0 
  tes
  set Particles
    int CoordSystem 0201402 
    double Mass[10] .100000 .100000 .100000 .100000 .100000 .100000 
      .100000 .100000 .100000 .100000 
    double PhaseSpace[10][2][3] .123714 .263689 -.853972 -.138099 
      0.0132169 -.299212 0.0674330 .637106 -.565022 .188248 .340942 
      -.117123 .569567 -.995788 1.01115 -.485795 -.443185 .138675 
      -1.57804 -.869878 -1.29730 -.185420 -.360114 0.0304421 -.709491 
      . . .
    int Raux[10] 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0
    double Iaux[10][2] 01 02 02 03 05 06 04 05 03 02 01 03 05 05
      02 03 04 07 05 06
    double Color[10] .117123 .569567 .995788 1.01115
      .117123 .569567 .995788 1.01115 .485795 .443185 .138675 
      . . .
  tes
  set Story
    int Nstories 02
    int IStories[1] 02 05
    char story_2[14] "SN at t=1.821"
    char story_5[16] "Just a Red Giant"
  tes
tes
\end{verbatim}\normalsize

As you can see, the data are present as hom. arrays, except the 
{\tt aux} arrays have integer variables in the header denoting
the second dimension.

\newpage
Next, to incoorporate the hierarchical data structure, an array
of pointers to the appropriate bodies is maintained. Initially
all bodies are stored as a linear array, if the
{\tt get\_snap()} type operation is to be performed. This is still
achieved by a simple:

\small\begin{verbatim}
    int nbody;      /* number of bodies to be read */
    Body *btab;     /* pointer to the linear array of bodies */
    
    btab = (Body *) malloc( nbody * sizeof(Body) );
\end{verbatim}\normalsize

exactly as in the older version of Nemo.


When the tree structure is being build, new {\it node}'s, defining
the tree structure of the original {\tt Nobj} particles, will
grow by a factor of about $\log_2{N}$. These {\it node}'s are added
to the list of particles when snap I/O is done with the request
of the tree structure to be written too. The full {\tt Body} structure
is available to each {\it node}.

Perhaps at I/O time a flag needs to be set if the tree part is to be
read, and hence initialized. This means all of the tree I/O could
be done hidden from the user.

\newpage
The snapshot could now look something like:

\small\begin{verbatim}
set SnapShot
  set Parameters
    int Nobj 014
    int Nbody 012
    double Time 0 
  tes
  set Particles
    int CoordSystem 0201402 
    double Mass[12] .100000 .100000 .100000 .100000 .100000 .100000 
      .100000 .100000 .100000 .100000 1.000000 0.40000
    double PhaseSpace[12][2][3] .123714 .263689 -.853972 -.138099 
      0.0132169 -.299212 0.0674330 .637106 -.565022 .188248 .340942 
      -.117123 .569567 -.995788 1.01115 -.485795 -.443185 .138675 
      -1.57804 -.869878 -1.29730 -.185420 -.360114 0.0304421 -.709491 
      . . .
  tes
  set Story
    int Nstories 02
    int IStories[1] 02 05 011 012
    char story_2[14] "SN at t=1.821"
    char story_5[17] "Just a Red Giant"
    char story_11[19] "I am just the root"
    char story_11[17] "some crazy node"
  tes
tes
\end{verbatim}\normalsize
with the convention that the first {\tt Nbody} objects are actually
proper N-body's or SPH-Body's, and the remainder the flattened
tree structure. {\it Hmmm, how do you rebuild it..., guess we need
a index pointer array for that...}.
