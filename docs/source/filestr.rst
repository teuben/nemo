.. _filestr:

Filestructure (*)
=================

.. note::
   NEMO stores its persistent data in binary files, which under most circumstances
   can also be used in a Unix pipe by using the ``-`` symbol.

Here we give an overview of the file structure of NEMO's persistent
data stored on disk.
The popular memory (object) models,
and how they interact with persistent data on disk, are discussed in
:ref:`progr`.
Most of the data handled by NEMO is in the form of a
specially designed 
XML-like binary format (well before XML was conceived)
although exceptions like ASCII files/tables will also be discussed. 
Ample examples illustrate creation, manipulation and data transfer.
We also mention a few examples of function descriptors, 
a dataformat that make use of
the native object file format of your operating system
(*a.out(5)* and *dlopen(3)*) that are dynamically loaded during runtime.

Binary Structured Files
-----------------------

.. note::
   There is also a program called **bsf**, which benchmarks a regression value
   of the floating point values in the file.

Most of the data files used by NEMO share a common low level binary file
structure, which can be viewed as a sequence of tagged data items.  Special
symbols are defined to group these items hierarchically into sets.  Data items
are typically scalar values or homogeneous arrays constructed from
elementary C data types, but the programmer can also add more complex
structures, such as C's ``struct`` structure definition, or any user
defined data structure. In this last case tagging by type is not 
possible anymore, and support for a machine independent format
is not guaranteed. Using such constructs is not recommended if the
data needs to be portable accross platforms.

The hierarchical structure of a binary file in this general format can
be viewed in human-readable format at the terminal using a special 
program, ``tsf`` ("*type structured file*").
Its counterpart, ``rsf`` ("*read structured file*"),
converts such human-readable files (in that special ASCII Structured
File format (ASF) into binary structured files (BSF).
In principle it is hence
possible to transfer data files between different types of computers
using ``rsf`` and ``tsf`` (see examples in Section~\ref{s:exch-data}).

Let us start with a small example: With the NEMO
program ``mkplummer`` we first create an
N-body realization of a spherical Plummer model:

.. code-block::

    1% mkplummer i001.dat 1024


Note that we made use of the shortcut that ``out=`` and ``nbody=``
are the first two *program keywords*, and they
were assigned their value by position rather than by associated name.
We can now display the contents of the binary file ``i001.dat`` with
``tsf``:

.. code-block::

    2% tsf i001.dat

    char Headline[33] "set_xrandom: seed used 706921861"
    char History[36] "mkplummer i001.dat 1024 VERSION=2.5"
    set SnapShot                                                            
      set Parameters                                                        
        int Nobj 01750
        double Time 0.00000                                                 
      tes                                                                   
      set Particles                                                         
        int CoordSystem 0201402                                             
        double Mass[1024] 0.00195313 0.00195313 0.00195313
          0.00195313 0.00195313 0.00195313 0.00195313 0.00195313
          0.00195313 0.00195313 0.00195313 0.00195313 0.00195313
          0.00195313 0.00195313 0.00195313 0.00195313 0.00195313
          . . .
        double PhaseSpace[1024][2][3] 4.92932 0.425103 -0.474249
          0.342025 -0.112242 4.60796 -0.00388599 -0.389558 -0.958787
          0.220561 0.213904 3.47561 0.0176012 1.22146 -0.903484
          -0.705422 4.26963 -0.263561 1.04382 -0.199518 -0.480749
          . . .                                                             
      tes                                                                   
    tes                                                                     


This is an example of a data-file from the N-body group, and consists of
a single *snapshot* at time=0.0.  This snapshot,
with 1024 bodies with double precision masses and full 6 dimensional
phase space coordinates, totals 57606 bytes, whereas a straight dump of
only the essential information would have been 57344 bytes, a mere 0.5%
overhead.  The overhead will be larger with small amounts of data,
*e.g.* diagnostics in an N-body simulation, or small N-body snapshots. 

Besides some parameters in the ``Parameters`` set, it consists
of a ``Particles`` set, where (along the type of coordinate system)
all the masses and phase space coordinates of all particles
are defined. Note the convention of integers starting with
a ``0`` in octal representation. This is done for portability
reasons.

A comment about **online** help:
NEMO uses the Unix *man(5)* format
for more detailed online help, 
although the **inline** help (system ``help=`` keyword)
is most of the times sufficient enough
to remind a novice user of the keywords and their meaning.
The ``man`` command is a last resort, if more detailed information
and examples are needed. 

.. code-block::

    3% man tsf


Note that, since the online manual page is a different file from the
source code, information in the manual page can easily get outdated, and
the inline (``help=``) help, although very brief,
is more likely to be up to date since it is generated from the source
code (executable) itself:


.. code-block::

    4% tsf help=h
    in               : input file name  [???]
    maxprec          : print nums with max precision  [false]
    maxline          : max lines per item  [4]
    allline          : print all lines (overrides maxline)  [false]
    indent           : indentation of compound items  [2]
    margin           : righthand margin  [72]
    item             : Select specific item []
    xml              : output data in XML format? (experimental) [f]
    octal            : Force integer output in octal again? [f]
    VERSION          : 29-aug-02 PJT  [3.1]



Pipes
-----

In the UNIX operating system pipes can be very
effectively used to pass information from one process to 
another. One of the well known textbook examples is how one
gets a list of misspelled (or unknown) words from a document:

.. code-block::

    % spell file | sort | uniq | more


NEMO programs can also pass data via UNIX pipes, although with a
slightly different syntax: a dataset that is going to be part of a pipe
(either input or output) has to be designated with  the ``-``
(*dash*) symbol for their filename.
Also, and this is very important, the receiving task
at the other end of the pipe should get data from only one source.
If the task at the sending end of the pipe wants to send binary data over
that pipe, but in addition the same task would also write *normal*
standard
output, the pipe would be corrupted with two incompatible sources of
data. An example of this is the program 
``snapcenter``. The keyword ``report`` must be set to
``false`` instead, which is actually the default now.
So, for example, the output of a previous N-body
integration is re-centered on it's center of mass, and subsequently
rectified and stacked into a single image as follows:

.. code-block::

    % snapcenter r001.dat . report=t  | tabplot - 0 1,2,3
    
    % snapcenter r001.dat - report=f       |\
        snaprect - - 'weight=-phi*phi*phi' |\
        snapgrid - r001.sum stack=t


If the keyword ``report=f`` would not have been set properly,
``snaprect``
would not have been able to process it's convoluted
input. Some other examples
are discussed in Section~\ref{ss:data}.



History of Data Reduction
-------------------------

Most programs
in NEMO will automatically keep track of the history of
their data-files in a self-describing and self-documenting
way. If a program modifies an input file and produces an
output file, it will prepend the
command-line with which it was invoked to its data history.  The
data history is normally located at the beginning of a data file. 
Comments entered using the frequently used program keyword
``headline=`` will also appear in the history section of your data file. 


A utility, ``hisf``
can be used to display the history of a data-file. 
This utility can also be used to create a pure history file (without any
data) by using the optional ``out=`` and ``text=`` keywords.  Of
course ``tsf``
could also be used by scanning its output for the string
``History`` or ``Headline``:

.. code-block::

    5% tsf r001.dat | grep History


which shows that ``tsf``, together with it's counterpart ``rsf`` has
virtually the same functionality as ``hisf``. 


Table format
------------

Many programs are capable of producing standard output in (ASCII)
tabular format.
The output can be gathered into a file using
standard UNIX I/O redirection.  In the example 

.. code-block::

    6% radprof r001.dat tab=true > r001.tab


the file ``r001.tab`` will contain (amongst others) columns with
surface density and radius from the snapshot ``r001.dat``.  These
(ASCII) *table* files can be used by various programs for further
display and analysis.  NEMO also has a few programs for this purpose
available (*e.g.*} ``tabhist`` for analysis and histogram
plotting, ``tablsqfit``
for checking correlations between two columns and
``tabmath`` for general table handling.
The manual 
pages of the relevant NEMO programs should inform you how to get nice
tabular output, but sometimes it is also necessary to write a shell/awk
script or parser to do the job.

A usefull (open source domain) program *redir(1NEMO)*
has been included in NEMO

.. code-block::

    7% redir -e debug.out tsf r001.dat debug=2


would run the ``tsf`` command, but redirecting the
*stderr* standard error output to a file ``stderr.out``. There are
ways in the C-shell to do the same thing, but they are
clumsy and hard to remember. In the bourne/bash
shell this is accomplished much easier:

.. code-block::

    7$ tsf r001.dat debug=2  2>debug.out

One last word of caution regarding tables: tables can also be used
very effectively in pipes, for example take the first example,
and pipe the output into ``tabplot`` to get a quick look 
at the profile:

.. code-block::

    8% snapprint r001.dat r | tabhist - 


If the snapshot contains more than 10,000 points, ``tabhist`` cannot
read the remainer of the file, since the default maximum number
of libes for reading from pipes
is set by a keyword ``nmax=10000``. To properly read all lines, you
have to know (or estimate) the number of lines. In 
the other case where the input is a regular file, table programs
are always able to find the correct amount to allocate for their
internal buffers by scanning over the file once. For very large tables
this does introduce a little extra overhead.

Dynamically Loadable Functions
------------------------------

A very peculiar data file format encountered in NEMO is that of the 
function descriptors. They present themselves to the user through
one or more keywords, and in reality point to a compiled
piece of code that will get loaded by NEMO (using *loadobj(3NEMO)*).
We currently have 4 of these in NEMO:



Potential Descriptors
~~~~~~~~~~~~~~~~~~~~~

The potential descriptor is used in orbit
calculations and a few N-body programs.  These are actually binary
object files (hence extremely system dependent!!), and 
used by the dynamic object loader
during runtime. Potentials are 
supplied to NEMO programs as an input variable (*i.e.* a set of 
keywords, normally called ``potname=``, ``potpars=`` and ``potfile=``.
For this, a mechanism is needed to dynamically load 
the code which calculates the potential. This is done by a
dynamic object loader that comes with NEMO. 
If a program needs a potential, and it is present in the
default repository (``$POTPATH`` or {``$NEMOOBJ/potential``), it is
directly loaded into memory by this dynamic object loader. 
If only a source file is present,
*e.g.* in the current directory, it is compiled on the fly 
and then loaded.  The source code can be written
in C or FORTRAN.  Rules and more information
can be found in *potential(3NEMO)* and *potential(5NEMO)*
The program *potlist(1NEMO)* 
can be used to test potential descriptors. 

Bodytrans Functions
~~~~~~~~~~~~~~~~~~~

Another family of object files used by the dynamic
object loader are the *bodytrans(5NEMO)* functions. These were
actually the first one of this kind introduced in NEMO.
They are functions generated from expressions containing body-variables
(mass, position, potential, time, ordinal number etc.).  They frequently occur
in programs where it is desirable to have an arbitrary
expression of body variables
*e.g.*  plotting and printing programs, sorting program etc.
Expressions which are not in the standard repository (currently 
``$BTRPATH`` or ``$NEMOOBJ/bodytrans``) will 
be generated on the fly and saved for later use. 
The program *bodytrans(1NEMO)* is available
to test and save new expressions. Examples are given in 
Section~\ref{s-dispanal}, a table of the 
precompiled ones are in Table~\ref{t:bodytrans}.


Nonlinear Least Squares Fitting Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The program *tabnllsqfit(1NEMO)* can fit (linear or non-linear, depending
on the parameters) a function to a set of datapoints from an ASCII table.
The keyword ``fit=`` describes the model (*e.g.* a line, plane, gaussian, circle,
etc.), of which a few common ones have been pre-compiled with the program.
In that sense this is different from the previous two function descriptors,
which always get loaded from a directory with precompiled object files.
The keyword ``load=`` can be used to feed a user defined function to
this program. The manual page has a lot more details.

Rotation Curves Fitting Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Very similar to the Nonlinear Least Squares Fitting Functions are the
Rotation Curves Fitting Functions, except they are peculiar to the
1- and 2-dimensional rotation curves one find in galaxies as the 
result of a projected circular streaming model. The program
*rotcurshape(1NEMO)* is the only program that uses these functions, the
manual page has a lot more details.
