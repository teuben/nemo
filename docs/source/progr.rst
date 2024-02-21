.. _progr:

Programmers Guide
=================

.. warning::
   There are still some relics of the old manual of features that
   have been deprecated.


In this section an introduction is given how to write programs within
the NEMO environment.  It is based on an original report *NEMO:
Elementary Mechanics Observatory* by Joshua Barnes (1986).

To the application programmer NEMO consists of a set of (C) macro
definitions and object libraries, designed for numerical work in general
and stellar dynamics in particular.  A basic knowledge how to program in
C, and use the local operating system to get the job done, is assumed. 
If you know how to work with *Makefile*'s, even better. 

We start by looking at the *Hello Nemo* program:

.. see also hello.src

.. code-block:: C
   
	#include <nemo.h>                         /* standard (NEMO) definitions */

	string defv[] = {       /* standard keywords and default values and help */
	    "n=10\n                    Number of interations",           /* key1 */
	    "VERSION=1.2\n             25-may-1992 PJT",                 /* key2 */
	    NULL,               /* standard terminator of defv[] vector          */
	};

	string usage = "Example NEMO program 'hello'";        /* usage help text */

	void nemo_main()                   /* standard start of any NEMO program */
	{
	    int n = getiparam("n");                          /* get n            */
    
	    printf("Hello NEMO!\n");                         /* do some work ... */
	    if (n < 0)                                       /* deal with fatal  */
	       error("n=%d is now allowed, need >0",n);      /* errors           */
	}
        


The NEMO Programming Environment
--------------------------------

The modifications necessary to your UNIX environment in order to
access NEMO are extensively described elsewhere.  This not only
applies to a user, but also to the application programmer.

In summary, the essential changes to your environment consist of
one simple additions to your local shell startup file
or you can do it interactively in your shell (be it sh or csh like).

.. code-block:: bash

      % source /opt/nemo/nemo_start.sh
   or
      % source /opt/nemo/nemo_start.csh

where the location of NEMO=/opt/nemo is something that will likely
be different for you.


The NEMO Macro Packages
-----------------------

We will describe a few of the most frequently used macro packages
available to the programmer.  They reside in the form of header include
files in a directory tree starting at ``$NEMOINC``.  Your application
code would need references like:

.. code-block:: C

    #include <nemo.h>              // this will include a few basic core NEMO include files
    #include <snapshot/body.h>
    #include <snapshot/get_snap.c>

Some of the macro packages are merely function prototypes, to
facilitate modern C compilers (we still use the C99 standard),
and  have associated object code in libraries in ``$NEMOLIB`` and programs
need to be linked with the appropriate libraries (normally controlled
via a ``Makefile``).


stdinc.h
~~~~~~~~

The macro package ``stdinc.h`` provides all basic
definitions that ALL of NEMO's code must include as the first include
file.  It also replaces the often used ``stdio.h`` include file in C
programs.  The ``stdinc.h`` include file will provide us with a way to
standardize on future expansions, and make code more
machine/implementation independent (*e.g.* POSIX.1).  In
addition, it defines a more logical standard for C notation.  For
example, the normal C practice of using pointers to character for
pointer to byte, or integer for bool, tends to encourage a degree of
sloppy programming, which can be hard to understand at a later date. 

A few of the basic definitions in this package:


- ``NULL``:
  macro for 0, used to distinguish null characters and null pointers.
  This is often already defined by ``stdio.h``.


- ``bool``:
  typedef for ``short int`` or ``char``,
  used to specify boolean data. See also next item.
  NOTE: the *curses* library also defines ``bool``, and this
  made us change from ``short int`` to the more space
  saving ``char``.

- ``TRUE, FALSE``:
  macros for 1 and 0, respectively, following normal C conventions,
  though a non-zero value can also be used as *true*.

- ``byte``:
  typedef for ``unsigned char``, used to specify byte-sized data.

- ``string``:
  typedef for ``char *``, used to point to strings. Don`t use
  ``string`` for pointers you increment, decrement or explicitly
  follow (using ``*``); such pointers are really ``char *``.


- ``real, realptr``:
  typedef for ``float`` or ``double`` (``float *`` or ``double *``,
  respectively), depending on the use of the ``SINGLEPREC`` flag. 
  The default is ``double``. 

- ``proc, iproc, rproc``: 
  typedefs for pointers to procedures (void functions), integer-valued
  functions and real-valued functions respectively. 
  % This always confusing issue will become clear later on.

- ``local, permanent``:
  macros for ``static``. Use ``local`` when declaring variables or
  functions within a file to be local to that file. They will not appear
  in the symbol table be usable as external symbols. Use 
  ``permanent`` within  a function, to retain their value upon
  subsequent re-entries in that function.

- ``PI, TWO_PI, FOUR_PI, HALF_PI, FRTHRD_PI``:]
  macros for :math:`\pi`, :math:`2\pi`,  :math:`4\pi`,  :math:`4\pi/3`, 
  respectively.

- ``ABS(x), SGN(x)``:
  macros for absolute value and sign of ``x``, irrespective of the 
  type of ``x``. Beware of side effects, it's a macro!

- ``MAX(x,y), MIN(x,y)``:
  macros for the maximum and minimum of ``x,y``, irrespective of the
  type of ``x,y``. Again, beware of side effects.

- ``stream``:
  typedef for ``FILE *``. They are mostly used with the NEMO
  functions ``stropen`` and ``strclose``, which are functionally
  similar to *fopen(3)* and *fclose(3)*, plus some added 
  NEMO quircks like (named) pipes and the dot (.) ``/dev/null`` file.


getparam.h
~~~~~~~~~~

The command line syntax 
is implemented by a small set of functions used by all conforming
NEMO programs.  A few function calls generally suffice to get the values
of the input parameters.  A number of more
complex parsing routines are also available, to be discussed below.

First of all, a NEMO program must define which 
**program keywords** it will 
recognize. For this purpose it must define an array of *string*s
with the names and the default values for the keywords, and optionally,
but STRONGLY recommended,
a one line help string for that keyword:

.. code-block:: C

    #include <nemo.h>       /*   every NEMO module needs this  */

    string defv[] = {       /* definitions of the keywords */
        "in=???\n          Input file (a snapshot)",
        "n=10\n            Number of particles to view",
        "VERSION=1.1\n     14-jul-89 - 200th Bastille Day - PJT",
        NULL,
    };

    string usage = "example program";   /* def. of the usage line */


The *keyword=value* and *help* part of the string must be
separated by a newline symbol (``\n``). If no 
newline is present,
as was the case in NEMO's first release, no help 
string is available.
ZENO uses a slightly different technique, where strings starting with ``;``
are the help string. This example in ZENO would look as follows:

.. code-block:: C

    string defv[] = {      "; example program",
        "in=???",          ";Input file (a snapshot)",
        "n=10",            ";Number of particles to view",
        "VERSION=1.1",     ";14-jul-89 - 200th Bastille Day - PJT",
        NULL,
    };

You can see the first string is actually the same as NEMO's *usage* string.
The current NEMO *getparam* package is able to parse both NEMO and ZENO style
``defv[]`` initializers.
   
  
The ``help=h`` command line option displays the *help* part of string during
execution of the program for quick inline reference.
The *usage* part defines a string that is used as a one
line reminder what the program does. It's used by the various
invocations of the user interface.
  
The first thing a NEMO program does, is comparing the command
line arguments of the program (commonly called
``string argv[]`` in a C program) 
with this default vector of *keyword=value*
strings (``string defv[]``), and replace
appropriate reset values for later retrieval. This is done by calling
``initparam`` as the first step in your MAIN program:

.. code-block:: C

    main (int argc, string argv[])
    {
        initparam(argv,defv);
        . . .


It also checks if keywords which do not have a default
value (*i.e.* were given ``???``) 
have really been given a proper 
value on the command line, if keywords are not specified twice, 
enters values of the system keywords etc.

There is a better alternative to define the main part of a NEMO program:
by renaming the main entry point ``main()`` to 
``nemo_main()``, without any arguments, and
calling the array of strings with default *key=val*'s
``string defv[]``, the linker will automatically include
the proper startup code (``initparam(argv,defv)``), the worker routine
``nemo_main()`` itself, and the stop code (``finiparam()``). 
The above section of code would then be replaced by a mere:

.. code-block:: C

    void nemo_main()
    {
        . . .

This has the obvious advantage that various NEMO related 
administrative details
are now hidden from the application programmers, and occur automatically.
Remember that standard ``main()`` already shields the application
programmer from a number of tedious setups (e.g. *stdio* etc.).
Within NEMO we have taken this one step further. A recent example that
was added to ``nemo_main`` is the management  of the number of processors
in an OpenMP enhanced computing mode.


Once the user interface has been initialized, keyword values may be obtained at
any point during execution of the program by calling ``getparam()``,
which returns a string

.. note::
   ANSI rules say you can't write to this location in memory if
   they are direct references ``string defv[]``;
   this is something that may well be fixed in a future release.

.. code-block:: C

    if (streq(getparam("n"),"0")
        printf(" You really mean zero or octal?\n");


There is a whole family of ``getXparam()`` functions which 
parse the
string in a value of one of the basic C types ``int, long, bool,``
and ``real``. It returns that value in that type: 

.. code-block:: C

    #include <getparam.h>     //   included by <nemo.h>
    . . .
    int nbody;
    . . .
    if ( (nbody = getiparam("n")) <= 0) {
        printf("Cannot handle %d particles\n",nbody);
        exit(0);
    }


Finally, there is a macro called ``getargv0()``, which returns
the name of the calling program, mostly used for identification:

.. code-block:: C

    if (getbparam("quit"))
        error("%s: early quit", getargv0());

This is very useful in library routines, who normally would not be
able to know who called them. Actually, NEMO's ``error`` function
already uses this technique, since it cannot know the name of the
program by whom it was called.
The ``error`` function prints a message, and exits the 
program.

More detailed information can also be found in the appropriate manual
page: *getparam(3NEMO)* and *error(3NEMO)*.


Advanced User Interface and String Parsing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we describe *setparam* to add some interactive capabilities in
a standard way to NEMO.  Values of keywords should only
be accessed and modified this way.  Since
keywords are initialized/stored within the source code, most compilers
will store their values in a read-only part of data area in the
executable image.  Editing them may cause unpredictable behavior. 

If a keyword string contains an array of items of the same type, 
one can use either ``nemoinpX`` or 
``getXrange``, 
depending if you know how many items to expect in the string.
The ``getXrange`` routines will allocate a new array which will contain
the items of the parsed string. If you do already have a
declared array, and know that all items will fit in there, the
``nemoinpX`` routines will suffice.

An example of usage:

.. code-block:: C

    double *x = NULL;
    double y[NYMAX];
    int nxret, nyret;
    int nxmax=0;

    nyret = nemoinpd(getparam("y"), y, NYMAX);

    nxret = getdrange(getparam("x"), &x, &nxmax);


   

In the first call the number of elements to be parsed
from an input keyword ``y=`` is limited to ``NYMAX``, and is useful
when the number of elements is expected to be small or more
or less known. The actual
number of elements returned in the array ``y[]`` is ``nyret``.

When the number of elements to be parsed is not known at all, or one needs
complete freedom, the dynamic allocation feature of *getdrange* can
be used. The pointer ``x`` is initialized to ``NULL``, as well as
the item counter ``nxmax``.
After
calling ``getdrange``, ``x`` will point to an array of length
``nxmax``, in which the first ``nxret`` element contain the parsed
values of the input keyword ``x=``. Proper re-allocation will be done
when a larger space is need on subsequent calls.

Both routines return negative error return codes, see *nemoinp(3NEMO)*.

More complex parsing is also done by calling ``burststring`` first
to break a string in pieces, followed by a variety of functions.

Alternatives to nemo_main
~~~~~~~~~~~~~~~~~~~~~~~~~

It is not required for your program to define with ``nemo_main()``. There
are cases where the user needs more control.  An example of this is the
``falcON`` N-body code in ``$NEMO/use/dehnen/falcON``.
A header file (see e.g. ``inc/main.h``) now defines main, instead of
through the NEMO library:

.. code-block:: C

    // in main.h:

    extern string defv[];
    extern string usage;

    namespace nbdy {
       char version [200];                             // to hold version info
       extern void main();                             // to be defined by user
    };

    int main(int argc, char *argv[])                  // ::main()
    {
       snprintf(nbdy::version,200,"VERSION=" /*..*/ ); // write version info
       initparam(argv,defv);                           // start NEMO
       nbdy::main();                                   // call nbdy::main()
       finiparam();                                    // finish NEMO
    }


and the application includes this header file, and defines the keyword list
in the usual way :

.. code-block:: C

    // in application.cc

    #include <main.h>

    string defv[] = { /*...*/, nbdy::version, NULL }; // use version info
    string usage  = /*...*/ ;

    void nbdy::main() { /*...*/ }                     // nbdy::main()



filestruct.h
~~~~~~~~~~~~

The *filestruct* package provides a direct and consistent way of passing
data between NEMO programs, much as *getparam* provides a way of passing
(command line) arguments to programs.
For reasons of economy and
accuracy, much of the data manipulated by NEMO is stored on disk in
binary form.  Normally, data stored this way is completely
unintelligible, except to specialized programs which create and access
it.  Furthermore, this approach to data handling tends to be very
brittle: a trivial addition or alteration to the data stored in such a
file can force the tedious and error-prone revision of many programs. 
To get around these problems and provide an explicit, flexible, and
structured method of storing binary data,  we developed a collection of
general purpose routines to access binary data files. 

From the programmers point of view, a structured binary file 
is a stream of tagged data objects. In XML parlor, you could
view this as binary XML.
These objects come in two classes.  An *item*
is a single instance or a regular array of one of the following C
primitive types: ``char, short, int, long, float`` or ``double``.  A
*set*
is an unordered sequence of *items* and *sets*. This definition is
recursive, so fully hierarchical file structures are allowed, and indeed
encouraged.  Every set or item has a name *tag* 
associated with it,
used to label the contents of a file and to retrieve objects from a set. 
Data items have a *type*
and array *dimension* attributed
associated with them as well. This of course means that there is a little
overhead, which may become too large if many small amounts of data
are to be handled. For example, a snapshot  with 128 bodies
(created by ``mkplummer``) with double
precision masses and full 6 dimensional phase space coordinates
totals 7425 bytes, whereas a straight dump of only the essential
information would be 7168 bytes, a mere 3.5% overhead. After an
integration, with 9 full snapshots stored and 65 snapshots with only
diagnostics output, the overhead is much larger: 98944 bytes of
data, of which only 64512 bytes are masses and phase space coordinates:
the overhead is 53% (of which 29% though are the diagnostics output,
such conservation of energy and angular momentum, cputime, center
of mass, etc.).


The filestruct package uses ordinary *stdio(3)* streams to access
input and output files; 
hence the first step in using filestruct is to
open the file streams.  For this job we use the NEMO library routine
``stropen()``,
which itself is not part of filestruct.
``stropen(name,mode)`` is much like ``fopen()`` of *stdio*, but
slightly more clever; it will not open an existing file for output,
unless the *mode* string is ``"w!"``. (append???)
An additional oddity to
``stropen`` is that it treats the dash filename ``"-"``
as standard in/output,
and ``"s"`` as a scratch file.
Since *stdio* normally flushes all buffers on exit, it is
often not necessary to explicitly close open streams, but if you do so,
use the matching routine ``strclose()``. This
also frees up the table
entries on temporary memory used by the *filestruct* package. As in most
applications/operating systems a task can have a limited set of
open files associated with it. Scratch files 
are automatically deleted from disk
when they are closed.

Having opened the required streams, it is quite simple to use the basic
data I/O routines.  For example, suppose the following declarations have
been made:

.. code-block:: C

    #include <stdinc.h>
    #include <filestruct.h>

    stream instr, outstr;
    int    nbody;
    string headline;

    #define MAXNBODY 100
    real    mass[MAXNBODY];


(note the use of the ``stdinc.h`` conventions).  And now suppose that, 
after some computation, results have been stored in the first ``nbody``
components of the ``mass`` array, and a descriptive message has been
placed in ``headline``.  The following piece of code will write the
data to a structured file:

.. code-block:: C

    outstr = stropen("mass.dat", "w");

    put_data(outstr, "Nbody", IntType, &nbody, 0);
    put_data(outstr, "Mass", RealType, mass, nbody, 0);
    put_string(outstr, "Headline", headline);

    strclose(outstr);    

Data (the 4th argument in ``put_data``,
is always passed by address, even if one element is written. This
not only holds for reading, but also for writing, as is apparent from
the above example.
Note that no error checking is needed when the file is opened for
writing. If the file ``mass.dat``
would already have existed, ``error()`` would 
have been called inside ``stropen()`` and aborted the
program. Names of tags are arbitrary, but we encourage you to use
descriptive names, although an arbitrary maximum of 64 is enforced
by chopping any incoming string.

The resulting contents of ``mass.dat`` can be viewed with the
``tsf`` utility:

.. code-block:: bash

    % tsf mass.dat
    int Nbody 8
    double Mass[8] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    char Headline[20] "All masses are equal"


It is now straightforward to read data from this file:

.. code-block:: C

    instr = stropen("mass.dat", "r");

    get_data(instr, "Nbody", IntType, &nbody, 0);
    get_data(instr, "Mass", RealType, mass, nbody, 0);
    headline = get_string(instr, "Headline");

    strclose(instr);

Note that we read the data in the same order as they were written.

During input, the filestruct routines normally perform strict type-checking;
the tag, type and dimension supplied to ``get_data()`` must match
the attributes of the data item, written previously,
exactly. Such strict checking helps
prevent many common errors in using binary data. Alternatively, you can use
``get_data_coerced()``, which is called exactly like ``get_data()``, but
interconverts ``float`` and ``double`` 
values.

To provide more flexibility in programming I/O, a series of related items
may be hierarchically wrapped into a set:

.. code-block:: C

    outstr = stropen("mass.dat", "w");

    put_set(outstr, "MySnapShot");
       put_data(outstr, "Nbody", IntType, &nbody, 0);
       put_data(outstr, "Mass", RealType, mass, nbody, 0);
       put_string(outstr, "Headline", headline);
    put_tes(outstr, "MySnapShot");

    strclose(outstr);    


Note that each ``put_set()``
must be matched by an equivalent  ``put_tes()``.
For input, corresponding routines ``get_set()`` and ``get_tes()``
are used. These also introduce a significant additional functionality:
between a ``get_set()``
and ``get_tes()``, the data items
of the set may 
be read in any order (tagnames must now be unique within an item,
as in a C ``struct``,
or not even read at all). For 
example, the following is also a legal way to access 
the ``MySnapShot``:

.. code-block:: C

    instr = stropen("mass.dat", "r");

    if (!get_tag_ok(instr,"MySnapShot"))
        error("File mass.dat is not a MySnapShot");

    get_set(instr,"MySnapShot");
    headline = get_string(instr, "Headline");
    get_tes(instr,"MySnapShot");

    strclose(instr);    

This method of *filtering* a data input stream clearly opens up many
ways of developing general-purpose programs. Also note that the 
``bool`` routine ``get_tag_ok()`` can be used to control the
flow of the program, as ``get_set()`` would call ``error()`` when
the wrong tag-name would be encountered, and abort the program.

The UNIX program ``cat`` can also be used to catenate
multiple binary data-sets into one,  *i.e.*

.. code-block:: bash

    % cat mass1.dat mass2.dat mass3.dat > mass.dat

The ``get_tag_ok`` routine can be used to handle such multi-set data
files.  The following example shows how loop through such a combined
data-file. 

.. code-block:: C

    instr = stropen("mass.dat", "r");

    while (get_tag_ok(instr, "MySnapShot") {
       get_set(instr, "MySnapShot");
       get_data(instr, "Nbody", IntType, &nbody, 0);
       if (nbody > MAXNBODY) {
            warning("Skipping data with too many (%d) items",nbody);
            get_tes(instr,"MySnapShot");
            continue;
       }
       get_data(instr, "Mass", RealType, mass, nbody, 0);
       headline = get_string(instr, "Headline");
       get_tes(instr,"MySnapShot");
       /*   process data   */
    }

    strclose(instr);    


The loop is terminated at either end-of-file, or if the next object in
``instr`` is not a ``MySnapShot``. 

It is easy to the skip for an item if you know if it is there:

.. code-block::

    while(get_tag_ok(instr,"MySnapShot"))     /*  ??????? */
        skip_item(instr,"MySnapShot");

The routine ``skip_item()`` is only effective, or for that matter
required, when doing input at the top level, *i.e.* not between
a ``get_set()`` and matching ``get_tes()``, since I/O at deeper
levels is random w.r.t. items and sets. In other words,
at the top level I/O is sequential, at lower levels random.

A relative new feature in data access is the ability to
do random and blocked access to the data. Instead of using
a single call to ``get_data`` and ``put_data``, the
access can be sequentially blocked using 
``get_data_blocked`` and ``put_data_blocked``, provided
it is wrapped using ``get_data_set`` and ``get_data_tes``,
for example:

.. code-block:: C

   get_data_set    (instr, "Mass", RealType, nbody, 0);
   real *mass = (real *) allocate((nbody/2)*sizeof(real));
   get_data_blocked(instr, "Mass", mass, nbody/2);
   get_data_blocked(instr, "Mass", mass, nbody/2);
   get_data_tes    (instr, "Mass");

would read in the ``Mass`` data in two pieces into a smaller sized ``mass``
array. A similar mode exists to randomly access data with an item. A current
limitation of this mode is that such access is only allowed on one item at a
time. In this mode an item must be closed before the next one can be opened
in such a mode.

vectmath.h
~~~~~~~~~~

The ``vectmath.h`` macro package
provides a set of macros to handle
some elementary operations on two, three or general **N** dimensional
vectors and matrices.  The dimension $N$ can be picked by providing the
package with a value for the preprocessor macro **NDIM**.  If this is
not supplied, the presence of macros **TWODIM** and **THREEDIM** will
be checked, in which case **NDIM** is set to 2 or 3 respectively.
The default of **NDIM**
when all of the above are absent,
is 3. Of course, the macro **NDIM** must be provided before
``vectmath.h`` is included to have any effect.
Resetting the value of **NDIM** after that,
if your compiler would allow it anyhow without an explicit ``#undef``,
may produce unpredictable results.

There are also a few of the macro's which can be used as a regular
C function, returning a real value, *e.g.*  ``absv()`` for the
length of a vector.

Operations such as **SETV** (copying a vector) are properly defined for
every dimension, but **CROSSVP** (a vector cross product) has a
different meaning in 2 and 3 dimensions, and is
absent in higher dimensions.

It should be noted that the matrices used here are true C
matrices, a pointer to an array of pointers (to 1D arrays), unlike
FORTRAN arrays, which only occupy a solid 2D block of memory. C
arrays take slightly more memory. For an example how to make C arrays and
FORTRAN arrays work closely together see e.g. *Numerical
Recipes in C* by *Press et al.* (MIT Press, 1988).
and NEMO's *mdarray(3NEMO)* package.

In the following example a 4 dimensional vector is cleared:

.. code-block:: C

    #define NDIM 4
    #include <vectmath.h>

    nemo_main()
    {
        vector a;           /* same as:   double a[4]  */

        CLRV(a);
    }



snapshots: get_snap.c and put_snap.c
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These routines exemplify an attempt to provide truly generic I/O of
N-body data.  They read and write structured binary data files
conforming to the overall form seen in earlier sections.
Internally they
operate on ``Body`` structures; 
A ``Body`` has components accessed
by macros such as ``Mass`` for the mass, ``Pos`` and ``Vel`` for
the position and velocity vectors,  etc.  Since ``get_snap.c``
and ``put_snap.c`` use only these macros to declare and access
bodies, they can be used with any suitable ``Body`` structure.  They
are thus provided as C source code to be included in the compilation of
a program. Definitions concerning Body's and snapshots are obtained 
by including the files ``snapshot/body.h`` and ``snaphot/snapshot.h``.

A program which should handle a large number of particles, may decide to
include a more simple ``Body`` structure, as is e.g. provided by
the ``snapshot/barebody.h`` macro file.  This body only includes the
masses and phase space coordinates, which would only occupy 28 bytes 
per particle (in single precision), as opposed to the 100/104
bytes per particle
for a double precision ``Body`` from the standard ``snapshot/body.h`` 
macro file.  This last one contains ``Mass, PhaseSpace, Phi, Acc, Aux``
and ``Key``.

In the example listed under Table~\ref{src:snapfirst}
the first snapshot of an input file is
copied to an output file. 

..   $NEMO/src/tutor/snap/snapfirst.c}

.. include:: snapfirst.src
   :literal:


Notice that the first argument of ``stropen()``, the filename, is
directly obtained from the user interface. The input file is
opened for reading, and the output file for writing.
Some history (see below)
is obtained from the input file (would we not have done
this, and the input file would have contained history, 
a subsequent ``get_snap()`` call would have failed to find
the snapshot), and the first snapshot is read into an array
of bodies, pointed to by ``btab``. Then the output file
has the old history written to it (although any command line
arguments were added to that), followed by that first snapshot.
Both files are formally closed before the program then returns.

history.h
~~~~~~~~~

When performing high-level data I/O, as is offered by a package such as
``get_snap.c`` and ``put_snap.c``, there is an automated way to
keep track of data history. 

When a NEMO program is invoked, the program name and command line
arguments are saved by the ``initparam()`` in a special history
database.  Most NEMO programs will write such history items to their
data-file(s) before the actual data.  Whenever a data-file is then
opened for reading, the programmer should first read these data-history
items. Conversely, when writing data, the history should
be written first. In case of the ``get/put_snap`` package:

.. code-block:: C

    get_history(instr);
    get_snap(instr, &btab, &nbody, &time, &bits);
       /*     process data     */
    put_history(outstr);
    put_snap(outstr,&btab, &nbody, &time, &bits);

Private comments should be added with the ``app_history()``
..   \footnote{The old name, {\tt add\_history} was already used by the GNU {\it readline} library}
When a series of snapshot is to be processed, it is recommended that
the program should only be output
the history once, before the first output of the snapshot, as in the
following example:

.. code-block:: C

    get_history(instr);
    put_history(outstr);
    for(;;) {
        get_history(instr);     /* defensive but in-active */
        get_snap(instr, &btab, &nbody, &time, &bits);
            /*     process data  and decide when done */
        put_snap(outstr,&btab, &nbody, &time, &bits);
    }

Note that the second call to ``get_history()``, within the for-loop, 
is really in-active. If there happen to be history items sandwiched
between snapshots, they will be read and added to the history stack,
but not written to the output file, since ``put_history()`` was only 
called before the for-loop. It is only a defensive call:
``get_snap()`` would fail since it expects only pure
``SnapShot`` sets (in effect, it calls ``get_set(instr,"SnapShot")``
first, and would call ``error()`` if no snapshot encountered).

Building NEMO programs
----------------------


Besides writing the actual code for a program, an application programmer
has to take care of a few more items before the software can be added 
and formally be accepted to NEMO.  This concerns 
writing
the documentation and possibly a Makefile, the former one preferably 
in the form of standard UNIX manual pages.
We have templates for both Makefile's and 
manual pages. Both these are discussed in detail in the next subsections. 

Because NEMO is a development package within which
a multitude of people are donating software and
libraries, linking a program can become cumbersome. In the most
simple case however (no graphics or mathematical libraries needed),
only the main NEMO library is needed, and the
following command should suffice to produce an executable:

.. code-block:: bash

    % gcc -g -o snapprint snapprint $NEMOLIB/libnemo.a -lm



Each user is given a subdirectory in ``$NEMO/usr``, under which
code may be donated which can be compiled into the running version
of NEMO. Stable code, which has been sufficiently tested and
verified, can be placed in one of the appropriate ``$NEMO/src`` 
directories. For proper inclusion of user contributed software
a few rules in the ``Makefile`` have to be adhered to.

The ``mknemo``
script should handle compilation and
installation of most of the standard NEMO cases. Some programs, like 
the N-body integrators, are almost like complicated packages themselves,
and require their own Makefile or install script. For most
programs you can compile it by

.. code-block:: bash

    % mknemo snapprint

Template
~~~~~~~~

If you need to write a new program in NEMO, you can always clone an existing
program and modify it, if it fits the workflow.   Another approach is to
use the ``template`` script that does all the initial tedious work of
writing your *nemo_main*.  It also can write the initial manual page. Here
is an example:

.. code-block:: bash

    % $NEMO/src/scripts/template foobar a=1 b=2.3 n=10 m=10
    % mknemo foobar
    % $NEMO/src/scripts/mkman foobar > $NEMO/man/man1/foobar.1

After some editing, compiling, testing those two files are ready for inclusion
in the package via a git commit!

Manual pages
~~~~~~~~~~~~

It is very important to keep a manual file (preferably in the UNIX man format)
online for every program. A program that does not have an accompanying manual
page is not complete. Of course there is always the inline help
``help=``) that every NEMO program has. We have a script
``checkpars.py`` that will check if the parameters listed in ``help=``
are the same (and in the same order) as the ones listed in the MAN page.

To a lesser degree this also applies to the public
libraries. A template roff sample can be found in ``example.8``.
We encourage authors to have a MINIMUM set of sections in a man-page as
listed below. The ones with a '*' are considered somewhat less 
important:

- **NAME** the name of the program, or whichever applies
- **SYNOPSIS**  command line format or function prototype, include
  files needed etc.
- **DESCRIPTION**
  maybe a few lines of what it does, or not does.
- **PARAMETERS**
  description of parameters, their meaning and default values.
  This usually applies to programs only.
- **EXAMPLES**
  in case non-trivial, but recommended anyhow
- **DEBUG**
  at what ``debug`` levels what output appears.
- **SEE ALSO**
  references to similar functions, more info
- **BUGS**
  one prefers not to have this of course
- **TIMING**
  performance, dependence on parameters if non-trivial
- **STORAGE**
  storage requirements - mostly of importance when programs
  allocate memory dynamically, or when applicable for the
  programmer.
- **LIMITATIONS**
  does it have any obvious limitations?
- **AUTHOR**
  who wrote it (a little credit is in its place) and/or who
  is responsible.
- **FILES**
  in case non-trivial
- **HISTORY**
  date, version numbers, why updated, by whom (when created)


Makefiles
~~~~~~~~~

Makefiles are scripts in which "the rules are defined to make targets", see
*make(1)* for many more details. In other words, the Makefile tells
how to compile and link libraries and programs. NEMO uses Makefiles
extensively for installation, updates and various other system
utilities. Sometimes scripts are also available
to perform tasks that can be done by a Makefile.

There are several types of Makefiles in NEMO:


- 1. The first (top) level Makefile. It lives in NEMO's root
  directory (normally referred to as ``$NEMO``)
  and can steer installation on any of a number of selected machines,
  it includes some import and export facilities (tar/shar)
  and various other system maintenance utilities.
  At installation it makes sure all directories are present,
  does some other initialization, and then calls Makefile's 
  one level down to do the rest of the dirty work. The top level
  Makefile is not of direct concern to an application programmer,
  nor should it be modified without consent of the NEMO
  system manager.

- 2. Second level Makefiles, currently in ``$NEMO/src``
  and ``NEMO/usr``,
  steer the building of libraries and programs by calling Makefiles in
  subdirectories one more level down.  Both this 2nd level Makefile and
  the one described earlier are solely the responsibility of NEMO system
  manager.  You don't have to be concerned with them, except to know of 
  their existence because your top level Makefile(s) must be callable
  by one of the second level Makefiles. This interface will be described 
  next.

- 3. Third level Makefiles live in source or user directories
  ``$NEMO/src/topic`` and ``$NEMO/usr/name`` (and possibly below). 
  They steer the installation of user specific 
  programs and libraries, they may update NEMO libraries too.
  The user writes his own Makefile, he usually splits up his directory 
  in one or more subdirectories, where the real work is done by what
  we could then call level 4 or even level 5 Makefiles. However, 
  this is completely the freedom of a user.
  The level 3 Makefiles normally have two kinds of entry points
  (or 'targets'): the user 'install' targets are used by the user, 
  and make sure this his sources, binaries, libraries, include files etc. are
  copied to the proper places under ``$NEMO``.
  The second kind of entry point
  are the 'nemo' targets and never called by you, the user;
  they are only called by Makefiles
  one directory level up from within ``$NEMO`` below during the
  rebuilding process of NEMO, i.e. a user never calls 
  a nemo target, NEMO will do this during its installation.
  Currently we have NEMO install itself in two phases, resulting in
  two 'nemo' targets: 'nemo_lib' (phase 1) and 'nemo_bin' (phase 2). 
  A third 'nemo' target must be present to create a lookup
  table of directories and targets for system maintenance. 
  This target must be called 'nemo_src', and must also call lower
  level Makefiles if applicable.

- **Testfile** : this is a Makefile for testing (the ``make test`` command will use it)

- **Benchfile** : this is a Makefile for benchmarks. Under development.
  
This means
that user Makefiles **MUST** have at least these three targets in order
to rebuild itself from scratch. In case a user
decides to split up his directories, the Makefiles must also
visit each of those directories and make calls through the same
entry points 'nemo_lib' and 'nemo_bin', 'nemo_src'; 
a sort of hierarchical install process.
   
For more details see the template Makefiles in NEMO's sec subdirectories
and the example below in section~\ref{ss-example}.

We expect a more general install mechanism with a few more strict rules for
writing Makefiles, in some next release of NEMO. 


An example NEMO program
~~~~~~~~~~~~~~~~~~~~~~~

Under Table~\ref{src:hello}
below you can find a listing of a very minimal NEMO program, 
``hello.c``:

..   $NEMO/src/tutor/hello/hello.c}

.. include:: hello.src
   :literal:

and a corresponding example Makefile to install by *user* and 
*nemo* could look like the one shown under Table~\ref{src:makefile}

Note that for this simple 
example the ``Makefile`` actually larger than
the source code, ``hello.c``, itself. Fortunately 
not every programs needs
their own Makefile,  in fact most programs can be compiled with
a default rule, via the ``bake``
script. This generic makefile is used by the ``bake`` command, 
and is normally installed in ``$NEMOLIB/Makefile``, but check
out your ``bake`` command or alias.

.. Sample makefile - cf. \$NEMOLIB/Makefile}

.. include:: makefile.src
   :literal:



Programming in C++
------------------

Most relevant header files from the NEMO C libraries have been made
entrant for C++. This means that all routines should be available 
through:

.. code-block::

    extern "C"  {    
        .....    
    }


The only requirement is of course that the ``main()`` be in C++.
For this you have to link with the NEMO++ library **before**
the regular NEMO library. So, assuming your header (-I) 
and library (-L) include flags have been setup, you should be able to
compile your C++ programs as follows:

.. code-block::

    % g++ -o test test.cc -L$NEMOLIB -lnemo++ -lnemo -lm


Programming in FORTRAN
----------------------


Programming in FORTRAN can also be done,
but since NEMO is written in C and there is no 'standard' way to link
FORTRAN and C code, such a description is always bound to be system
dependent (large differences exist between UNIX, VMS, MSDOS, and
UNICOS is somewhat of a peculiar case).
Even within a UNIX environment there are a number of ways how the
industry has solved this problem (cf. Alliant). Most comments that
will follow, apply to the BSD convention of binding FORTRAN and C.

In whatever language you program, 
we do suggest that the startup of the program is done in C,
preferably through the ``nemo_main()`` function.
As long as file I/O is
avoided in the FORTRAN routines, character and boolean 
variables are avoided in arguments of C callable FORTRAN functions, 
all is relatively simple. Some care is also needed for multidimensional
arrays which are not fully utilized.
The only thing needed are C names of the FORTRAN routines to be called
from C.  This can be handled automatically by a macro package. 

Current examples can be found in the programs
``nbody0`` and ``nbody2``. In both cases data
file I/O is done in C in NEMO's  *snapshot(5NEMO)* format,
but the CPU is used in the FORTRAN code. 

Examples of proposals for other FORTRAN interfaces can be found
in the directory ``$NEMOINC/fortran``}.

Again this remark: the *potential(5NEMO)* assumes for now a BSD type
f2c interface, because character variables are passed.  This has not
been updated yet.  You would have to provide your own f2c interface to
use FORTRAN potential routines on other systems. 

Simple FORTRAN interface workers within the snapshot interface are
available in a routine ``snapwork(n,m,pos,vel,...)``.

Calling NEMO C routines from FORTRAN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The NEMO user interface, with limited capabilities, is 
also available to FORTRAN
programmers. First of all,
the keywords, their defaults and a help string must be made
available.
This can be done by supplying them as comments
in the FORTRAN source code, as is show in the following 
example listed under Table~\ref{src:testf2c}


..  $NEMO/src/kernel/fortran/test.f}

.. include:: testf2c.src
   :literal:


The documentation section between ``C+`` and ``C-`` can be extracted
with a NEMO utility, ``ftoc``, to the appropriate C module as follows:

.. code-block:: bash

    % ftoc test.f test_main.c

after which the new ``test_main.c`` file merely has to be included on
the commandline during compilation. To avoid having to include
FORTRAN libraries explicitly on the commandline, easiest is
to use the ``gfortran`` command, instead of ``gcc``:

.. code-block::

    % gfortran -o test test.f test_main.c -I$NEMOINC -L$NEMOLIB -lnemo

This only works if your operating supports mixing
C and FORTRAN source code on one commandline. Otherwise try:

.. code-block:: bash

    % gcc -c test_main.c
    % gfortran -o test test.f test_main.o -L$NEMOLIB -lnemo

where the NEMO library is still needed to resolve the
user interface of course.


Calling FORTRAN routines from NEMO C
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


No official support 
is needed, although for portablility it would
be nice to include a header file that maps the symbol names
and such. 


Debugging
---------


Apart from the usual debugging methods (gdb, ddd etc.),  NEMO
programs usually have the following additional properties which can cut
down in debugging time. If not conclusive during runtime, you can either
decide to compile the program with debugging flags turned on, and run
the program through the debugger, or add more ``dprintf``, ``warning``
or ``error`` function calls:

- During runtime you can set the value for the ``debug=``
  (or use the equivalent ``DEBUG`` environment variable)
  system keyword to increase the amount of output. Note that only levels
  0 (the default) through 9 are supported. 9 should produce a lot of output.
  In bash the following two constructs are equivalent:

.. code-block:: bash

    % DEBUG=2 mkplummer . 10
    % mkplummer . 10 debug=10		

  
- During runtime you can set the value for the ``error=``
  (or use the equivalent ``ERROR`` environment variable)
  system keyword to bypass a number of fatal error messages that you
  know are not important. For example, to overwrite an existing file
  you would need to increase ``error`` by 1.


Examples
--------

Two examples, one with existing code (easier), and one about adding new code.

Existing Code
~~~~~~~~~~~~~
To hack an existing code, no new directories, Makefiles or the like need to be added
or modified.

Lets say we want to add a new feature to a program, lets say 
``ccdprint``.   The first step would be to find the code and confirm it still compiles.
For this the ``mknemo`` script is probably the easiest:

.. code-block:: bash

   % mknemo ccdprint
   MKNEMO> Searching ccdprint.c: 
   found one: /home/teuben/NEMO/nemo/src/image/io/ccdprint.c
   gcc -g ...

where it found the code in the ``src/image/io`` directory of NEMO, after which
it got compiled with `gcc`.
In that directory
there will also be Makefile, so an alternative is to do the usual edit-compile-debug
cycle in that directory:

.. code-block:: bash

   %1 cd $NEMO/src/image/io
   %2 edit ccdprint.c
   %3 make ccdprint
   %4 ccdgen junk.ccd
   %5 ccdprint junk.ccd ...

   %6 mknemo ccdprint

Where the final ``mknemo`` was meant to copy the new binary to ``$NEMOBIN`` for
general usage.

Now lets say the this version of ``ccdprint`` need a new function in the image library.
This is located in ``$NEMO/src/image/io/image.c``, which was conveniently in the same
directory as ``ccdprint.c``. Apart from editing the library code, the corresponding
``$NEMOINC/image.h`` probably also needs to know about this new function.  So now the
debug cycle is a two-step proccess:

.. code-block:: bash

   %1 edit $NEMOINC/image.h
   %2 edit image.c
   %3 make install

   %4 edit snapprint.c
   %5 make snapprint
   %6 ccdprint junk.ccd ...

   %7 mknemo ccdprint		

where the new feature is being tested out, and finally
submitted for general usage.  For any important new features, it might be useful
to add this to the self-explanatory ``Testfile`` in this directory.


.. code-block:: bash

   %1 edit Testfile
   %2 make -f Testfile ccdprint
		
Finally, some man pages may need an update

.. code-block:: bash

   %1 edit $NEMO/man/man1/ccdprint.1
   %2 edit $NEMO/man/man3/image.3

New Code
~~~~~~~~

A new program, or a new code for the library is something that can be discovered from looking
at the library. Lets take the example of adding ``xyzio.c`` to the NEMO library. Looking at
the library, this will likely be adding new xyz related references to the Makefile

.. code-block::

   INCFILES = image.h matdef.h xyio.h                            xyzio.h
   SRCFILES = image.c ccddump.c xyio.c                           xyzio.c
   OBJFILES = image.o xyio.o wcsio.o                             xyzio.o
   LOBJFILES= $L(image.o) $L(xyio.o) $L(wcsio.o)                 $L(xyzio.o)
   BINFILES = ccddump ccdprint ccdslice sigccd ccdspec ccdhead   ccdxyz

and testing/compiling

.. code-block:: bash

   %1 edit Makefile
   %2 edit $NEMOINC/xyzio.h
   %3 edit xyzio.c
   %4 make install

   %5 edit ccdxyz.c
   %6 mknemo ccdxyz

and some man page updates:

   %7 $NEMO/src/scripts/mknemo ccdxyz > $NEMO/man/man1/ccdxyz.1
   %7 edit $NEMO/man/man1/ccdxyz.1
   %7 edit $NEMO/man/man3/xyzio.3


Extending NEMO environment
~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now summarize the steps to follow to add and/or create new software to
NEMO.  The examples below are suggested steps taken from adding
Aarseth's ``nbody0`` program to NEMO, and we assume him to have his
original stuff in a directory ``~/nbody0``. 


- 1: Create a new directory, ``"cd $NEMO/usr ; mkdir aarseth"``
        and inform the system manager of NEMO that a new user should be
        added to the user list in ``$NEMO/usr/Makefile``. You can also
        do it yourself if the file is writable by you.

- 2: Create working subdirectories in your new user directory,
        ``"cd aarseth ; mkdir nbody0"``.

- 3: Copy a third level Makefile from someone else, and substitute
        the subdirectory names to be installed for you, i.e. your new
        working subdirectories (``'nbody0'`` in this case):
        ``"cp ../pjt/Makefile . ; edit Makefile"``.

- 4: Go 'home' and install, ``"cd ~/nbody0 ; make install"``, assuming
        the Makefile there has the proper install targets. Check the
        target Makefile in the directory 
        ``$NEMO/usr/aarseth/nbody0`` what this last command
        must have done.

Actually, only step 1 is required. If a user cannot or does not want 
to confirm to the level 3/4 separation, he may do so, as long as the 
Makefile in level 3 (e.g. ``$NEMO/usr/aarseth/Makefile``) contains the 
nemo_lib, nemo_bin and nemo_src install targets. 
which has it's own internal structure. 
