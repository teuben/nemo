.. _aiface:

The User Interface (Reference)
==============================

.. todo::  section needs to be merged with the previous one
   
This Appendix overviews the basic command-line interface of NEMO programs. 

Every NEMO program accepts input through a user supplied parameter list
of *keyword=value* arguments. In addition to these program specific
**program keywords**,
there are a number of system wide  defined **system keywords**,
known to every NEMO program.

Program keywords
----------------

Program keywords are unique to a program, and need to be
looked up in the online manual page or by using the 
``help=`` system keyword (dubbed the **inline** help). Parsing of
*values* is usually done, though sometimes primitive. Program
keywords also have the ability to read the value(s) of a keyword from a file
through the ``keyword=@file`` construct. This is called the 
**include keyword file**, and is very handy for long keyword values,
not having to escape shell characters etc.  Newlines are replaced by blanks.

System keywords
---------------

The 'hidden' system keywords, although overridden by 
any program defined counterpart, can also be set by an 
equivalent environment variable (in upper case).



- **help=** Sets the help level to a program. As with all
  system keywords, their value can be fixed for
  a session by setting the appropriate
  environment variable in upper case, *e.g.* ``expor HELP=5``.


  By using the keyword form, the value of the environment variable 
  will be ignored.

  The individual help levels are numeric and add up to combine
  functionality, and are hence powers of 2:

  - ``1`` Remembers previous usage of a program, by
    maintaining a keyword file from program to program.  These files are
    normally stored in the current directory, but can optionally be stored
    in one common directory if the environment variable 
    {\bf NEMODEF}\footnote{mirtool also uses this environment variable} is
    set.  The keyword files have the name {{\it "progname"}{\bf.def}},
    {\it e.g.} {\tt snapshot.def}\footnote{This may result in long
    filenames, Unix SYS5 allows only 14 characters - a different solution is
    needed here}.  When using this lowest help-level it is still possible to
    use UNIX I/O redirection.  This help level
    reads, as well as writes the keyword file during the program execution;
    hence the user needs both 
    read and write permission in the keyword directory.  As can also
    be seen, programs cannot run in parallel while using this help-level: they
    might compete for the same keyword file.
    Within the simple commandline interface it is not possible
    to maintain a global keyword database, as is {\it e.g.}  the case in AIPS;
    you would have to use the {\tt miriad} shell.

  - ``2`` prompts the user for a (new) value for every
    keyword; it shows the default (old) value on the prompt line, which can
    then be edited.  
    It is not possible to combine this level with UNIX I/O redirection. 
    By combining the previous helplevel with this one, previous
    values and modified ones are maintained in a keyword file.

  - ``4`` provides a simple fullscreen menu interface, by having
    the user edit the keyword file. The environment variable
    {\bf EDITOR} can be used to set any other editor than good old 
    {\it vi(1)}.
    It is not possible to combine this level with UNIX I/O redirection. 

  - ``8,16,...`` although not processed, higher powers of 2 are reserved for
    future options


    Example: ``help=3`` will remember old keywords in a local keyword file,
    prompt you with new values, and puts the new values in the keyword file
    for the next time.  The ``help=5`` option happen to be 
    somewhat similar to the way ``AIIPS`` and ``IRAF`` appear to the user. 

    Help levels can also include an alpha-string, which generally display
    the values of the keyword, their default values or their help strings.

  - ``?``
    lists all these options, as a reminder. It also displays the
    version \index{version, user interface} of the 
    {\tt getparam} user interface package.

  - ``h``
    list all the keywords, plus a help string what the keywords does/expects.
    This is really what we call the inline manual or inline 
    help. \index{inline, help} \index{manual, inline} \index{help, inline}

  - ``a``
    list all arguments in the form {\it keyword=value}.

  - ``p,k``
    list parameters (keywords) of all arguments in the form {\it keyword}.

  - ``d,v``
    list defaults (values) of all arguments in the form {\it value}.

  - ``n``
    add a newline to every {\it keyword/value} string on output.
    In this way a keyword file could be build manually by redirecting this
    output.

  - ``t``
    output a documentation file according to the
    \%N,\%A specifications \index{mirtool} of 
    {\tt miriad}\footnote{Both {\tt mirtool} and {\tt miriad} need such a doc-file
    \index{doc file, miriad} to lookup keywords and supply help}.
    Is mainly intended to be used by scripts such as {\tt mktool}. 
    The procedure in NEMO to update a {\tt .doc} file would be:

    .. code-block::

         % program help=t > $NEMODOC/program.doc

  - ``q``
    quit, do not start program. Useful when the helpstring contains
    options to print.

    Example: **key=val help=1q** redefines a keyword in the keywordfile,
    but does not run the program. This is also a way to 'repair' a keyword
    file, when the program has been updated with new keywords.
    **key=val help=1aq** redefines the keyword,
    shows the results but does still not run the program. 
    Finally, **key=val help=1a** redefines a keyword, shows
    the result and then runs the program.


- **debug=**  Changes the debug output level.  
  The higher the debug
  level, the more output can appear on the standard error output device
  ``stderr``.  The default value is either 0 or the value set by the
  **DEBUG** environment variable.  The use of the ``debug=`` keyword
  will override your default setting.  A value of '0' for debug 
  may still show some warning messages.  Setting debug to 
  -1 will prevent even those warning/debug messages.  Legal values are 0
  through 9.  Values of **DEBUG** higher than 9 are not used, or
  you may get some weird screen output. Values larger than
  5 cause an error to coredump, which can then be used with debug utilities
  like *abd(1)* and *gdb(1)*.

- **error=** Specifies how many times the fatal error routine can be
  bypassed. The **ERROR** environment
  variable can also be set for this. The default, if neither of them
  present, is 0.

- **yapp=** Defines the device to which graphics output is send. 
  Currently only interpreted for a limited number of yapp devices.  
  Some yapp
  devices do not even listen to this keyword.  Check *yapp(5NEMO)* or
  your local NEMO guru which one is installed.  The default device is
  either 0 or the value set by the **YAPP** environment variable.

- **np=**  Defines the number of processors (e.g. in an OpenMP setting)
  that can be used. This would override the OMP_NUM_THREADS environment
  variable, if it was present.

- **outkeys=**  TBD

- **argv=**  TBD

yapp_ps
~~~~~~~

By default NEMO is compile with a very simple PostScript device driver, as
specified in yapp_ps. This YAPP interface  produces a simple PS
(supposedly correctly calibrated to be 20 x 20 cm), and 
the yapp= keyword value specifies the PS filename.

yapp_pgplot
~~~~~~~~~~~

The YAPP interface to the common PGPLOT library is the most used
interface, and allow one to select from a variety of graphics output
devices without having to recompile the program.

A graphics device in PGPLOT
is defined by preceding it with a slash
Optional parameters (e.g. filename, X device etc.)
can be supplied before the slash. The following
list gives an overview of some of the available devices
(your list may be a lot shorter (see ``?``) in list below):

.. code-block::

       ?           Get a list of all currently defined graphics devices   
       /XTERM     (XTERM Tek terminal emulator)
       /XWINDOW   (X window window@node:display.screen/xw)
       /XSERVE    (A /XWINDOW window that persists for re-use)
    Non-interactive file formats:
       /NULL      (Null device, no output)
       /PNG       (Portable Network Graphics file)
       /TPNG      (Portable Network Graphics file - transparent background)
       /PS        (PostScript file, landscape orientation)
       /VPS       (PostScript file, portrait orientation)
       /CPS       (Colour PostScript file, landscape orientation)
       /VCPS      (Colour PostScript file, portrait orientation)
       /EPS       (Encapsulated Postscript, colour)


See also manual pages such as *getparam(3NEMO)* and
*yapp(5NEMO)*

The REVIEW section
------------------

.. warning::
   This option will likely be deprecated. By default it is not enabled.

By setting the **REVIEW** environment variable a NEMO program is 
always put into the REVIEW
section just before the start of the actual execution of
the program (the end of the {\it initparam(3NEMO)} routine). 
This functionality is quite similar to using the helplevel
{\tt help=4} (see previous Section).

A NEMO program can also be interrupted, using the quit signal
(see {\it signal(2)}), \index{signal(2)} into the 
REVIEW section, although the program must be adapted to get
keyword information through {\it getparam(3NEMO)} and not through it's
own local database, in order for modified keywords to take effect.
This does not hold for the system keywords, whose new value is always
correctly interpreted by the program.

In the REVIEW section the prompt is {\bf ``REVIEW''} and
the following commands are understood:

- \item{{\bf exit, quit, end}}
  Exit the program (ungracefully).

- \item{{\bf stop}}
  Gracefully end the program, but first goes through {\tt finiparam()} (see
  {\it getparam(3NEMO)}) to update the keyword file if the
  helplevel includes 1.

- \item{{\bf set [key=[value]]}}
  Set a new value for a program keyword ({\tt set key=value}), where
  {\tt value} may also be blank, or display the contents of a 
  program keyword ({\tt set key}).

- \item{{\bf show key}}
  Show the value of a program keyword.

- \item{{\bf keys}}
  Show the values of all program keywords.

- \item{{\bf syskeys}}
  Show the values of all system keywords.

- \item{{\bf set syskey[=value]}}
  Set a new value for a system keyword {\tt set syskey=value}
  or display its current contents {\tt set syskey}. 

- \item{{\bf time}}
  Show the cputime (in minutes) used so far. \index{cputime}

- \item{{\bf !cmd}}
  Pass a command {\tt cmd} to the shell.

- \item{{\bf go,continue}}
  Continue execution of the program.

- \item{{\bf version}}
  Display version of {\tt initparam()} compiled into program.

- \item{{\bf ?, help}}
  Displays all commands and their format.



When the system keyword ``debug`` is non-zero, the ``REVIEW`` prompt also
includes the process identification number of the process.
