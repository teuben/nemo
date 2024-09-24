Summary
=======

Here we provide a brief reminder of NEMO's features.

* NEMO uses your login shell in a terminal. To use NEMO you will need
  to source the "nemo_start.sh" (or .csh) file to load NEMO in your
  shell environment. This will add about 300 programs to your shell.
  For some work this could cause side-effects for other
  packages, as they may share the name, or cause conflicting shared
  libraries.

* Programs have a consistent command line interface, featuring --help, --version, --man.
  For a full overview use "help=?".   Each program should have a standard Unix "man" page.
  There are program keywords and system keywords.
  See also the getparam(3NEMO) man page on https://teuben.github.io/nemo/man_html/getparam.3.html

* Special filenames are "-" (stdin/stdout) used in pipes, but a http:// style reference
  can be also used to obtain files directory over a network connection.
  See also stropen(3NEMO) on  https://teuben.github.io/nemo/man_html/stropen.3.html

* Data files in NEMO are generally binary, of a special hierarchically name and type tagged
  format, almost like binary XML, and similar to HDF5. The "tsf" program shows the content
  of such binary files in human readable format.

* ASCII tables are also a popular format in NEMO, and a large set of programs support these
  files, e.g. tabmath, tabplot.

* Programs that plot use the yapp= system keyword to determine the plot output. At installation
  a yapp interface is choosen (e.g. pgplot). A simple postscript driver can be activated
  using "--with-yapp=ps" at configure time.

* (new) programs can be compiled using "mknemo", adding the "-u" flag will also pull the
  most recent version via git. mknemo is also used to install 3rd party packages in $NEMO/opt.

* Most NEMO programs work fully in memory, and are single-core.

* The "Code with Papers" project collects bibcodes for papers that have a corresponding code in NEMO.
  See https://teuben.github.io/nemo/man_html/bibcode.html

* Dynamic Object Loader will provide body variable parsers, potentials for orbit integrations, fitting
  functions etc. See $NEMO/obj

* NEMO follows a fairly traditional Unix source code hierarchy, with typical directories such as
  src, bin, lib, and man, with additional peculiar ones for nemo: obj, opt, local data, text, docs, etc.

* The build system uses Makefile's, and are locally stored per directory and called hierarchically. Common
  parameters are stored in $NEMOLIB/makedefs, which is generated during the configure stage.

* NEMO includes a number of "legacy" codes. There is a "run" interface for legacy codes, to make them
  act like NEMO programs. An example is runbody1, which is the front-end for Aarseth's  NBODY1 code.
   
