Research Software Engineering
=============================

NEMO was rooted with a Unix philosopy, and by now can be called
legacy. There has been a constant change in the software engineering
(SE) aspects it was using. This page hopes to summarize those,
especially useful for those with no domain knowledge. We also highlight
a few peculiar aspects.

* Like many research programming projects, NEMO grew organically but over
  the years has tried to absorb good modern SE habbits to make the project
  sustainable and organically grow.

* NEMO uses a Unix like directory tree, with a few familiar 3 and 4 letter
  names (src, usr, man, etc, bin, lib, opt, obj, tmp, docs, data, local)

* We use git (on github):
   * since its inception in 1986 we have been sharing the code and absorbed new codes
     from collaborators
   * We started out with email (1987) and *shar* files, then moved to *CVS* in 2000, and finally
     to *git* in 2017, hosted on github for the moment
   * issues (for bug reports) and pull requests are the preferred way to collaborate
   * projects (in github) discuss a few ideas we have where changes could be expected

* installing NEMO occurs in place, we use **configure** to track down the system
  dependent portions (or force them), and a set of **make build** steps will compile
  and build NEMO. There is no **install** step, so users will need to modify their
  *shell environment* to take advantage of the new tools.  There is no support for
  a build tree that's different from the source tree.

* In terms of documentation there are the unix manual (``man``) pages, very hyperlinked,
  and the old *Users and Programmers Guide* (in latex).  The man pages are still used,
  but also made available via html (see https://teuben.github.io/nemo/man_html/index1.html),
  and the newer sphinx based documentation shared on https://astronemo.readthedocs.io/en/latest.
  There is also the experimental TLDR suport and the classic "motd" via the **fortune**
  program.

* Q: Could you explain the software development lifecycle that you use in this project?


* Q: Is the package modularized enough to be improved by other software engineers without knowing the domain knowledge?
  maybe

* Q: Are there sufficient unit tests for the developers?
  some

* Q: Are there benchmark data to test the new implementations?
  yes

* Q: Could you highlight the software engineering aspects of the project? 

* Q: What is your advice for early-career software engineers who want to work on this project?


* Q: If you want to redesign the package, what would you do?
   * writing the core in C was probably good. Portability of C++ is now better. Even our conversion to C99 is not complete.
   * add unit tests, now they are distributed in local Testfile and "make testXXX"
   * better autoconf and/or cmake based install built from the ground up
   * integration in other than bash, eg. python or other high level scripting (julia?)

* a number of peculiar data (e.g. rotcur.5), and benchmark (CPU wise, but also "24+1" body)

* Q: What are lessons learned to keep this code going?

* Q: How many users are using it? There is probably a small core group, and there is good competition from python based
  environments (yt, amuse, galpy, gala, pynbody). In philosphy. `AMUSE <https://www.amusecode.org/>`_ comes closest
  to NEMO.

* Q: How many citations are there?
  This is the sad story about legacy software. Only recently ADS and ASCL make it possible to track software.

* Some innovative aspects of NEMO perhaps not seen widely used:
   * connect code with papers via an `ADS bibcode <https://ui.adsabs.harvard.edu/help/actions/bibcode>`_.
     See our `bibcode <https://teuben.github.io/nemo/man_html/bibcode.html>`_ table
   * help= vs. --man vs. man and **checkpars** 
   * loadobj:   body, potential, ...
   * data are binary, but name and type tagged and hierarchical, a bit like binary XML if you like. The **tsf** program
     will show the contents of such files in a more human readable format
   * most NEMO programs work fully in memory (not so innovative, but it's fast)
   * scattered Testfile to pick up tests for "make check" (like a Makefile)
   * scattered Benchfile to pick up benchmarks for "make bench" (like a Makefile)
   * the use of the **bsf** program to check reproducability to N digits.
   * the user of - (dash) for files in a pipe, and . (dot) for the last pipe (like /dev/null)
   * *surely there must be more*  
