.. _table:

Tables are obiquitous in astronomy in general and in NEMO specifically. Often
they are simple ASCII tables, with columns separated by a space or a comma (CSV).
NEMO handles both, though generally output tables are space separated.

Intelligent tables (e.g. IPAC tables, the astropy ECSV format) add
some header type information to annotate the table with
column information (name, type, units) and sometimes
even FITS-like provenance (e.g. IPAC tables).   NEMO has limited
ways to process these, currently loosing much of this header and provenance.

There is no good consistency if an **out=** keyword is needed,
or if the table is written to *stdout* by default,
or if **out=-** is the last parameter.
Otherwise table programs are well suited for working in pipes,
as most NEMO programs are. Some examples are given below, with
the caveat that the output pipe can be inconsistent.

Columns are typically given through the **xcol=** and **ycol=** keywords,
and are 1-based. Column 0 has a special meaning, it is the row number
(again 1 being the first row).

Manipulation
~~~~~~~~~~~~

Here are some programs that manipulate tables:

- **tabgen** was designed to produce tables from scratch, with some type of *random*
  values. Mostly to test performance and scalability of tables. Here a small table
  of 4 (rows) by 3 (columns) is written to  *stdout:

.. code-block::

  % tabgen - nr=4 nc=3
  0.552584 0.126911 0.520753
  0.0930232 0.563683 0.258931
  0.0369577 0.965763 0.634585
  0.981766 0.334841 0.988963

- **tabcols** selects columns from a table

- **tabrows** selects rows from a table

- **tabtab** takes two or more tables, and will paste (increasing the columns)
  or catenate (increasing the rows)

- **tabtranspose** transposes columns and rows.

.. code-block::

  % tabgen - nr=4 nc=3 | tabtranspose  - -
  0.552584 0.0930232 0.0369577 0.981766 
  0.126911 0.563683 0.965763 0.334841 
  0.520753 0.258931 0.634585 0.988963 
   
- **tabcomment**  comments out things that look like comments. Sometimes this filter
  is needed before another table program can be used, as they can get confused with
  non-numeric information. For example picking out two
  columns from a dirty table:

.. code-block::

  % tabcomment mytable.tab | tabcols - 2,3
  ...

- **tabcsv** converts a table  to some XSV, where X can be choosen from a set of
  characters

- **tabmath** is a column calculator, but columns are referenced with a % sign, so
  %2 refers to column 2.. By default it will add the new columns
  to the old columns, unless some selection, or **all**, is/are removed.
  It can also select rows based on the values in a row, for example, only
  select rows where column 2 is positive.

  Here is an example of adding two columns



.. code-block::

  % tabgen -  nr=4 nc=3 | tabmath - - %1+%2 all
  0.679495 
  0.656706 
  1.00272 
  1.31661 

- **txtpar** extracts parameters from a text file. The often complex operations involving
  Unix tools such as grep/awk/sed/cut/head/tail can often be simplified with **txtpar**.
  Here is an example :

.. code-block::

  % cat example.txt
  Worst fractional energy loss dE/E = (E_t-E_0)/E_0 = 0.00146761 at T = 1.28125
  QAC_STATS: - 0.039966 0.0274195 0.00185505 0.135854  0.799319 1  20

  % txtpar in=example.txt expr="log(abs(%1)),log(abs(%3/%2))" format=%.3f \
           p0=Worst,1,9  p1=QAC,1,3 p2=QAC,1,4
  -2.833 -0.164
   

Generating
~~~~~~~~~~

Some other NEMO programs that produce tables:

- **snapprint** tabulates choosen body variables in ascii. The output
  format can also be choosen.

.. code-block::

  % mkplummer - 3 seed=123 | snapprint - comment=t format=%8.5f
  #       x       y       z       vx      vy      vz
  -2.19118 -0.01225  0.18687 -0.21951 -0.14248  0.28165 
   3.22756 -0.27674 -0.88792  0.27564  0.18469  0.11529 
  -1.03638  0.28899  0.70105 -0.05613 -0.04221 -0.39694


- **orblist** tabulates an orbit (*should be renamed to orbprint*)

- **ccdprint** tabulates values of an image


Plotting
~~~~~~~~

- **tabplot** creates a scatterplot of one or more columns

- **tabhist** creates a histogram

- **tabview, tabzoom** dynamic queries (cf. glueviz, topcat)


Analysis
~~~~~~~~

- **tabtrend** computes differences between successive rows

- **tabint** integrate a sorted table

- **tabsmooth**

- **tablsqfit** and **tabnllslqfit** do linear and non-linear fitting of functions

- **tabstat**
