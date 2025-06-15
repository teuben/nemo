# Plotting in NEMO

Plotting in NEMO is generally done by using the *yapp* API. At installation time
NEMO will be compiled and linked with a specific choice, e.g.

      ./configure --with-yapp=pgplot

but most current yapp implementations (in $NEMO/src/kernel/yapp) are either old
or not very flexible. Currently pgplot is the most used, since it can handle
a variety of output formats, including plotting in different window panels. However,
many people now use matplotlib, which has some big advantages:

* can plot to different screens
* plot can be zoomed/panned (but still needs a connection to the program)
* saving in many formats
* support pipeline batch mode


What options do we have?

* abandon all plotting and depend on 3rd party plotting (matplotlib comes to mind).
  This seems to be what ZENO and AMUSE have done.
* write a new yapp.   SVG ?  https://github.com/teuben/nemo/issues/166
* write a matplotlib based yapp, so we can use dynamic pan/zoom/save
* fix plplot?   https://github.com/teuben/nemo/issues/21

Guiding principles:

* do one thing, and do it well
* low entry barrier, one or two command line arguments get you going

## tabplot

The most common plotting program in NEMO is **tabplot**, but cannot easily combine
data from different tables. However, it can plot multiple columns, paired up as well.
It also has a large number (>30) of keywords, which
makes it hard to use, and prone to bugs. 

To overcome the single table approach, the NEMO *getparam* command line user interface
would need to be abandoned (short of rewriting it). The standard *parseargs* module
in python - with a minor hack - will allow parsing the commandline in sections
identified with the input table

Below is that new approach to tabplot, currently being implemented in **tabplot3.py**
and (soon) available for testing.

The originating issue is written up here:  https://github.com/teuben/nemo/issues/87


```
tabplot3.py -i table1 [options] [-i table2 ...] ...

One (or more) tables can be plotted in a scatter diagram using matplotlib

* options per input table

-x --xcol     column number(s) for X axis. Some tables might allow colnames?
-y --ycol     column number(s) for Y axis
-X --dxcol    column(s) with errors in X. Must match the X set. 0 to skip
-Y --dycol    column(s) with errors in Y. Must match the Y set. 0 to skip
              -y 2,3 -Y 4,0   but same as   -y 2,3 -Y 4
-c --color    color (default should be black? or should we follow MPLs default order)
-l --line     line type (default should be no connecting lines)
-p --point    point type (default should be a dot)
-L --legend   add legend

* global options

--xlab        Label along X axis
--ylab        Label along Y axis
--xrange      Min/Max along X
--yrange      Min/Max along Y
--title       label along Y axis
--bigtitle    label along Y axis
--out         plotting file, default is on screen

Other thoughts:
? logarithmic axes
? object oriented plotting vs. pyplot
? allow parameter file for other global things (fonts, ...)
  or rely on rcParams
? need a generic class to read table, like panda's ?
? need a generic class to contain the scatterplot parameters (-x...-L)
? template generator (current tabplot has the pyplot= keyword for this)
? make classes and functions available for derivate products?
  (e.g. a spectrum plotter where one can switch between freq/wave/velo)

Examples:

Generate some tables with 10 rows and 4 columns, sorted in the first columns:

   tabgen - 10 4 seed=1 | sort -n > tab1
   tabgen - 20 6 seed=2 | sort -n > tab2
   tabgen - 30 8 seed=3 | sort -n > tab3


./tabplot3.py -i tab1 -c red \
              -i tab2 -y 2,3 -l 1 -c green \
              -i tab3 -c blue -p 10 \
	      --out fig1.png 

```

###   Table and Scatter classes?

Python classes, one for reading tables (or is pandas good enough?) and a scatter plot
(since each table can do one) seems a good way to modularize the code.

## Other codes and papers 

* stilts has some interesting command line approaches:
  [plot2d](https://www.star.bris.ac.uk/mbt/stilts/sun256/plot2d.html)
  and the more generic
  [plot2plane](https://www.star.bris.ac.uk/mbt/stilts/sun256/plot2plane.html) -
  they refer to it as "old style" plotting.


* https://www.nature.com/articles/s41556-025-01684-z

* https://hyperfit.readthedocs.io/en/latest/
