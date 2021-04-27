Examples
========

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


N-body experiments
------------------

In this section we will describe how to set up an N-body experiment, run,
display and analyze it.  In the first example, we shall set up a head-on
collision between two spherical "galaxies" and do some simple analysis.

Setting it up
~~~~~~~~~~~~~

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

.. code-block::

    1% mkommod in=$NEMODAT/k7isot.dat out=tmp1 nbody=100 seed=280158
              


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

.. code-block::

    2% snapscale in=tmp1 out=tmp1s rscale=2 "vscale=1/sqrt(2.0)"


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

.. code-block::

    3% mkommod in=$NEMODAT/k7isot.dat out=tmp2 nbody=100 seed=130159


This second galaxy needs to be rescaled too, if you want virial units:


.. code-block::

    4% snapscale in=tmp2 out=tmp2s rscale=2 "vscale=1/sqrt(2.0)"


We then set up the collision by stacking the two snapshots, albeit with
a relative displacement in phase space.  The program {\tt snapstack} was exactly
written for this purpose:\index{snapstack(1)}


.. code-block::

    5% snapstack in1=tmp1s in2=tmp2s out=i001.dat deltar=4,0,0 deltav=-1,0,0


The galaxies are initially separated by 4 unit length and approaching
each other with a velocity consistent with infall from infinity
(parabolic encounter). The particles assembled in the data file
{\tt i001.dat} are now ready to be integrated.

To look at the initials conditions we could use:

.. code-block::

    6% snapplot i001.dat xrange=-5:5 yrange=-5:5

which is displayed in Figure X

Integration using hackcode1
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We then run the collision for 20 time units, with the standard
N-body integrator based on the Barnes 
\index{Barnes J.} \& Hut \index{Hut P.} "hierarchical tree" 
algorithm\footnote{see also their paper in: Nature, Vol. 324, pp 446 (1986).}:
\index{hackcode1(1)}

Images
------

more coming


Tables
------

more coming

Potential
---------

more coming

Orbits
------

more coming

Exchanging data
---------------

more coming
