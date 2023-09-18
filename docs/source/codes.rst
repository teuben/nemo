.. _codes:

Related Codes (*)
=================

Here we summarize some codes used in stellar dynamics that are similar to NEMO.
We only list codes that are (publically) available. SPH/hydro codes are currently
not included.   See also `ASCL <https://ascl.net>`_ to find more codes.

First we start off with some expanded examples on a few specific codes
that have a tighter connection to NEMO:

.. include:: amuse.rst

.. include:: martini.rst

.. include:: clustertools.rst
   
.. include:: zeno.rst

.. include:: starlab.rst   

.. include:: yt.rst



List of Related Codes
---------------------

.. todo::   this list needs to be annotated and spiced up with links. 

   
.. glossary::


  ACS
      *The Art of Computational Science - How to build a computational lab.* In C++ and ruby.
      With ideas from NEMO and StarLab.
      | http://www.artcompsci.org/

  agama
      *Action-based galaxy modeling framework*.
      Also usable via ``$NEMO/usr/agama``.
      | http://ascl.net/1805.008
  
  AMIGA
      *Adaptive Mesh Investigations of Galaxy Assembly*.
      | http://ascl.net/1007.006
  
  AMUSE
      *Astrophysical Multipurpose Software Environment*. Also usable via ``$NEMO/usr/amuse``.
      | http://ascl.net/1107.007

  arepo
      *Cosmological magnetohydrodynamical moving-mesh simulation code*.
      | http://ascl.net/1909.010

  bhint
      High-precision integrator for stellar systems
      | https://ascl.net/1206.005
  
  bonsai
      N-body GPU tree-code, in :term:`AMUSE`.
      | http://ascl.net/1212.001
  
  brutus
      See also :term:`AMUSE`

  Catena
      Ensemble of stars orbit integration
      | https://ascl.net/1206.008
      

  CGS
      Collisionless Galactic Simulator. 
      Also usable via :term:`NEMO` with the
      `runCGS <https://teuben.github.io/nemo/man_html/CGS.1.html>`_
      interface.
      | http://ascl.net/1904.003
  
  ChaNGa
      Charm N-body GrAvity solver
      | http://ascl.net/1105.005 	

  clustertools
      A Python package with tools for analysing star clusters. 
      | https://github.com/webbjj/clustertools
  
  DICE
      Disk Initial Conditions Environment
      | https://ascl.net/1607.002
   
  Fewbody
      Numerical toolkit for simulating small-N gravitational dynamics
      | http://ascl.net/1208.011 	
   
  fractal
      *A parallel high resolution Poisson solver for an arbitrary distribution of particles*.
      | https://github.com/jensvvillumsen/Fractal 

  gadgetX
      *A Code for Cosmological Simulations of Structure Formation*.
      Several versions available, X=1,2,3,4.   gadget2 also available via ``$NEMO/usr/gadget``.
      | http://ascl.net/0003.001 	

  Gala
      *Galactic astronomy and gravitational dynamics*.
      | http://ascl.net/1302.011
  
  galaxy
      *N-body simulation software for isolated, collisionless stellar systems*.
      The older version still usable via :term:`NEMO`  with the **rungalaxy** interface.
      | http://ascl.net/1904.002
  
  galpy
      Galactic dynamics package (python)
      | http://ascl.net/1411.008
  
  GalPot
      Galaxy potential code
      | http://ascl.net/1611.006
      
  GANDALF
      Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
      | http://ascl.net/1602.015
  
  GENGA
      Gravitational ENcounters with Gpu Acceleration
      | http://ascl.net/1812.014

  Glnemo2
      Interactive Visualization 3D Program
      | http://ascl.net/1110.008

  gnuastro
      GNU Astronomy Utilities
      | https://www.gnu.org/software/gnuastro/

  GraviDy
      Gravitational Dynamics. Also usable via :term:`NEMO` with the **rungravidy** interface.
      | http://ascl.net/1902.004 	
  
  gsf
      galactic structure finder
      | http://ascl.net/1806.008 	
      
  gyrfalcON
      Dehnen's code. Included in :term:`NEMO`
      | http://ascl.net/1402.031 	
  
  hermite
      In :term:`AMUSE`.

  HiGPUs
      *Hermite's N-body integrator running on Graphic Processing Units*.  Part of :term:`AMUSE`.
      | https://ascl.net/1207.002
   
  
  HUAYNO
      Hierarchically split-Up AstrophYsical N-body sOlver N-body code. Part of :term:`AMUSE`.
      | http://ascl.net/2102.019
  
  Hydra
      A Parallel Adaptive Grid Code
      | http://ascl.net/1103.010
  
  hnbody
      also
      | http://ascl.net/1201.010

  ICICLE
      Initial Conditions for Isolated CoLlisionless systems      https://ascl.net/1703.012
      | http://ascl.net/1703.012

  identikit
      1: A Modeling Tool for Interacting Disk Galaxies.
      | https://ascl.net/1011.001
      2: An Algorithm for Reconstructing Galactic Collisions.
      | https://ascl.net/1102.011

  InitialConditions
      Website with a collection of programs for integrating the equations of motion for N objects,
      implemented in many languages, from Ada to Swift.
      | http://www.initialconditions.org/codes
   
  JSPAM
      Interacting galaxies modeller
      | http://ascl.net/1511.002
  
  limepy
      Lowered Isothermal Model Explorer in PYthon.
      | https://ascl.net/1710.023

  
  MARTINI
      Mock spatially resolved spectral line observations of simulated galaxies
      Also usable via ``$NEMO/usr/martini``, see `example <https://teuben.github.io/nemo/examples/eagle.html>`_.
      | http://ascl.net/1911.005
      
  mcluster
      Make a plummer.
      Also usable via :term:`NEMO`
      | http://ascl.net/1107.015
  
  McScatter
      Three-Body Scattering with Stellar Evolution
      | http://ascl.net/1201.001
      
  mercury
      *A software package for orbital dynamics*. In :term:`AMUSE`.
      | http://ascl.net/1209.010

  montage
      An Astronomical Image Mosaicking Toolkit.
      | https://github.com/Caltech-IPAC/Montage

  MYRIAD
      N-body code for simulations of star clusters
      | https://ascl.net/1203.009

  nbodyX
      Where X=0,1,2,3,4,5,6,6++,7
      Also usable via :term:`NEMO`  with the
      `runbodyX <https://teuben.github.io/nemo/man_html/runbody1.1.html>`_
      interface.
      | http://ascl.net/1904.027
  
  nbody6tt
      Tidal tensors in N-body simulations
      | http://ascl.net/1502.010 	
  
  nbodykit
      Massively parallel, large-scale structure toolkit
      | http://ascl.net/1904.027
      
  nbody6xx
      Alias for nbody6++
      Also usable via :term:`NEMO`  
      | http://ascl.net/1502.010

  N-BodyShop
      Simulation codes of astrophysical phenomenon with particle methods.
      (tipsy, changa, gasoline)
      | https://github.com/N-BodyShop
      
  nemesis
      Another code to document.
      | http://ascl.net/1010.004 	
  
  NEMO
      Stellar Dynamics Toolbox. Our current version is 4.
      | http://ascl.net/1010.051

  NIGO
      Numerical Integrator of Galactic Orbits
      | https://ascl.net/1501.002
  
  N-MODY
      A  code for Collisionless N-body Simulations in Modified Newtonian Dynamics
      | http://ascl.net/1102.001
      
  octgrav
      in :term:`AMUSE`.
      | http://ascl.net/1010.048

  partiview
      Immersive 4D Interactive Visualization of Large-Scale Simulations 
      | https://ascl.net/1010.073
   
  PENTACLE
      Large-scale particle simulations code for planet formation
      | http://ascl.net/1811.019
      
  petar
      Another code to document.
      | http://ascl.net/2007.005 	
  
  plumix
      another
      | http://ascl.net/1206.007
  
  pNbody
      A python parallelized N-body reduction toolbox    https://ascl.net/1302.004

  
  pycola
      N-body COLA method code  
      | http://ascl.net/1509.007
  
  pyfalcon
      Python interface for :term:`gyrfalcON`  
      | https://github.com/GalacticDynamics-Oxford/pyfalcon
  
  pynbody
      N-Body/SPH analysis for python. 
      | http://ascl.net/1305.002

  QYMSYM
      A GPU-accelerated hybrid symplectic integrator
      | https://ascl.net/1210.028

  RAMSES
      A new N-body and hydrodynamical code 
      | https://ascl.net/1011.007
  
  Raga
      Monte Carlo simulations of gravitational dynamics of non-spherical stellar systems
      | http://ascl.net/1411.010

  rebound
      Multi-purpose N-body code for collisional dynamics
      | https://ascl.net/1110.016
      Also usable via :term:`NEMO`

  
  SecularMultiple
      Hierarchical multiple system secular evolution model
      | http://ascl.net/1909.003
  
  sidm-nbody
      Monte Carlo N-body Simulation for Self-Interacting Dark Matter
      | http://ascl.net/1703.007 	
  
  slimplectic
      Discrete non-conservative numerical integrator
      | http://ascl.net/1507.005
  
  smalln
      in :term:`AMUSE`.
      | http://ascl.net/1106.012
  
  smile
      orbits?
      | http://ascl.net/1308.001
  
  SpaceHub
      High precision few-body and large scale N-body simulations  
      | http://ascl.net/2104.025
  
  SpheCow
      Galaxy and dark matter halo dynamical properties
      | http://ascl.net/2105.007
  
  Starlab
      Also usable via :term:`NEMO`  
      | https://ascl.net/1010.076

  Swarm-NG
      Parallel n-body Integrations
      | https://ascl.net/1208.012
  
  Torch
      Coupled gas and N-body dynamics simulator
      | http://ascl.net/2003.014

  TPI
      Test Particle Integrator
      | http://ascl.net/1909.004

  UNSIO
      Universal Nbody Snapshot I/O - 
      See `examples <https://teuben.github.io/nemo/examples/uns.html>`_.
      
  VINE
      A numerical code for simulating astrophysical systems using particles
      | http://ascl.net/1010.058
  
  yt
      A Multi-Code Analysis Toolkit for Astrophysical Simulation Data
      | https://ascl.net/1011.022
  
  ZENO
      Barnes version that was derived from NEMO V1.
      | https://ascl.net/1102.027
      Also usable via :term:`NEMO` , but watch out for duplicate names of programs


A large number of these codes can also be found by searching on ASCL, for example:
https://ascl.net/code/search/dynamics
and
https://ascl.net/code/search/hermite
and
https://ascl.net/code/search/orbit
and
https://ascl.net/code/search/nbody.
The last time this list was cross-checked was ... 16-jul-2021.


.. include:: python.rst

Categories
----------

Such a niche list of codes made me wonder what kind of meta-data we could use to
categorize such dynamics codes, but then perhaps along the lines of the
`Unified Astronomy Thesaurus <https://github.com/astrothesaurus>`_ project.


dynamics - nbody, orbit, integrator, sph, hydro, analysis, integrator


