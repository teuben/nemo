      program galaxy
c program for simulation of an isolated galaxy written by:
c
c                Jerry Sellwood and Richard James
c
c It is a direct descendant of that first used by James & Sellwood (1978, MNRAS
c   v182, p331) with only minor refinements since.  The Poisson solver was
c   written exclusively by Richard James and uses the algorithm descibed in
c   James (1977, J Comp Phys v25, p71).  This algorithm determines the
c   gravitational potential of an isolated mass distribution using FFTs on a
c   3-D Cartesian grid without the need to remove images by doubling the grid
c   size in each dimension.  The motion of particles within the grid volume is
c   integrated forwards in time according to accelerations determined by
c   differencing the grid potential values.  Particles that leave the grid
c   volume are discarded.
c
c Very large numbers of particles can be employed at little cost since most cpu
c   time is taken up by the determination of the potential on the grid.  [See
c   Sellwood (1997, in "Computational Astrophysics" eds Clarke & West, ASP Conf
c   series v123, p215) for a performance comparison with other codes.]
c
c The version of the algorithm used here requires the individual grid cells to
c   be cubic, but the overall grid need not be cubic.  The FFTs supplied
c   require the number of mesh spaces in each direction to be (2**n + 1), where
c   the exponent n may be chosen independently for each coordinate direction.
c   As large grids require a great deal of memory (c400 MB for a 257**3 grid),
c   it is recommended that the parameters set in the include file 'rjinclude.h'
c   be no larger than necessary, although the code will function correctly as
c   long as the actual dimensions used do not exceed those set by the parameter
c   statement at compile time.
c
c Time integration follows the standard 2nd order time-centered leap-frog, with
c   the velocities one half a time step out of synchrony with positions.  This
c   difference is maintained in the internally stored coordinates and is
c   created, and can be undone, by a call to subroutine TMCENT.  For output of
c   the particle coordinates at a particular instant, the velocities need to be
c   the average of those before and after the time for the positions.
c
c Results, in this public version, are simply the phase space coordinates of
c   all the particles as often as requested, which can create a very large
c   file.  The authors therefore do not employ this scheme themselves,
c   preferring instead to measure and save properties of the model as the
c   simulation evolves.  An example of this "on the fly" analysis is provided
c   in the procedure to determine the global integrals (energy, momentum, etc).
c
c The grid is set up in subroutine GRDSET using data read in from a short ASCII
c   input file (galaxy.dat)
c The positions and velocities of the particles are read in subroutine LOADUP
c The gravitational field is determined by a call to FINDF
c The model is integrated forward by a call to STEP
c After the desired evolution is completed, the positions and velocities of
c   the particles are saved by a call to UNLOAD
c
c The main files associated with the run are:
c   galaxy.dat - ASCII input: grid parameters, length and time scales
c   galaxy.lis - ASCII output: a brief summary of progress (appended)
c   galaxy.ini - ASCII input: initial coordinates of all the particles
c   galaxy.fin - ASCII output: final coordinates of all the particles
c   galaxy.res - binary output: coords and potentials at intervals (appended)
c
c There are also two scratch files:
c   galaxy.aux - short ASCII file (deleted when closed)
c   galaxy.tmp - large binary file (deleted when closed)
c
c common block
c
      include 'admin.h'
c
c make sure block data segment is linked
      external blkdat
c
c local variables
      integer i, ilast, j, l
      real tend

c say hello and remind the user of the max- parameters
      print *, 'galaxy V1.3 '
      print *, 'Maximum number of particles (mbuff): ',mbuff
c
c store the sizes of the main common arrays
      call dimens
c read galaxy.dat file (ASCII input) and se up grid
      call grdset
      read( ni, * )tend
      close( ni )
c create Green's function for Poisson solver
      call greenm( .false. )
c read the galaxy.ini file
      call loadup
      print *, 'Run (re)started at time', ts * real( istep )
      ilast = nint( tend / ts ) + 1
      print *, '  and will stop after the step at time ',
     +         ts * real( ilast )
      if( istep .ge. ilast )call crash( 'GALAXY',
     +                                  'End time already passed' )
c assign masses for the current particle positions
      call masset
c find gravitational field
      call findf
c back up velocities one half step
      call tmcent( .true. )
c main time-step sequence
      do while ( istep .lt. ilast )
        print *, 'Starting step', istep
        write( no, * )'Starting step', istep
c move particles
        call step
c increment step number
        istep = istep + 1
c calculate new gravitational field
        call findf
      end do
      write( no, * )'Run stopped at time', ts * real( istep )
c step forward velocities one half step
      call tmcent( .false. )
c update the .dmp file
      call unload
      stop
      end

      subroutine dimens
c written by Jerry Sellwood for the galaxy simulation code
c
c Sets the dimensions of the two adjustable arrays in the common blocks
c   / grids / and / ptcls /  They are dimensioned to have size 1 in every
c   other routine except for this in order to avoid having to recompile
c   the entire library whenever the size of the problem is changed
c
c The actual sizes of these arrays are stored in the first variable in
c   each of these blocks
c
      integer mbod, par1, par2
      parameter ( mbod = 10000 )
c
      include 'rjinclude.h'
c
c common blocks
c
c size needs to be large enough to store all the particles and pointers
      parameter ( par1 = 8 * mbod )
      integer lptcls
      real ptcls( par1 )
      common / ptcls / lptcls, ptcls
c
c size needs to be large enough to store coordinates plus potential
      real coor( 7 * mbod )
      common / results / coor
c
c size needs to be two whole meshes (for masses and potentials) plus 6 planes
c (one pair for each acceleration component)
      parameter ( par2 = 2 * rjn1 * rjn2 * rjn3 + 6 * rjn2 * rjn3 )
      integer lgrids
      real grids( par2 )
      common / grids / lgrids, grids

c
c say hello and remind user to the size
      print *, 'Maximum gridsize: ',rjn1,rjn2,rjn3
c
c store sizes in common blocks
      lptcls = par1
      lgrids = par2
      return
      end
