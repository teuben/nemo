      subroutine step
c written by Jerry Sellwood for the galaxy simulation code
c
c Advances all the particles forward by one step
c
c The potential over the entire grid should have been determined previously
c   by a call to FINDF
c
c The particles are processed in groups, each group being worked through the
c   following sequence of calls:
c         GETGRP:  collect the next group from the main particle storage area
c         GETACC:  look up acceleration components to be applied
c         STPGRP:  advance both velocities and positions
c         STOGRP:  replace new coordinates in the main particle storage area
c         MASSGN:  assign mass to the grid using new positions
c
c The routine begins with calls to MASSCL to zero out the mass array in
c   readiness for the new values
c After all particles have been processed, a call to SCALED multiplies the
c   mass array by the appropriate normalisation factor
c
c In the 3-D Cartesian code, the particles are sorted by grid planes to enable
c   all those in one plane to be processed before beginning the next.  This
c   is so that acceleration components can be determined from the potential
c   more efficiently by a single call to DIFPOT for each plane.
c The particles in each plane are located through a linked list with the initial
c   location stored in the array ISLIST and the last particle has a negative
c   pointer.  A new linked list is created for the new positions, and new the
c   initial list origins are copied over the original set after all particles
c   have been processed.
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer i, jst
      logical check, out
c
      check = .false.
      if( chkstp .gt. 0 )check = mod( istep, chkstp ) .eq. 0
c initialize global integrals
      if( check )call measure( .false. )
c set mass arrays to zero
      call masscl
c determine whether particles are to be output at this step
      out = .false.
      if( outstp .gt. 0 )out = mod( istep, outstp ) .eq. 0
c initialize new lists and preserve pointer to particles off
      do i = 1, nlists - 1
        islist( 2, i ) = -1
      end do
      islist( 2, nlists ) = islist( 1, nlists )
c step forward all particles except those off the grid
      do ilist = 1, nlists - 1
c set zone for this list
        izone = ilist - 1
        inext = islist( 1, ilist )
        if( inext .ge. 0 )call difpot( ilist - 1 )
c work through groups of particles
        do while ( inext .ge. 0 )
          call getgrp( jst )
c get force components
          call getacc( jst )
c move particles
          call stpgrp( jst )
c determine integrals
          if( check )call icheck( jst )
c prepare particles for output if requested
          if( out )call cenvel( jst )
c store new coordinates
          call stogrp( jst )
c assign mass of particles
          call massgn( jst )
        end do
      end do
c update list origin table
      do i = 1, nlists
        islist( 1, i ) = islist( 2, i )
      end do
c output particles if requested
      if( out )call putout
c evaluate and output global integrals
      if( check )call measure( .true. )
      return
      end
