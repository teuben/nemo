      subroutine masset
c written by Jerry Sellwood for the galaxy simulation code
c
c Assigns the masses of all the particles to the grid
c
c Because mass assignment is part of a normal time step, this routine is not
c   called at every step.  It is needed only during the start-up phase.
c
c Note that the coordinates retrieved by GETGRP are put into the oldc array,
c   while MASSGN assumes they are in the newc array.  The coordinates have
c   to be copied from oldc to newc between calls to GETGRP and MASSGN
c
c The final mass array is rescaled by the mass of each particle (in internal
c   program units) rather than scaling each particle as it is added.
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer i, is, jst
c
c set mass array to zero
      call masscl
c assign mass of all lists except that of particles off the grid
      do ilist = 1, nlists - 1
        izone = ilist - 1
        inext = islist( 1, ilist )
c work through groups of particles
        do while ( inext .ge. 0 )
          call getgrp( jst )
c copy data to new coordinates array
          do is = 1, jst
            do i = 1, ncoor
              newc( i, is ) = oldc( i, is )
            end do
          end do
c assign mass of particles
          call massgn( jst )
        end do
      end do
      return
      end
