      subroutine stpgrp( jst )
c written by Jerry Sellwood for the galaxy simulation code
c
c advances the motion of the current group of particles
c
c Despite its appearance, this algorithm is properly time symmetric and is
c   exactly equivalent to the usual symplectic leap-frog integrator.  This
c   is because the velocities are initially backed up by half a step and
c   the stored values are always those appropriate for half a step before
c   the values stored for positions.
c
c calling argument
      integer jst
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer i, is
      logical left
c
c advance first velocity and then position for all particles in this group
      do is = 1, jst
        do i = 1, 3
          newc( i + 3 , is ) = oldc( i + 3 , is ) + acc( i, is )
          newc( i, is ) = oldc( i, is ) + newc( i + 3, is )
        end do
      end do
      do is = 1, jst
c set new plane number
        iz( is ) = newc( 3, is ) + zm + 1.
c check for particles leaving the grid
        left = ( abs( newc( 1, is ) ) .gt. xm ) .or.
     +         ( abs( newc( 2, is ) ) .gt. ym ) .or.
     +         ( abs( newc( 3, is ) ) .gt. zm )
        if( left )then
          iz( is ) = nlists
          iflag( is ) = -istep
          noff = noff + 1
        end if
      end do
      return
      end
