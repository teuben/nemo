      subroutine tmcgrp( jst, start )
c written by Jerry Sellwood for the galaxy simulation code
c
c Time centres velocities (or undoes this process)
c
c If the logical variable start is .TRUE. then the velocities are backed up
c   half a step.  The acceleration components for all particles in this
c   group should have been determined by a call to GETACC before entering
c   this routine.
c
c calling argument
      integer jst
      logical start
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer i, is
c
      do is = 1, jst
c copy positions and adjust velocities
        do i = 1, 3
          newc( i, is ) = oldc( i, is )
          if( start )then
            newc( i + 3, is ) = oldc( i + 3, is ) - .5 * acc( i, is )
          else
            newc( i + 3, is ) = oldc( i + 3, is ) + .5 * acc( i, is )
          end if
        end do
       end do
      return
      end
