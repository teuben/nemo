      subroutine stogrp( jst )
c written by Jerry Sellwood for the galaxy simulation code
c
c Scatters the current group of particles back to / ptcls /
c
c Scatters the current group of particles in / buffer / back to the appropriate
c   locations in / ptcls / and makes links for the new list to which each
c   belongs.
c
c See the comments in getstr for the meaning of the / ptcls / variables
c
c calling argument - input the number of particles to be processed
      integer jst
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer i, is, j
c
c scatter particles back to storage area
      do is = 1, jst
        j = loc( is )
        do i = 1, ncoor
          ptcls( j + i ) = newc( i, is )
        end do
c save flag
        iptcls( j + ncoor + 1 ) = iflag( is )
c make links for new list - particles in intensive care
        if( iz( is ) .eq. 0 )then
          i = 1
c particles off the grid
        else if( iz( is ) .eq. nlists )then
          i = nlists
        else
c list number determined by current plane
          i = iz( is ) + 1
        end if
        iptcls( j + nwpp ) = islist( 2, i )
        islist( 2, i ) = loc( is )
      end do
      return
      end
