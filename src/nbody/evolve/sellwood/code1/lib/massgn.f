      subroutine massgn( jst )
c written by Jerry Sellwood for the galaxy simulation code
c
c Assigns the masses of the current group of particles to the grid
c
c It is necessary to test for particles off the grid, since this routine is
c   normally called after the particles have been stepped forward, and they
c   could have left the grid
c
c calling argument
      integer jst
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer is, j, k, l
c
c determine cell nos and weights
      call weight( jst, .false. )
c flag massless particles as off the grid
      do is = 1, jst
        if( iflag( is ) .eq. 3 )nskip( is ) = .false.
      end do
c grid methods
      do is = 1, jst
c skip particles off grid
        if( nskip( is ) )then
c set pointers
          j = ipt( 1 ) + ngxy * iz( is ) + ncl( is )
          k = j + ngx
c distribute mass
          grids( j ) = grids( j ) + wt( 1, is )
          grids( k ) = grids( k ) + wt( 3, is )
          grids( j + 1 ) = grids( j + 1 ) + wt( 2, is )
          grids( k + 1 ) = grids( k + 1 ) + wt( 4, is )
          j = j + ngxy
          k = k + ngxy
          grids( j ) = grids( j ) + wt( 5, is )
          grids( k ) = grids( k ) + wt( 7, is )
          grids( j + 1 ) = grids( j + 1 ) + wt( 6, is )
          grids( k + 1 ) = grids( k + 1 ) + wt( 8, is )
        end if
      end do
      return
      end
