      subroutine weight( jst, old )
c written by Jerry Sellwood for the galaxy simulation code
c
c Determines the weights for linear interpolation between grid points
c
c Computes cell numbers and weights for the current group of particles from
c   either the new or the old coordinates.  Particles off the grid are flagged.
c
c calling arguments
      integer jst
      logical old
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer is, jx, jy, jz
      real dx, dy, dz, x, y, z
c
c flag particles on the curent grid
      do is = 1, jst
        nskip( is ) = ( iz( is ) .lt. nlists )
      end do
c
      do is = 1, jst
c skip particles off grid
        if( nskip( is ) )then
c compute cell number in complete ( ngx x ngy ) grid plane
          if( old )then
            x = oldc( 1, is ) + xm
            y = oldc( 2, is ) + ym
            z = oldc( 3, is ) + zm
          else
            x = newc( 1, is ) + xm
            y = newc( 2, is ) + ym
            z = newc( 3, is ) + zm
          end if
          jz = z
          jy = y
          jx = x
          ncl( is ) = ( jy + 1 ) * ngx + jx + 2
c compute weights - CIC scheme
          dx = x - real( jx )
          dy = y - real( jy )
          dz = z - real( jz )
          wt( 1, is ) = ( 1. - dx ) * ( 1. - dy ) * ( 1. - dz )
          wt( 2, is ) =        dx   * ( 1. - dy ) * ( 1. - dz )
          wt( 3, is ) = ( 1. - dx ) *        dy   * ( 1. - dz )
          wt( 4, is ) =        dx   *        dy   * ( 1. - dz )
          wt( 5, is ) = ( 1. - dx ) * ( 1. - dy ) *        dz
          wt( 6, is ) =        dx   * ( 1. - dy ) *        dz
          wt( 7, is ) = ( 1. - dx ) *        dy   *        dz
          wt( 8, is ) =        dx   *        dy   *        dz
        end if
      end do
      return
      end
