      subroutine difpot( iplane )
c written by Jerry Sellwood for the galaxy simulation code
c
c Tabulates the three Cartesian acceleration components on two adjacent planes
c   of grid points
c
c The acceleration components are determined by differencing the pre-calculated
c   potential array.  In order to avoid storing each of these components at
c   every grid point, which would require three full grids of extra memory, the
c   potential differences are determined for just two planes at a time.  (Two
c   planes are needed for linear interpolation.)  This strategy requires all
c   particles in one grid plane to be accelerated before any others are
c   processed, and it is this requirement that forces the linked list structure.
c
c The 10 point difference operator, devised by Andrew May and Richard James,
c   reduces force anisotropies well below those yielded by the obvious 2-point
c   formula.  [See Appendix A of Sellwood & Merritt (1994, ApJ v425, p530) for
c   a study of force quality.]
c
c calling argument
      integer iplane
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer i, j, k, l, m, n
      real f
c
      l = ipt( 2 ) + iplane * ngxy + ngx - 1
      n = ngx - 1
c work over plane
      do j = 2, ngy - 1
        l = l + 2
        n = n + 2
        do i = 2, ngx - 1
          l = l + 1
          m = l + ngxy
          n = n + 1
c x force components
          f = grids( l  + ngx + 1 ) - grids( l  + ngx - 1 ) +
     +        grids( l  - ngx + 1 ) - grids( l  - ngx - 1 ) +
     +        grids( l + ngxy + 1 ) - grids( l + ngxy - 1 ) +
     +        grids( l - ngxy + 1 ) - grids( l - ngxy - 1 ) +
     + 2. * ( grids( l        + 1 ) - grids( l        - 1 ) )
          k = ipt( 3 ) + n
          grids( k ) = f / ( 12. * dh( 3 ) )
          f = grids( m  + ngx + 1 ) - grids( m  + ngx - 1 ) +
     +        grids( m  - ngx + 1 ) - grids( m  - ngx - 1 ) +
     +        grids( m + ngxy + 1 ) - grids( m + ngxy - 1 ) +
     +        grids( m - ngxy + 1 ) - grids( m - ngxy - 1 ) +
     + 2. * ( grids( m        + 1 ) - grids( m        - 1 ) )
          k = ipt( 6 ) + n
          grids( k ) = f / ( 12. * dh( 3 ) )
c y force components
          f = grids( l    + 1 + ngx ) - grids( l    + 1 - ngx ) +
     +        grids( l    - 1 + ngx ) - grids( l    - 1 - ngx ) +
     +        grids( l + ngxy + ngx ) - grids( l + ngxy - ngx ) +
     +        grids( l - ngxy + ngx ) - grids( l - ngxy - ngx ) +
     + 2. * ( grids( l        + ngx ) - grids( l        - ngx ) )
          k = ipt( 4 ) + n
          grids( k ) = f / ( 12. * dh( 2 ) )
          f = grids( m    + 1 + ngx ) - grids( m    + 1 - ngx ) +
     +        grids( m    - 1 + ngx ) - grids( m    - 1 - ngx ) +
     +        grids( m + ngxy + ngx ) - grids( m + ngxy - ngx ) +
     +        grids( m - ngxy + ngx ) - grids( m - ngxy - ngx ) +
     + 2. * ( grids( m        + ngx ) - grids( m        - ngx ) )
          k = ipt( 7 ) + n
          grids( k ) = f / ( 12. * dh( 2 ) )
c z force components
          f = grids( l   + 1 + ngxy ) - grids( l   + 1 - ngxy ) +
     +        grids( l   - 1 + ngxy ) - grids( l   - 1 - ngxy ) +
     +        grids( l + ngx + ngxy ) - grids( l + ngx - ngxy ) +
     +        grids( l - ngx + ngxy ) - grids( l - ngx - ngxy ) +
     + 2. * ( grids( l       + ngxy ) - grids( l       - ngxy ) )
          k = ipt( 5 ) + n
          grids( k ) = f / ( 12. * dh( 1 ) )
          f = grids( m   + 1 + ngxy ) - grids( m   + 1 - ngxy ) +
     +        grids( m   - 1 + ngxy ) - grids( m   - 1 - ngxy ) +
     +        grids( m + ngx + ngxy ) - grids( m + ngx - ngxy ) +
     +        grids( m - ngx + ngxy ) - grids( m - ngx - ngxy ) +
     + 2. * ( grids( m       + ngxy ) - grids( m       - ngxy ) )
          k = ipt( 8 ) + n
          grids( k ) = f / ( 12. * dh( 1 ) )
        end do
      end do
      return
      end
