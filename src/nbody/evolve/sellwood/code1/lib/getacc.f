      subroutine getacc( jst )
c written by Jerry Sellwood for the galaxy simulation code
c
c Returns acceleration acting on each particle in the current group arising
c   from the mass distribution on the grid.
c
c The acceleration components are the weighted sums of the values the eight
c   nearest cell corners.  The appropriate weights are obtained by a call to
c   subroutine WEIGHT, which assumes that the "old" coordinates contain the
c   appropriate position.
c
c Linear interpolation between mesh points is hard-wired.
c
c No supplementary or external components are added in this routine.
c
c calling argument
      integer jst
c
c common blocks
c
      include 'admin.h'
c
c local array
      integer ic( 8, 3 )
c
c local variables
      integer i, is, j, k
c
c compute grid cell number and weights for interpolation
      call weight( jst, .true. )
c re-initialize
      do is = 1, jst
        do i = 1, 3
          acc( i, is ) = 0
        end do
        gpot( is ) = 0
      end do
c
      do is = 1, jst
        if( nskip( is ) )then
c tabulate locations of cell corners in each array
          do i = 1, 3
            ic( 1, i ) = ipt( i + 2 ) + ncl( is )
            ic( 2, i ) = ic( 1, i ) + 1
            ic( 3, i ) = ic( 1, i ) + ngx
            ic( 4, i ) = ic( 1, i ) + ngx + 1
            ic( 5, i ) = ipt( i + 5 ) + ncl( is )
            ic( 6, i ) = ic( 5, i ) + 1
            ic( 7, i ) = ic( 5, i ) + ngx
            ic( 8, i ) = ic( 5, i ) + ngx + 1
          end do
c evaluate accelerations
          do i = 1, 3
            do j = 1, 8
              k = ic( j, i )
              acc( i, is ) = acc( i, is ) + grids( k ) * wt( j, is )
            end do
          end do
c potential value
          do i = 1, 4
            j = ic( i, 1 ) - ipt( 3 ) + ipt( 2 ) + iz( is ) * ngxy
            gpot( is ) = gpot( is ) - grids( j ) * wt( i, is )
            k = j + ngxy
            gpot( is ) = gpot( is ) - grids( k ) * wt( i + 4, is )
          end do
        end if
      end do
      return
      end
