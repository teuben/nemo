      subroutine findf
c written by Jerry Sellwood for the galaxy simulation code
c
c Calculates the gravitational potential on a 3-D Cartesian grid
c
c It simply copies the pre-set density array to the target area, calls
c   Richard James's Poisson solver and copies the resulting potential array
c   back to the area expected by the N-body code
c
c This version is rather wasteful of memory since Richard's Poisson solver
c   uses arrays that are totally independent of those in the N-body code
c
c common blocks
c
      include 'admin.h'
c
      include 'rjinclude.h'
c
c local variable
      integer i
c
c copy mass distribution to / rjscm / and rescale, including extra fudge factor
      do i = 1, mesh
        w( scmbuf - 1 + i ) = grids( ipt( 1 ) + i ) * c3pmf * pmass
      end do
c solve for the potential
      call poiss
c copy potential array back to / grids /
      do i = 1, mesh
        grids( ipt( 2 ) + i ) = w( scmbuf - 1 + i )
      end do
      return
      end
