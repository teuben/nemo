      subroutine masscl
c written by Jerry Sellwood for the galaxy simulation code
c
c Resets the mass array to zero
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer j, k, l
c
c set pointers
      k = ipt( 1 ) + 1
      l = ipt( 1 ) + mesh
c set mass array to zero
      do j = k, l
        grids( j ) = 0
      end do
      return
      end
