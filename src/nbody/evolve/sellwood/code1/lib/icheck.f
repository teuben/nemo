      subroutine icheck( jst )
c written by Jerry Sellwood for the galaxy simulation code
c
c Adds contributions from current group of particles to the global integrals
c   working in internal program units
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
      real v( 3 )
c
c work through all particles in this group
      do is = 1, jst
        ncheck = ncheck + 1
        do i = 1, 3
          com( i ) = com( i ) + oldc( i, is )
c time centered velocity components
          v( i ) = .5 * ( oldc( i + 3, is ) + newc( i + 3, is ) )
          lmom( i ) = lmom( i ) + v( i )
          ke = ke + v( i )**2
          claus = claus + oldc( i, is ) * acc( i, is )
        end do
c the half here is because the summation includes every pair of particles twice
        pe = pe + .5 * gpot( is )
c angular momentum components
        amom( 1 ) = amom( 1 ) + v( 2 ) * oldc( 1, is ) -
     +                          v( 1 ) * oldc( 2, is )
        amom( 2 ) = amom( 2 ) + v( 1 ) * oldc( 3, is ) -
     +                          v( 3 ) * oldc( 1, is )
        amom( 3 ) = amom( 3 ) + v( 3 ) * oldc( 2, is ) -
     +                          v( 2 ) * oldc( 3, is )
      end do
      return
      end
