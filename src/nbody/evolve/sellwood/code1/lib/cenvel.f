      subroutine cenvel( jst )
c written by Jerry Sellwood for the galaxy simulation code
c
c Prepares positions and time-centred velocities of the current group of
c   particles for output
c
c The time-centred velocity in simple 2nd-order leap-frog is the average of the
c   velocities before and after acceleration.  This routine must be called at
c   at the the group is accelerated because that is the only moment when the
c   velocities at both times and the potential are available.
c
c calling arguments
      integer jst
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer i, is, j
c
      do is = 1, jst
        j = 7 * ( loc( is ) / 8 )
c save rescaled position
        do i = 1, 3
          coor( j + i ) = oldc( i, is ) / lscale
        end do
c compute and save rescaled time-centered velocities
        do i = 4, 6
          coor( j + i ) = .5 * ( oldc( i, is ) + newc( i, is ) ) /
     +                                                  ( lscale * ts )
        end do
c the gravitational potential at the position of the particle in external units
        coor( j + 7 ) = gpot( is ) / ( lscale * ts )**2
      end do
      return
      end
