      subroutine unload
c written by Jerry Sellwood for the galaxy simulation code
c
c outputs the particles to the galaxy.fin file
c
c Converts the current particle positions and velocities to the original
c   supplied units and writes them to the file in exactly the order and
c   format in which they were read in.  Those particles that are outside the
c   grid volume will have remained motionionless from the moment they left!
c
c common blocks
c
      include 'admin.h'
c
c local arrays
       real x( 3 ), v( 3 )
c
c local variables
       integer i, is, izon, j, nx
       real a
c
c open file and write header record
      ndistf = 10
      call opnfil( ndistf, 'fin', 'formatted', 'unknown', 'seq', i )
      a = ts * real( istep )
      write( ndistf, * )a, tmass, nbod
c process particles
      do j = 1, lpf, nwpp
        if( iptcls( j + ncoor ) .lt. 0 )then
c particle off the grid
          do i = 1, 3
            x( i ) = 0
            v( i ) = 0
          end do
        else
c rescale coordinates
          do i = 1, 3
            x( i ) = ptcls( j + i - 1 ) / lscale
            v( i ) = ptcls( j + i + 2 ) / ( lscale * ts )
          end do
        end if
        write( ndistf, '( 6e14.6 )' )x, v
      end do
      close( ndistf )
      return
      end
