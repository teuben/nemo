      subroutine putout
c written by Jerry Sellwood for the galaxy simulation code
c
c Appends the current particle coordinates to the galaxy.res file
c
c Appends the pre-prepared array of particles to the .res file in exactly the
c   order in which they were read in.  All 6 phase space coordinates of a
c   particle outside the grid volume will be zero.  In order to ensure this
c   the entire array is reset to zero after it is written.  The actual values
c   to be output are prepared in subroutine CENVEL
c
c common blocks
c
      include 'admin.h'
c
c local variables
       integer i
       real t
c
c (re)-open the file
      call opnfil( nres, 'res', 'unformatted', 'unknown', 'append', i )
      if( i .ne. 0 )call crash( 'PUTOUT', 'Error opening .res file' )
c output the current time and the entire array
      t = ts * real( istep )
      write( nres )t, nbod
      write( nres )( coor( i ), i = 1, 7 * nbod )
c close the file after every write for safety (in case of a system crash)
      close ( nres )
      write( no, * )'Particles output to .res file at step and time',
     +              istep, t
c reset the array to zero
      do i = 1, 6 * nbod
        coor( i ) = 0
      end do
      return
      end
