      subroutine loadup
c written by Jerry Sellwood for the galaxy simulation code
c
c Reads, rescales and stores the initial particle positions and velocities
c
c The data are read from the file galaxy.ini, converted into internal program
c   units and stored in memory.  A separate linked list is constructed for
c   particles lying between each pair of planes in order to enable them to
c   be found efficiently when applying accelerations.
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
c initialise linked list origins
      nx = 0
      do i = 1, nlists
        islist( 2, i ) = -1
      end do
c open file and read header record
      ndistf = 10
      call opnfil( ndistf, 'ini', 'formatted', 'old', 'seq', i )
c
c a = time now
c tmass is the total mass of the model (G=1)
c nbod is the number of particles
c
      read( ndistf, * )a, tmass, nbod
c check / ptcls /  space
      lpf = nwpp * nbod
	write(*,*) nbod,nwpp,lpf,lptcls
      if( lpf .gt. lptcls )call crash( 'LOADUP',
     +                                 'Particle array too small' )
c set step number
      istep = nint( a / ts )
c process particles from file
      is = 1
      do j = 1, nbod
c get another particle
        read( ndistf, * )x, v
c clear buffer when full
        if( is .gt. mbuff )then
          call stogrp( mbuff )
          is = 1
        end if
c convert to internal units and store
        do i = 1, 3
          newc( i, is ) = lscale * x( i )
          newc( i + 3, is ) = lscale * ts * v( i )
        end do
        iz( is ) = newc( 3, is ) + zm + 1.
c mark particles outside the grid 
        if( ( abs( newc( 1, is ) ) .gt. xm ) .or.
     +      ( abs( newc( 2, is ) ) .gt. ym ) .or.
     +      ( abs( newc( 3, is ) ) .gt. zm ) )then
          iz( is ) = nlists
          noff = noff + 1
        end if
        loc( is ) = nx
        nx = nx + nwpp
        is = is + 1
      end do
      close( ndistf )
c clear remaining particles from buffer
      is = is - 1
      call stogrp( is )
c
      if( noff .gt. 0 )then
        write( no, * )noff, ' particles outside the grid at the start'
        print *, noff, ' particles outside the grid at the start'
      end if
c preserve linked list origins
      do i = 1, nlists
        islist( 1, i ) = islist( 2, i )
      end do
c set mass of each particle
      write( no,
     +       '( '' Total mass of all particles is'', f10.6 )' )tmass
      pmass = tmass * lscale**3 * ts**2 / real( nbod )
      return
      end
