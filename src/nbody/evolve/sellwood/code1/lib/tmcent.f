      subroutine tmcent( start )
c written by Jerry Sellwood for the galaxy simulation code
c
c Sets, or undoes, a half time step offset between the positions and velocities
c
c Steps velocities only half a step backwards (forwards) in order to set up
c   (undo) time centered velocities for leap-frog type integration schemes.
c   The positions are unaffected.
c
c If the logical variable start is .TRUE. then the velocities are backed
c   up half a step.
c This routine should be preceded by calls to MASSET and FINDF in order that
c   the correct acceleration components can be determined in the usual way.
c The routine process all the particles.
c
c calling argument
      logical start
c
c common block
c
      include 'admin.h'
c
c local variables
      integer i, jst
c
      if( start )then
        print *, 'Time-centering velocities'
        write( no, * )'Time-centering velocities'
      else
        print *, 'Un-centering velocities'
        write( no, * )'Un-centering velocities'
      end if
c initialize new linked list table
      do i = 1, nlists - 1
        islist( 2, i ) = -1
      end do
c work through all lists of particles except those off the grid
      do ilist = 1, nlists - 1
        izone = ilist - 1
        inext = islist( 1, ilist )
        if( inext .ge. 0 )call difpot( izone )
c work through groups of particles
        do while ( inext .ge. 0 )
          call getgrp( jst )
c get acceleration components
          call getacc( jst )
c time centre velocities
          call tmcgrp( jst, start )
c restore coordinates
          call stogrp( jst )
        end do
      end do
c save new list origins
      do i = 1, nlists - 1
        islist( 1, i ) = islist( 2, i )
      end do
      return
      end
