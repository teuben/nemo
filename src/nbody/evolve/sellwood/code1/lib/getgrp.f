      subroutine getgrp( jst )
c written by Jerry Sellwood for the galaxy simulation code
c
c Gathers the next group of particles belonging to the current list into the
c   arrays in / buffer /
c
c Each particle has the appropriate number of phase space coordinates, ncoor
c   (= 4 or 6), stored in floating point plus an additional two integers.  The
c   first may be used as a flag.  The second is the link to the next particle
c   in the list (or a negative value if the list ends with this particle)
c
c calling argument - input value ignored, returned value is the actual
c   number of particles gathered.  The maximum is set by the parameter
c   mbuff which is used to dimension the workspace arrays
      integer jst
c
c common blocks
c
      include 'admin.h'
c
c local variable
      integer i
c
      jst = 0
c gather particles from linked list
      do while ( jst .lt. mbuff )
        jst = jst + 1
        do i = 1, ncoor
          oldc( i, jst ) = ptcls( inext + i )
        end do
        iflag( jst ) = iptcls( inext + ncoor + 1 )
c save additional info
        ncl( jst ) = 0
        loc( jst ) = inext
        iz( jst ) = izone
c pick up next pointer
        inext = iptcls( inext + nwpp )
c check for end of list
        if( inext .lt. 0 )return
      end do
      return
      end
