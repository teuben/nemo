      subroutine crash( routine, reason )
c utility routine written by Jerry Sellwood
c
c Writes an error message to standard output and to unit number no and then
c   halts execution
c
c calling arguments
      character*(*) routine, reason
c
c common block
c
      include 'admin.h'
c
      print *, 'Run terminated'
      print '( 5x, a )', reason
      print 200, routine
  200 format( 'Error in ', a )
c
      write( no, * )'Run terminated'
      write( no, '( 5x, a )' )reason
      write( no, 200 )routine
      stop
      end
