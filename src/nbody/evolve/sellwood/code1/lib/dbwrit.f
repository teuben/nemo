      subroutine dbwrit(ich, is, nrec, n, isl)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to check the current output transfer on a file and to
c        initiate a new one.
c
c     parameters
c
c     ich   =  channel number
c
c     is    -  a pointer to the beginning of information in /scm/.
c
c     nrec  =  the number of the record to be written to the file.
c              the routine updates this number after the transfer.
c
c     n     =  the transfer length
c
c     isl   -  the identifier for the transfer in progress on entry.
c
      include 'rjinclude.h'
c
c local variables
      integer ich, iii, is, isl, n, nrec, ntran
c
c     calculate block count, check transfer length
c
      ntran = n/512
      if((n - 512*ntran).ne.0) then
        write(s2, 200 )n
 200    format( ' transfer length', i10, 5x,
     1          'is not an integral number of blocks' )
        write(s2, '('' run stopped by subroutine dbwrit(wn)'')')
        call crash( 'dbwrit', '1' )
      end if
c
c     initiate new write and finish
c
 1    write( ich ) ( w( iii ), iii = is, is + n - 1 )
c
      nrec = nrec + 1
      return
      end
