      subroutine newfil
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to attach, create and close files as required using
c        ansi fortran statements, but leaving the possibility of
c        incorporating direct system calls if required.  actions taken
c        are reported.  multiple entry points are used.
c
      include 'rjinclude.h'
c
      include 'admin.h'
c
c local variables
      integer i, ich, istat, llen
c
      character*8 fil, messge*80
      logical advanc, form
c
c     define and attach unformatted file
c
      entry definu( ich, fil, llen )
      write( messge, 300 )fil, ich
 300  format( 'define unformatted file ', a8, ' on channel', i10 )
c      write( no, '( 1x, a80 )' )messge
    1 continue
      open( ich, file = fil, form = 'unformatted', status = 'unknown',
     +     iostat = istat )
      if( istat .ne. 0 )then
        write( no, 200 )ich, fil, istat
 200    format( ' Error in opening file number', i10,
     +          ' name = ', a8, ' iostat =', i10 )
        call crash( 'newfil', 'wq 1' )
        stop
      end if
      fnamtb( ich ) = fil
      return
c
c     define and attach formatted file
c
      entry definf( ich, fil, llen )
      write( messge, 301 )fil, ich
 301  format( 'define formatted file ', a8, ' on channel', i10 )
c      write( no, '( 1x, a80 )' )messge
    2 continue
      open( ich, file = fil, status = 'old', iostat = istat )
      if( istat.ne.0 ) then
        write( no, 200 )ich, fil, istat
        call crash( 'newfil', 'wq 2' )
        stop
      end if
      fnamtb( ich ) = fil
      return
c
c     swap old file for new
c
      entry swpkpu( ich, fil, llen, advanc )
      form = .false.
      go to 3
      entry swpkpf( ich, fil, llen, advanc )
      form = .true.
    3 continue
      close( ich )
      write( messge, 302 )ich
 302  format( 'file number', i10, '  released' )
c      write( no, '( 1x, a80 )' )messge
      if( advanc ) then
        read( fil( 7:8 ), '( i2 )' ) i
        i = mod( i + 1, 100 )
        write( fil( 7:8 ), '( i2 )' ) i
        if( fil( 7:7 ).eq.' ' ) fil( 7:7 ) = '0'
      end if
      if( form ) then
        write( messge, 303 )ich, fil
 303    format( 'new formatted file on channel', i10, ' has name ', a8 )
      else
        write( messge, 304 )ich, fil
 304    format( 'new unformatted file on channel', i10, ' has name ',
     +           a8 )
      end if
c      write( no, '( 1x, a80 )' )messge
      if( form ) go to 2
      go to 1
c
c     open unformatted scratch file
c
      entry scfilu( ich, fil, llen )
      write( messge, 305 )fil, ich
 305  format( 'open unformatted scratch file ', a8, ' on channel', i10 )
c      write( no, '( 1x, a80 )' )messge
      if( ( fil( 8:8 ).eq.'7' ).and.( fil( 7:7 ).ne.'0' ) )call crash(
     +    'NEWFIL', 'Exit placed by Richard' )
    4 open( ich, file = fil, form = 'unformatted', status = 'unknown',
     +     iostat = istat )
      if( istat.ne.0 ) then
        write( no, 200 )ich, fil, istat
        call crash( 'newfil', 'wq 3' )
        stop
      end if
      fnamtb( ich ) = fil
      return
c
c     attach and open existing file for unformatted read
c     access only
c
      entry openfu( ich, fil )
      write( messge, 306 )fil, ich
 306  format( 'open unformatted read only file ', a8, ' on channel',
     +        i10 )
c      write( no, '( 1x, a80 )' )messge
      open( ich, file = fil, form = 'unformatted', status = 'old',
     +     iostat = istat )
      if( istat.ne.0 ) then
        write( no, 200 )ich, fil, istat
        call crash( 'newfil', 'wq 4' )
        stop
      end if
      fnamtb( ich ) = fil
      return
c
c     close current file and re-open as an unformatted scratch file
c
      entry swplsu( ich, fil, llen )
      write( messge, 307 )fil, ich
 307  format( 'reopen file ', a8,
     +        ' as unformatted scratch file on channel', i10 )
c      write( no, '( 1x, a80 )' )messge
      close( ich )
      go to 4
c
c     release current file.
c
      entry releas( ich )
      write( messge, 308 )fnamtb( ich ), ich
 308  format( 'release file ', a8, ' channel', i10 )
c      write( no, '( 1x, a80 )' )messge
      close( ich, status = 'delete' )
      return
      end
