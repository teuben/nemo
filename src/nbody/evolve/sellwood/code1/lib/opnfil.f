      subroutine opnfil( iunit, typ, frm, stat, acc, istat )
c written by Jerry Sellwood for the galaxy simulation code
c
c Opens a new file with name galaxy.typ having the specified attributes
c
c calling arguments
      character typ*(*), frm*(*), stat*(*), acc*(*)
      integer istat, iunit
c
c set default return flag
      istat = 0
      if( acc .eq. 'append' )then
c open in append mode
        open( iunit, err = 1, file = 'galaxy.'//typ(1:3),
     +        status = stat, form = frm, access = acc )
      else
c standard sequential open
        open( iunit, err = 1, file = 'galaxy.'//typ(1:3),
     +        status = stat, form = frm )
      end if
      return
c
c error occured - print a warning message and return an error flag
    1 istat = 1
      print *, 'Failed to open file galaxy.',typ, ' on unit', iunit
      print *, '         with status = ', stat, ', form = ', frm,
     +         ', access = ', acc
      return
      end
