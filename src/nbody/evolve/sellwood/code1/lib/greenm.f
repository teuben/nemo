      subroutine greenm
c written by Jerry Sellwood for the galaxy simulation code but based on
c   Richard James's own start up program
c
c Creates the Greens function for Richard James's Poisson solver
c
c Input data in the appropriate format must be already set up in the file
c    galaxy.aux by a previous call to grdset
c
c common blocks
c
      include 'rjinclude.h'
c
      include 'admin.h'
c
c local variables
      character instr*6
      integer i, il, i1, j, k, lenbuf
      parameter ( lenbuf = 8192 )
      real u
c
c open auxiliary input file for Richard's software
      close( ni )
      call opnfil( ni, 'aux', 'formatted', 'old', 'seq', i )
      if( i .ne. 0 )call crash( 'GREENM',
     +                          '.aux file apparently not available' )
c open temporary output file to discard output file for Richard's software
      close( no )
      call opnfil( no, 'tmp', 'formatted', 'unknown', 'seq', i )
      if( i .ne. 0 )call crash( 'GREENM', 'Error opening .tmp file' )
c
      k = 0
      do i = 1, 20
        bufer( 1, i ) = k + 1
        k = k + lenbuf
        bufer( 2, i ) = k
      end do
c
c set up control parameters for potential solver
c
      call psnset( lev3, lnfac, 1 )
c
c skip Greens function calculation if not required
c
      if( gren )then
        call greenf( .false. )
        call poiss
c
        istack = 0
        call greenl( lnfac, .false. )
        call psnset( lev3, lnfac, 1 )
      end if
c
c     adjust point selection array.  negative values imply the top,
c     north or west plane as appropriate
c
      do i = 1, nverif
        if( iverif( 1, i ) .lt. 0 )iverif( 1, i ) = n1 - 1
        if( iverif( 2, i ) .lt. 0 )iverif( 2, i ) = n2 - 1
        if( iverif( 3, i ) .lt. 0 )iverif( 3, i ) = n3 - 1
        iverif( 1, i ) = max( iverif( 1, i ), 0 )
        iverif( 2, i ) = max( iverif( 2, i ), 0 )
        iverif( 3, i ) = max( iverif( 3, i ), 0 )
        iverif( 1, i ) = min( iverif( 1, i ), n1 - 1 )
        iverif( 2, i ) = min( iverif( 2, i ), n2 - 1 )
        iverif( 3, i ) = min( iverif( 3, i ), n3 - 1 )
      end do
c
c If verification requested for potential solver, copy direct
c   Greens function file to channel 4.  This procedure avoids
c   opening the channel if verification is not required.
c
      read( ni, '( a6 )' )instr
c      if( instr( 1:3 ) .ne. 'end' )then
c
c        rewind s6
c        call scfilu( s4, 'tempft04', skip / 512 + 5 )
c        rewind s4
c        read( s6 )i, j, k, ( w( i ), i = 1, 3 ), i1
c        write( s4 )j, k, i, ( w( i ), i = 1, 3 ), i1
cc
cc read and transpose direct Greens function.
cc
c        do k = 1, n3
c          do i = 1, n1
c            il = ( i - 1 ) * n2 * n3 + k - n3
c            read( s6 )( w( il + j * n3 ), j = 1, n2 )
c          end do
c        end do
c        do i = 1, n1 * n2
c          il = ( i - 1 ) * n3
c          write( s4 )( w( il + j ), j = 1, n3 )
c        end do
c        print *, 'Direct Greens function transposed'
c      end if
c
c re-organise Greens function
c
      call grnorg
      print *, 'Greens function creation complete'
c
c check poisson solver if called for
c
c      do while ( instr( 1:3 ) .ne. 'end' )
c        if( instr .eq. 'chpois' )then
c          print *, 'Calling chpois'
c          call chpois
c        end if
c        if( instr .eq. 'challs' )then
c          print *, 'Calling challs'
c          call challs
c        end if
c        if( instr .eq. 'challp' )then
c          print *, 'Calling challp'
c          call challp
c        end if
c        read( ni, '( a6 )' )instr
c      end do
      close( ni, status = 'delete' )
c extra factor for mass normalization
      u = dh( 1 )**(-2) + dh( 2 )**(-2) + dh( 3 )**(-2)
      u = sqrt( 2. * u )
      if( scfact .eq. 0. )then
        c3pmf = u / factor
      else
        c3pmf = u / scfact
      end if
c close and discard scratch output file and re-open main output file'
      close( no, status = 'delete' )
      call opnfil( no, 'lis', 'formatted', 'unknown', 'append', i )
      if( i .ne. 0 )call crash( 'GREENM', 'Error opening .lis file' )
      write( no, * )'greenm completed'
      return
      end
