      subroutine grnorg
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to re-organise the greens function data, and (if present)
c        the  pop2  potentials from file  s3  to double buffered
c        format on files  s3, s6
c
      include 'rjinclude.h'
c
c local variables
      integer i, ir, is, is3, is6, i1, i2, j, jd, jr, ll, l1, l2, l23,
     +        l3, nrec
c
c     dummy addresses for initial 'last transfer'
      data is3, is6 / 1, 1 /
C
c recover greens function
      rewind s3
      read( s3 )l1, l2, l3, w( 1 ), w( 2 ), w( 3 ), i
      l23 = l2 * l3
      jd = 0
      do i = 1, l1
        read( s3 )( w( j + jd ), j = 1, l23 )
        jd = jd + l23
      end do
c 3 pairs of boundary planes
      i1 = l2
      i2 = l3
      do i = 1, 3
        if( i .eq. 3 )i2 = l2
        ll = ( i1 - 2 ) * i2
        read( s3 )( w( jd + j ), j = 1, ll )
        jd = jd + ll
        read( s3 )( w( jd + j ), j = 1, ll )
        jd = jd + ll
        i1 = l1
      end do
c calculate number of values in greens function and decide destination
      ir = n1 * n2 * n3
      jr = ir
      ir = ir + 2 * ( n1 * ( n2 + n3 ) + n2 * n3 )
      grnmem = ir .le. lngrnb
      if( grnmem )then
        write( s2, '( '' Storing Greens function in memory'' )' )
c release files s3 and s6
        call releas( s3 )
        call releas( s6 )
c copy greens function to storage area
        do  i = 1, ir
          grnbuf( i ) = w( i )
        end do
        return
      end if
c lose and re-create files s3, s6 to receive Greens function
      write( s2, 200 )s3, s6
 200  format( ' Storing Greens function on files', i10,
     1        ' and', i10 )
      i = bufer( 2, 1 ) - bufer( 1, 1 ) + 1
      i = ( skip + i - 1 ) / i
      j = i / 2
      i = nspio * max0( j, i - j )
      call swplsu( s3, 'tempft03', i )
      call swplsu( s6, 'tempft06', i )
      fil3 = .true.
      nrec3 = 0
      nrec6 = 0
      ll = bufer( 2, 1 )- bufer( 1, 1 ) + 1
      nrec = jd/ll
c write complete records to files
      is = 1
      do i = 1, nrec
        if( fil3 )then
          call dbwrit( s3, is, nrec3, ll, is3 )
          is3 = is
        else
          call dbwrit( s6, is, nrec6, ll, is6 )
          is6 = is
        end if
        fil3 = .not. fil3
        is = is + ll
      end do
      rewind s3
      rewind s6
      return
      end
