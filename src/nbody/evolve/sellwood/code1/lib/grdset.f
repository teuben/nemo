      subroutine grdset
c written by Jerry Sellwood for the galaxy simulation code
c
c Reads parameters for 3-D Cartesian grid from galaxy.dat file
c
c The input parameters required are the grid size, two scale variables and
c   the time between outputs.
c The scale variables are lscale and ts: lscale is the number of mesh spaces
c   per unit of distance in the particle input file and ts is the time step.
c The time between outputs is converted to an integral number of steps.
c
c It also sets up the arcane input data file for Richard James's Poisson solver
c
c common block
c
      include 'admin.h'
c
c local arrays
      integer ichan( 9 ), ngsize( 3  )
      real*8 wgt( 3 )
c
c local variables
      character line*80
      integer i, j, k, naux
      real a
      real*8 den
c
      data ichan / 7, 8, 20, 21, 1, 23, 9, 10, 11 /
c
c logical unit nos for input/output
      ni = 7
      no = 8
      call opnfil( ni, 'dat', 'formatted', 'old', 'seq', i )
      if( i .ne. 0 )call crash( 'GRDSET', '.dat file not found' )
c open output file
      call opnfil( no, 'lis', 'formatted', 'unknown', 'append', i )
      if( i .ne. 0 )call crash( 'GRDSET', 'Error opening .lis file' )
c read in mesh size
      call getline( ni, line )
      read( line, *, err = 1 )ngx, ngy, ngz
c convert grid dimensions to Richard's input format
      ngsize( 1 ) = ngz
      ngsize( 2 ) = ngy
      ngsize( 3 ) = ngx
      do i = 1, 3
        j = 0
        k = 1
        do while ( k .lt. ngsize( i ) )
          j = j + 1
          k = 2**j + 1
        end do
        if( k .ne. ngsize( i ) )then
          print *, ngsize
          call crash( 'GRDSET', 'Invalid mesh dimensions' )
        end if
        ngsize( i ) = j
      end do
c mesh cell shapes
      do i = 1, 3
        dh( i ) = 1
      end do
c convert shape parameters to Richard's format
      a = max( dh( 1 ), dh( 2 ), dh( 3 ) )
      den = 2. * ( dh( 1 )**(-2) + dh( 2 )**(-2) + dh( 3 )**(-2) )
      do i = 1, 3
        wgt( i ) = dh( i )**(-2) / den
        dh( i ) = dh( i ) / a
      end do
c create auxiliary input file for Richard's software
      naux = 4
      call opnfil( naux, 'aux', 'formatted', 'unknown', 'seq', i )
      i = 5
      write( naux, '( ''mesh      '', 3i10 )' )i, i, i
      write( naux, '( ''channels'' )' )
      write( naux, '( 9i6 )' )ichan
      write( naux, '( ''params'' )' )
      write( naux, '( 3f20.11 )' )wgt
      write( naux, '( ''green'' )' )
      write( naux, '( ''end'' )' )
c
      write( naux, '( ''mesh      '', 3i10 )' )ngsize
      write( naux, '( ''channels'' )' )
      write( naux, '( 9i6 )' )ichan
      write( naux, '( ''params'' )' )
      write( naux, '( 3f20.11 )' )wgt
      write( naux, '( ''green'' )' )
      write( naux, '( ''end'' )' )
c
      write( naux, '( ''mesh      '', 3i10 )' )ngsize
      write( naux, '( ''channels'' )' )
      write( naux, '( 9i6 )' )ichan
      write( naux, '( ''params'' )' )
      write( naux, '( 3f20.11 )' )wgt
      write( naux, '( ''poisson'' )' )
      write( naux, '( ''source'' )' )
      write( naux, '( ''potential'' )' )
      write( naux, '( ''end'' )' )
c
      write( naux, '( ''end'' )' )
      close( naux )
c grid size
      write( no, '( //'' Mesh size'', 3i6 )' )ngx, ngy, ngz
      ngxy = ngx * ngy
      mesh = ngx * ngy * ngz
      nplanes = ngz - 3
c coordinate boundaries in mesh spaces - active mesh is smaller because
c   accelerations have to be derived from potential differences
      xm = .5 * real( ngx - 3 )
      ym = .5 * real( ngy - 3 )
      zm = .5 * real( ngz - 3 )
c separate list for each plane plus one for hot spot and one for particles off
      nlists = nplanes + 2
      noff = 0
      ncoor = 6
      nwpp = ncoor + 2
c set pointers for / scm /
      ipt( 1 ) = 0
      ipt( 2 ) = ipt( 1 ) + mesh
      j = ipt( 2 ) + mesh
      do i = 1, 6
        ipt( i + 2 ) = j
        j = j + ngxy
      end do
c number of grid spaces per length unit
      call getline( ni, line )
      read( line, *, err = 1 )lscale
      write( no, * )'lscale value read', lscale
      if( lscale .le. 0. )call
     +                crash( 'GRDSET', 'Value of lscale nonsensical' )
c length of time step
      call getline( ni, line )
      read( line, *, err = 1 )ts
      write( no, * )'time step value read', ts
      if( ts .le. 0. )call crash( 'GRDSET', 'Value of ts nonsensical' )
c time between particle outputs
      call getline( ni, line )
      read( line, *, err = 1 )a
      if( a .le. 0. )then
        outstp = 0
      else
        outstp = nint( a / ts )
      end if
      write( no, * )'Particles will be output every',
     +               outstp, ' time steps'
c time between integral checks
      call getline( ni, line )
      read( line, *, err = 1 )a
      if( a .le. 0. )then
        chkstp = 0
      else
        chkstp = nint( a / ts )
      end if
      write( no, * )'Integrals will be checked every',
     +               chkstp, ' time steps'
c logical unit number for .res file
      nres = 9
      return
c error reading line
    1 write( no, '( a )' )line( 1:60 )
      call crash( 'GRDSET', 'Error reading input data line' )
      stop
      end
