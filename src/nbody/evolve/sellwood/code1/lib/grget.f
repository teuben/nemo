      subroutine grget(x, n)
c routine written by Richard James for 3-D Poisson solver
      integer n
      real x( n )
c
c     purpose
c
c        to recover the next plane of the greens function or pop2
c        potential from buffer, and if necessary to organise a new
c        file transfer
c
c     parameters
c
c     x  -  a real array to receive the data.
c
c     n  =  the number of words of data to be transferred.
c
      include 'rjinclude.h'
c
c local variables
      integer i, id, liend, il, nn, nnn
c
c      real x(n)
c
c     if available from memory, use memory version.
c
      if(grnmem) then
        do 8 i = 1, n
 8      x(i) = grnbuf(iptgrn + i)
        iptgrn = iptgrn + n
        return
      end if
c
c     check data availability
c
      nn = n
      id = 0
      if(iptgrn.gt.limgrn) go to 3
c
c     transfer available data, return if complete
c
 1    liend = min0(limgrn, iptgrn + nn - 1)
      do 2 i = iptgrn, liend
 2    x(id + i - iptgrn + 1) = area(i)
      nnn = liend - iptgrn + 1
      nn = nn - nnn
      id = id + nnn
      iptgrn = liend + 1
      if(nn.eq.0) return
c
c     find buffer if on file s6
c
 3    if(fil3) go to 6
      irec6 = irec6 + 1
      if(irec6.le.nrec6) go to 4
c
c     searching for a non-existent record
c
      write(s2, 200) irec6, s6, nrec6
200   format('0seeking record', i10, 5x, 'on channel', i10,
     1 5x, 'limit =', i10)
      stop
c
c     check file
c
 4    i = bufer(1, icr6)
      iptgrn = bufer(1, icr6)
      limgrn = bufer(2, icr6)
      fil3 = .true.
c
c     swap buffers
c
      i = icr6
      icr6 = inx6
      inx6 = i
c
c     initiate new transfer if not finished
c
      if(irec6.eq.nrec6) go to 1
      il = bufer(1, icr6)
      read( s6 )( area( i ), i = il, bufer( 2, icr6 ) )
      go to 1
c
c     find buffer on file s3
c
 6    irec3 = irec3 + 1
      if(irec3.le.nrec3) go to 7
c
c     searching for a non-existent record
c
      write(s2, 200) irec3, s3, nrec3
      stop
c
c     check file
c
 7    i = bufer(1, icr3)
      iptgrn = bufer(1, icr3)
      limgrn = bufer(2, icr3)
      fil3 = .false.
c
c     swap buffers
c
      i = icr3
      icr3 = inx3
      inx3 = i
c
c     initiate new transfer if not finished
c
      if(irec3.eq.nrec3) go to 1
      il = bufer(1, icr3)
      read( s3 )( area( i ), i = il, bufer( 2, icr3 ) )
      go to 1
      end
