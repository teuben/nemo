      subroutine poiss
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to coordinate the calculation of the potential.  the routine
c        chooses the boundary calculation routine according to whether
c        it is calculating a regular potential or a greens function.
c        when calculating a regular potential, it sets up files
c        s3, s6 for double-buffered reading by grget.
c
c     note:
c
c        this routine uses the vector syntax facilities peculiar
c     to cyber fortran 200.  to run on a conventional machine,
c     it will be necessary to remove the descriptor statement and
c     the assign statements for  bdy.  the calls for subroutines
c     using vector facilities (vftsna, vftsns, vslvpt, vbdpot)
c     need not be removed, provided that  cal(1)  in  /caluse/
c     is set  .false.  the explicit vector syntax statements used
c     in this routine are delimited by a line of  *  symbols.
c
      include 'rjinclude.h'
      integer n(3)
      equivalence (n, n1)
c
c local variables
      integer i, il, in, ip, is, istck0, j, jp, k, nn, n12, n123, n13,
     1        n23
c
      logical par
      integer bdy(6, 2), stage, dir
c
c     set up pointers, advance and check stack pointer.
c
      j = istack
      bdy(1, 2) = nom(3)
      bdy(2, 2) = nom(3)
      bdy(3, 2) = nom(1)
      bdy(4, 2) = nom(1)
      bdy(5, 2) = nom(2)
      bdy(6, 2) = nom(2)
      do 16 i = 1, 6
      bdy(i, 1) = j
16    j = j + bdy(i, 2)
      istck0 = istack
      istack = j
      if(istack.gt.maxstk) then
        maxstk = istack
        mxstid = 'poiss '
      end if
      if(istack.gt.lstack) then
        write(s2, 200 )istck0, istack, lstack
 200    format( ' stack overflow in poiss (tg)' /
     1          ' old value =', i10, '  new value =', i10,
     2          ' limit =', i10 )
        write(s2, '('' boundary pointer array'', 2(/1x, 6i10))') bdy
        call crash( 'poiss', 'tg 1')
        stop
      end if
c
c
c     skip file initialisation if calculating greens function
c
      if(gren) go to 2
c
c     if greens function from memory, set pointer and skip file
c     initialisation.
c
      if(grnmem) then
        iptgrn = 0
        go to 2
      end if
c
c     initiate greens function transfers and set up pointers
c
c      include 'inc/rjzztg1.f'
c
      rewind s3
c
      il = bufer(1, 1)
c      include 'inc/rjzztg2.f'
c      read( s3, id = il ) area( il )...area( bufer( 2, 1 ) )
c
      read( s3 )( area( i ), i = il, bufer( 2, 1 ) )
c
      icr3 = 1
      inx3 = 2
c
c     we ensure that the greens function occupies more than
c     one buffer before initiating the second transfer
c
      i = (n1 - 2)*(n2 + n3) + (n2 - 2)*n3 + skip
      if(i.gt.(512*nspio)) then
c      include 'inc/rjzztg3.f'
c
        rewind s6
c
        il = bufer(1, 3)
c      include 'inc/rjzztg4.f'
c        read( s6, id = il ) area( il )...area( bufer( 2, 3 ) )
c
        read( s6 )( area( i ), i = il, bufer( 2, 3 ) )
c
        icr6 = 3
        inx6 = 4
      end if
      fil3 = .true.
      limgrn = bufer(2, 2)
      iptgrn = limgrn + 1
      irec3 = 0
      irec6 = 0
 2    n23 = n2*n3
      n13 = n1*n3
      n12 = n1*n2
c
c     preserve boundary densities for smoothing calculation
c
      if(gren) go to 6
c
c
        ip = bufer(1, 5) - 1
c
c       1-boundaries
c
        do 23 i = 1, 2
        is = bd(i) - 1
        do 1 j = 1, n23
 1      area(ip + j) = w(is + j)
23      ip = ip + n23
c
c       2-boundaries
c
        do 4 i = 3, 4
        is = bd(i) + n23 - 1
        do 4 j = 2, n1 - 1
        do 3 k = 1, n3
 3      area(ip + k) = w(is + k)
        is = is + n23
 4      ip = ip + n3
c
c       3-boundaries
c
        n123 = n1*n23
        nn = n12 - 2*n2
        do 7 i = 5, 6
        is = bd(i) + n23
        do 8 j = 1, nn
 8      area(ip + j) = w(is + j*n3 - n3)
 7      ip = ip + nn
        go to 6
c
c     scan stages of solution process
c
 6    do 5 stage = 1, 2
c
c     scan over transform directions
c
      do 5 in = 1, 3
      dir = in
      if(stage.eq.2) dir = 4 - in
      par = dir.eq.3
c
c     sine analysis or synthesis
c
      if(stage.eq.1) then
c
c         fourier analysis and (if appropriate) mesh rotation
c
          call vftsna(n(dir), nom(dir), w(scmbuf))
          if(dir.ne.3) then
            call transp(w(scmbuf), n(dir), nom(dir))
          else
c
c           potential in transform space, with screening charges
c
            call slvptf(bdy)
c
c           boundary potential
c
            if(gren) then
              call vbdgrn(bdy)
            else
              call vbdpot(bdy)
            end if
          end if
      else
c
c         final potential in transform space
c
          if(dir.eq.3) call slvptt(bdy)
c
c         fourier synthesis and (if appropriate) mesh rotation)
c
          if(dir.ne.3) call transp(w(scmbuf), nom(dir), n(dir))
          call vftsns(n(dir), nom(dir), w(scmbuf))
      end if
 5    continue
c
c     smoothing for boundary values if appropriate
c
      if(gren) return
        ip = bufer(1, 5) - 1
c
c       1-boundaries - this section performs full correction for
c       2- and 3- edges
c
        do 11 i = 1, 2
        jp = bd(i) - 1
        do 9 j = 1, n23
 9      w(jp + j) = w(jp + j) - smooth*area(ip + j)
11      ip = ip + n23
c
c       2-boundaries - this section does not affect edge values
c
        do 12 i = 3, 4
        jp = bd(i) + n23
        do 12 j = 2, n1 - 1
        do 10 k = 1, n3 - 2
10      w(jp + k) = w(jp + k) - smooth*area(ip + k + 1)
        jp = jp + n23
12      ip = ip + n3
c
c       3-boundaries - the 1-edges are corrected here
c
        do 13 i = 5, 6
        jp = bd(i) + n23 - n3
        do 14 j = 1, nn
14      area(ip + j) = w(jp + j*n3) - smooth*area(ip + j)
        do 15 j = 1, nn
15      w(jp + j*n3) = area(ip + j)
13      ip = ip + nn
c
c       restore stack pointer.
c
        istack = istck0
c
c     finish
c
      return
      end
