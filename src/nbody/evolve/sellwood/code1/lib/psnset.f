      subroutine psnset(cc, lim1, lim2)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to set up tables required by the potential solver
c
      include 'rjinclude.h'
      integer m(3)
      equivalence (m(1), m1)
c
c     local variables
      integer i, ib, j, l, lim1, lim2, ll
c
      integer p(3), cc(lim2), is(9)
      logical set(11), test
      character*3 a(12), order
      data a/'end', 'mes', 'cha', '***', '***', 'par', 'gre',
     1 'inh', 'sel', 'poi', 'sou', 'pot'/
      data set/2*.false., 2*.true., .false., 6*.true./
      write( 8, 220)
220   format('1/ctrl/ set up as follows')
c
c     initialise record of setting up
c
      gren = .false.
      permit = .true.
c
c     read and decode instructions
c
 2    read( 7, 100) order, p
100   format(a3, 7x, 3i10)
      do 3 i = 1, 12
      if(order.eq.a(i)) go to 4
3     continue
c
c     nonsensical instruction
c
      write( 8, 200) order
200   format('0request to obey order', 10x, a3)
      write( 8, 205)
205   format('0run terminated by subroutine psnset(tj)')
      stop
4     if(i.gt.1) set(i - 1) = .true.
      gren = (i.eq.7).or.gren
      permit = (i.ne.8).and.permit
      if(i.ne.1) write( 8, 202) order
202   format(' order read', 10x, a3)
      go to (81, 86, 88, 2, 2, 91, 2, 2, 6, 7, 7, 7), i
c
c     return section
c
81    test = .true.
      do 5 i = 1, 11
5     test = test.and.set(i)
      if(.not.test) go to 14
      skip = n1*n2*n3
      scmbuf = 1
      base = 1
      brick2 = 1
      brick3 = 1
      call fftset(3, m, cc, lim1, lim2)
c
c     set up control integers for /factor/
c
      i = max0(n1, n2, n3)
      ctrlc = (i + 1)/2
      ib = ctrlc + i
      if((i + ib - 1).le.lim1) go to 19
      write( 8, 212) lim1, ctrlc, ib
212   format('0space available is', i10/'0control integers are', 2i10)
      write(s2, 205)
      stop
19    bd(7) = 2*n3*(n1 + n2)
      do 18 j = 9, 18
18    bd(j) = bd(j - 1) + n3
c
c     check that space is adequate in /scm/
c
      if((scmbuf + skip - 1).le.lnscm) go to 1
c
c     insufficient space - terminate run
c
      write(s2, 204) scmbuf, skip, lnscm
204   format('0work space starts at', i10, 5x, 'length required is',
     1 i10, 5x, 'but area ends at', i10)
      write(s2, 205)
      stop
c
c     put mesh area at end of /scm/
c
 1    scmbuf = lnscm - skip + 1
c
c     set pointers to boundary planes
c
      bd(1) = scmbuf
      bd(2) = scmbuf + (n1 - 1)*n2*n3
      bd(3) = scmbuf
      bd(4) = scmbuf + (n2 - 1)*n3
      bd(5) = scmbuf
      bd(6) = scmbuf + n3 - 1
c
c     calculate conversion factors
c
      j = ctrlc
      ll = i
      if(ll.le.lnscm) go to 21
      write( 8, 213) lnscm, ll
213   format('0insufficient space for calculating conversion factors'
     1/'0space available =', i10, 5x, 'dimension is', i10)
      write( 8, 205)
      stop
21    do 22 l = 1, ll
22    w(l) = 0.0
      w(2) = 1.0
      call fftcsa(ll, 1, 1,1, .true.)
      do 20 l = 3, ll
      fact(j) = 2.0 - w(l)
20    j = j + 2
c
c     calculate sine values
c
      do 24 l = 1, ll
24    w(l) = 0.0
      w(2) = 0.5
      call fftsna(ll, 1, 1, 1, .true.)
      j = ctrlc + 1
      do 23 l = 2, ll
      fact(j) = w(l)
23    j = j + 2
      call markrs
      write( 8, 214)
214   format('0setting up complete')
      permit = permit.and.(.not.gren)
      return
14    write( 8, 201) (a(i + 1), set(i), i = 1, 11)
      call crash( 'psnset', 'tj 1')
201   format('0not all options set'/'0', 7(/1x, a3, 5x, l10))
      write( 8, 205)
      stop
c
c     mesh size
c
86    write( 8, 206) p
206   format(' mesh indices are', 3x, 3i10)
      do 15 i = 1, 3
      m(i) = p(i)
      if(p(i).lt.0) stop
15    continue
      n1 = 2**p(1) + 1
      n2 = 2**p(2) + 1
      n3 = 2**p(3) + 1
      twodim = .false.
      write( 8, 207) n1, n2, n3
207   format(' mesh size is', 7x, 3i10)
      go to 2
c
c     read channel numbers
c
88    read( 7, *) schan, s7, s8, s9
      do 8 i = 1, 6
 8    is(i) = schan(i)
      is(7) = s7
      is(8) = s8
      is(9) = s9
      write( 8, 209) (i, is(i), i = 1, 9)
209   format(' channel allocation', 9(/1x, 2i10))
      go to 2
c
c     read mesh parameters
c
91    read( 7, *) wt1, wt2, wt3
      write( 8, 216) wt1, wt2, wt3
216   format(' mesh geometry parameters are', 5x, 1p3e20.7)
      go to 2
c
c     replace default set of test points
c
 6    nverif = p(1)
      if((nverif.gt.50).or.(nverif.le.0)) call crash( 'psnset', 'tj 2')
      read( 7, *) ((iverif(i, j), i = 1, 3), j = 1, nverif)
      write( 8, 203) nverif, ((iverif(i, j), i = 1, 3), j = 1, nverif)
203   format('0verify points changed, new set contains', i10//
     1 ' new points'//100(/1x, 3i10))
      go to 2
c
c     set a mode for check on potential solver
c
 7    verify(i - 9) = .true.
      go to 2
      end
