      subroutine vftsna(n, no, x)
c routine written by Richard James for 3-D Poisson solver
      integer n, no
      real x( n * no )
c
c     purpose
c
c        to perform fourier sine analysis along the vertical axis of
c        a 3-dimensional mesh.
c
c     parameters
c
c     n   =  transform length
c
c     no  =  number of transforms to be evaluated in parallel
c
c     x   -  a real array holding the data for transformation.
c            the store direction is orthogonal to the transform
c            direction.
c
      include 'rjinclude.h'
c
c local variables
      integer i, inc, ip, iplane, iq, is1, iw, j, jump, lim, na, nb
      real cc, c1, sc, ss, u
c
c      real x(n*no)
      integer p, q, pi, qi
      logical rev
c
c     check stack size.
c
      inc = no
      i = 2*inc
      j = istack + i
      if(j.gt.maxstk) then
        maxstk = j
        mxstid = 'vftsna'
      end if
      if(j.gt.lstack) then
        write( s2, 200 )istack, i, lstack
 200    format( 'Stack overflow in vftsna(xb)'/
     1          ' stack pointer =', i10, ' space required =', i10,
     2          ' limit =', i10 )
        call crash( 'vftsna', 'xa 1')
        stop
      end if
c
c     set work space pointers and control integers.
c
      ip = istack
      iq = ip + inc
      na = n - 1
      lim = na*inc
c
c     fold real section
c
 1    p = inc
      jump = na*inc
      q = jump - inc
      nb = na/2
c
c     cycle over folds in real section
c
      if(nb.eq.1) go to 3
      do 2 iplane = 2, nb
      do 11 i = 1, inc
11    wstack(ip + i) = x(p + i) - x(q + i)
      do 12 i = 1, inc
12    wstack(iq + i) = x(p + i) + x(q + i)
      do 13 i = 1, inc
13    x(p + i) = wstack(ip + i)
      do 14 i = 1, inc
14    x(q + i) = wstack(iq + i)
      p = p + inc
 2    q = q - inc
c
c     double central values
c
      do 15 i = 1, inc
15    x(q + i) = 2.0*x(q + i)
c
c     no action is required to convert real group to complex form
c
c     test scan complete
c
 3    is1 = 2*jump
      if(is1.ge.lim) go to 9
c
c     fold first complex group
c
c     initialise and set first terms
c
      pi = is1
      q = pi + jump
      qi = q
      p = q + jump
      do 16 i = 1, inc
16    wstack(ip + i) = 1.414213562373*x(qi + i) - x(pi + i)
      do 44 i = 1, inc
44    wstack(iq + i) = 1.414213562373*x(qi + i) + x(pi + i)
      do 17 i = 1, inc
17    x(qi + i) = wstack(ip + i)
      do 18 i = 1, inc
18    x(pi + i) = wstack(iq + i)
      if(nb.eq.1) go to 5
c
c     fold intermediate terms
c
      do 4 j = 2, nb
      pi = pi + inc
      qi = qi - inc
      p = p - inc
      q = q + inc
      do 19 i = 1, inc
19    wstack(ip + i) = 0.707106781187*(x(q + i) - x(qi + i))
      do 20 i = 1, inc
20    wstack(iq + i) = 0.707106781187*(x(q + i) + x(qi + i))
      do 45 i = 1, inc
45    vect(i) = wstack(iq + i) - x(pi + i)
      do 21 i = 1, inc
21    x(q + i)  = vect(i)
      do 22 i = 1, inc
22    x(pi + i) = wstack(iq + i) + x(pi + i)
      do 46 i = 1, inc
46    vect(i) = x(p + i) - wstack(ip + i)
      do 23 i = 1, inc
23    x(qi + i) = vect(i)
      do 24 i = 1, inc
24    x(p + i)  = x(p + i) + wstack(ip + i)
 4    continue
c
c     calculate end terms
c
 5    p = p - inc
      pi = pi + inc
      do 25 i = 1, inc
25    wstack(ip + i) = 0.707106781187*(x(pi + i) + x(p + i))
      do 47 i = 1, inc
47    vect(i) = wstack(ip + i) - x(pi + i)
      do 26 i = 1, inc
26    x(p + i)  = vect(i)
      do 27 i = 1, inc
27    x(pi + i) = wstack(ip + i) + x(pi + i)
c
c     test scan complete
c
      is1 = is1 + 2*jump
      if(is1.eq.lim) go to 9
c
c     fold remaining groups
c
      iw = 1
      rev = .true.
 6    rev = .not.rev
      if(rev) then
        u = cc
        cc = ss
        ss = u
      else
        iw = iw + 2
        cc = fact(iw)
        ss = fact(iw + 1)
      end if
c
c     pointers and first term
c
      pi = is1
      qi = pi + jump
      q = qi
      p = q + jump
      c1 = 1.0/cc
      sc = c1*ss
      do 28 i = 1, inc
28    wstack(ip + i) = c1*x(qi + i)
      do 48 i = 1, inc
48    vect(i) = wstack(ip + i) - x(pi + i)
      do 29 i = 1, inc
29    x(qi + i) = vect(i)
      do 30 i = 1, inc
30    x(pi + i) = wstack(ip + i) + x(pi + i)
      if(nb.eq.1) go to 8
c
c     fold intermediate terms
c
      do 7 j = 2, nb
      pi = pi + inc
      qi = qi - inc
      q = q + inc
      p = p - inc
      do 31 i = 1, inc
31    wstack(ip + i) = x(q + i) - sc*x(qi + i)
      do 32 i = 1, inc
32    wstack(iq + i) = x(qi + i) + sc*x(q + i)
      do 33 i = 1, inc
33    vect(i)  = cc*wstack(iq + i) - x(pi + i)
      do 49 i = 1, inc
49    x(q + i) = vect(i)
      do 34 i = 1, inc
34    x(pi + i) = cc*wstack(iq + i) + x(pi + i)
      do 35 i = 1, inc
35    vect(i) = x(p + i) - cc*wstack(ip + i)
      do 50 i = 1, inc
50    x(qi + i) = vect(i)
      do 36 i = 1, inc
36    x(p + i)  = x(p + i) + cc*wstack(ip + i)
 7    continue
c
c     calculate end terms
c
 8    pi = pi + inc
      p = p - inc
      do 37 i = 1, inc
37    wstack(ip + i) = x(pi + i) + sc*x(p + i)
      do 38 i = 1, inc
38    vect(i) = cc*wstack(ip + i) - x(pi + i)
      do 51 i = 1, inc
51    x(p + i) = vect(i)
      do 39 i = 1, inc
39    x(pi + i) = cc*wstack(ip + i) + x(pi + i)
c
c     next group if not end of scan
c
      is1 = is1 + 2*jump
      if(is1.lt.lim) go to 6
c
c     jump back for next fold
c
 9    na = nb
      if(na.ne.1) go to 1
c
c     final reduction
c
      is1 = 2*inc
      iw = 1
      na = (n - 1)/2
      nb = 2*inc
      do 10 j = 2, na
      iw = iw + 1
      p = is1 + inc
      q = is1
      c1 = 1.0/fact(iw)
      do 40 i = 1, inc
40    wstack(ip + i) = c1*x(p + i)
      do 41 i = 1, inc
41    vect(i) = wstack(ip + i) - x(q + i)
      do 52 i = 1, inc
52    x(p + i) = vect(i)
      do 42 i = 1, inc
42    x(q + i) = wstack(ip + i) + x(q + i)
10    is1 = is1 + nb
      p = inc
      do 43 i = 1, inc
43    x(p + i) = 2.0*x(p + i)
c
c     finish
c
      return
      end
