      subroutine vftsns(n, no, x)
c routine written by Richard James for 3-D Poisson solver
      integer n, no
      real x( n * no )
c
c     purpose
c
c        to perform fourier sine synthesis along the vertical axis of
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
      integer i, inc, ip, iq, is1, iw, j, jump, lim, na, nb
      real cc, cs, c1, fcc, sc, ss, u
c
c      real x(n*no)
      integer p, q, pi, qi
      logical rev
c
c     set up pointers and check stack.
c
      inc = no
      na = n - 1
      lim = na*inc
      ip = istack
      iq = ip + inc
      i = 2*inc
      j = istack + i
      if(j.gt.maxstk) then
        maxstk = j
        mxstid = 'vftsns'
      end if
      if((istack + i).gt.lstack) then
        write( s2, 200 )istack, i, lstack
 200    format( 'Stack overflow in vftsna(xc)' /
     1          ' stack pointer =', i10, ' space required =', i10,
     2          ' limit =', i10 )
        call crash( 'vftsns', 'xa 1')
        stop
      end if
c
c     preliminary reduction
c
      p = inc
      is1 = 2*inc
      do 14 i = 1, inc
14    x(p + i) = 2.0*x(p + i)
      iw = 1
      nb = 2*inc
      do 2 j = 5, n, 2
      iw = iw + 1
      p = is1 + inc
      q = is1
      c1 = fact(iw)
      do 15 i = 1, inc
15    wstack(ip + i) = x(q + i) - x(p + i)
      do 13 i = 1, inc
13    wstack(iq + i) = x(q + i) + x(p + i)
      do 16 i = 1, inc
16    x(q + i) = wstack(ip + i)
      do 17 i = 1, inc
17    x(p + i) = c1*wstack(iq + i)
 2    is1 = is1 + nb
      na = 1
c
c     unfold real section
c
 3    p = inc
      nb = na
      na = na + na
      jump = na*inc
      q = jump - inc
      if(nb.eq.1) go to 5
      do 4 j = 2, nb
      do 18 i = 1, inc
18    wstack(ip + i) = x(p + i) + x(q + i)
      do 19 i = 1, inc
19    vect(i)  = x(q + i) - x(p + i)
      do 47 i = 1, inc
47    x(q + i) = vect(i)
      do 20 i = 1, inc
20    x(p + i) = wstack(ip + i)
      p = p + inc
 4    q = q - inc
c
c     test scan complete
c
 5    is1 = jump
      if(is1.eq.lim) go to 12
c
c     convert complex group to real form
c
      p = is1 + jump - inc
      do 6 j = is1, p, inc
      do 6 i = 1, inc
 6    x(j + i) = 2.0*x(j + i)
c
c     test scan complete
c
      is1 = is1 + jump
      if(is1.ge.lim) go to 12
c
c     unfold first complex group
c
c     initialise and set first terms
c
      pi = is1
      q = pi + jump
      qi = q
      p = q + jump
      do 21 i = 1, inc
21    wstack(ip + i) = 0.707106781187*(x(pi + i) + x(qi + i))
      do 22 i = 1, inc
22    vect(i) = x(pi + i) - x(qi + i)
      do 48 i = 1, inc
48    x(pi + i) = vect(i)
      do 23 i = 1, inc
23    x(qi + i) = wstack(ip + i)
      if(nb.eq.1) go to 8
c
c     unfold intermediate terms
c
      do 7 j = 2, nb
      pi = pi + inc
      qi = qi - inc
      p = p - inc
      q = q + inc
      do 24 i = 1, inc
24    wstack(ip + i) = 0.707106781187*(x(q + i) + x(pi + i))
      do 25 i = 1, inc
25    wstack(iq + i) = 0.707106781187*(x(p + i) - x(qi + i))
      do 26 i = 1, inc
26    vect(i) = x(p + i)  + x(qi + i)
      do 49 i = 1, inc
49    x(p + i) = vect(i)
      do 27 i = 1, inc
27    vect(i) = x(pi + i) - x(q + i)
      do 50 i = 1, inc
50    x(pi + i) = vect(i)
      do 28 i = 1, inc
28    x(q + i)  = wstack(ip + i) + wstack(iq + i)
      do 29 i = 1, inc
29    x(qi + i) = wstack(ip + i) - wstack(iq + i)
 7    continue
c
c     calculate middle terms
c
 8    p = p - inc
      pi = pi + inc
      do 30 i = 1, inc
30    wstack(ip + i) = 1.414213562373*(x(p + i) + x(pi + i))
      do 31 i = 1, inc
31    vect(i) = x(pi + i) - x(p + i)
      do 51 i = 1, inc
51    x(pi + i) = vect(i)
      do 32 i = 1, inc
32    x(p + i) = wstack(ip + i) - vect(i)
c
c     test scan complete
c
      is1 = is1 + jump + jump
      if(is1.eq.lim) go to 12
c
c     unfold remaining groups
c
      iw = 1
      rev = .true.
 9    rev = .not.rev
      if(rev) then
        u = cc
        cc = ss
        ss = u
      else
        iw = iw + 2
        cc = fact(iw)
        ss = fact(iw + 1)
      end if
      sc = ss/cc
c
c     pointers and first term
c
      pi = is1
      qi = pi + jump
      q = qi
      p = q + jump
      do 33 i = 1, inc
33    wstack(ip + i) = cc*(x(pi + i) + x(qi + i))
      do 34 i = 1, inc
34    vect(i) = x(pi + i) - x(qi + i)
      do 52 i = 1, inc
52    x(pi + i) = vect(i)
      do 35 i = 1, inc
35    x(qi + i) = wstack(ip + i)
      if(nb.eq.1) go to 11
c
c     unfold intermediate terms
c
      do 10 j = 2, nb
      pi = pi + inc
      qi = qi - inc
      q = q + inc
      p = p - inc
      do 36 i = 1, inc
36    wstack(ip + i) = cc*(x(p + i) - x(qi + i))
      do 37 i = 1, inc
37    wstack(iq + i) = cc*(x(q + i) + x(pi + i))
      do 38 i = 1, inc
38    vect(i) = x(p + i) + x(qi + i)
      do 53 i = 1, inc
53    x(p + i) = vect(i)
      do 39 i = 1, inc
39    vect(i) = x(pi + i) - x(q + i)
      do 54 i = 1, inc
54    x(pi + i) = vect(i)
      do 40 i = 1, inc
40    x(q + i)  = wstack(ip + i) + sc*wstack(iq + i)
      do 41 i = 1, inc
41    x(qi + i) = wstack(iq + i) - sc*wstack(ip + i)
10    continue
c
c     calculate middle terms
c
11    pi = pi + inc
      p = p - inc
      c1 = 1.0/ss
      cs = c1*cc
      do 42 i = 1, inc
42    wstack(ip + i) = c1*(x(p + i) + x(pi + i))
      do 43 i = 1, inc
43    vect(i) = x(pi + i) - x(p + i)
      do 55 i = 1, inc
55    x(pi + i) = vect(i)
      do 44 i = 1, inc
44    x(p + i)  = wstack(ip + i) - cs*vect(i)
c
c     next group if not end of scan
c
      is1 = is1 + jump + jump
      if(is1.lt.lim) go to 9
c
c     jump back for next fold
c
12    if((na + na).lt.n) go to 3
c
c     normalise
c
      fcc = 0.5/float(n - 1)
      do 45 i = inc + 1, lim
45    x(i) = fcc*x(i)
c
c     finish
c
      return
      end
