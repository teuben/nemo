      subroutine fftset(n, m, c, lw, lc)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c       to set up in /factor/ a table of multipliers for the fft package
c       and (optionally) to set up a re-ordering table for each distinct
c       length of transform required
c
c     parameters
c
c     n   -  the number of lengths of transform required
c
c     m   -  an integer array of length  m,  holding the indices for the
c            transfer lengths on entry.  its values do not change.
c
c     c   -  an integer array to receive the re-ordering tables.  the
c            first  n  elements of  c  hold pointers to the
c            corresponding re-ordering tables on exit.  only one table
c            is generated for each distinct length.  it may be referred
c            to by any number of pointers.
c
c     lw  -  the length of the array in /factor/ to receive the
c            multipliers
c
c     lc  -  the length of the array  c.  if the value of  lc  is  1,
c            the re-ordering tables are not generated.
c
c     note:  only one multiplier table is set up, its length being that
c            for the longest transform.  shorter transforms use a subset
c            of this table, starting at its first entry.
c
      integer n, m(n), i, j, k, l, p, q, r, lw, lc, st, pt, c(lc)
c
      include 'rjinclude.h'
c
c local variable
      real x
c
c     find maximum m and tag repeated values
c
      pt = 1
      if(n.eq.1) go to 1
      do 2 i = 2, n
      if(m(pt).lt.m(i)) pt = i
      do 2 j = i, n
 2    if((m(j).eq.m(i - 1)).and.(m(j).gt.0)) m(j) = 1 - i
 1    st = n + 1
c
c     cycle over required index arrays
c
      do 3 i = 1, n
      if(m(i).gt.0) go to 4
c
c     plant link if already available
c
      j = - m(i)
      m(i) = m(j)
      if(lc.eq.1) go to 3
      c(i) = c(j)
      go to 3
c
c     check space is available
c
 4    r = 2**m(i)
      if(lc.eq.1) go to 3
      if((r + st + 2).lt.lc) go to 5
      r = r + 1
      write(6, 200) i, m(i), r, st, lc
200   format('0overflow in integer area for case', i10, 5x, 'index ='
     1, i10, 5x, 'number of points =', i10
     2/'0start address =', i10, 5x, 'limit =', i10)
      stop
c
c     plant link and initial entries
 5    c(i) = st
      l = m(i)
      c(st) = l
      c(st + 1) = r + 1
      c(st + 2) = 0
      c(st + 3) = r
      q = st + 4
      c(q) = r/2
c
c     set up index array
c
      st = q + 1
      p = 1
      do 6 j = 2, l
      do 7 k = 1, p
      c(st) = c(q)/2
      c(st + 1) = r - c(st)
      st = st + 2
 7    q = q + 1
 6    p = 2*p
      p = c(i)
 3    continue
c
c     check space for multipliers
c
10    l = m(pt)
      r = 2**(l - 1)
      if(r.le.lw) go to 8
      write(6, 201) pt, l, r, lw
201   format('0no room for multipliers'/'0', 4i10)
      stop
c
c     set up multiplier array
c
 8    p = 1
      q = 2
      st = 3
      fact(1) = 2**(m(pt) - 2)
      fact(2) = fact(1)
      x = 2.0*fact(1)
      do 12 j = 3, l
      do 9 k = 1, p
      fact(st) = 0.5*fact(q)
      fact(st + 1) = x - fact(st)
      st = st + 2
 9    q = q + 1
12    p = 2*p
      x = 1.57079632679/x
      st = st - 1
      do 11 j = 1, st
11    fact(j) = cos(x*fact(j))
      return
      end
