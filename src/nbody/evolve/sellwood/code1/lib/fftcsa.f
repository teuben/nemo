      subroutine fftcsa(n, no, group, is, unity)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to perform a fourier cosine analysis along one axis of a
c        2 or 3 dimensional mesh
c
      integer group, st, st1, step
     1, iw, i, j, p, q, pi, qi
c
      include 'rjinclude.h'
c
c local variables
      integer lii, incr, is, it, jump, lim, n, na, nb, no
      real c, s, u, v
c
      logical unity, rev, halve
      halve = .false.
      go to 16
      entry fftcha(n, no, group, is, unity)
      halve = .true.
c
c     set initial values
c
16    if(unity) then
        incr = n
        step = 1
      else
        incr = 1
        step = group
      end if
      st = is
c
c     cycle over transforms
c
      do 15 lii = 1, no, group
      it = min0(group, no - lii + 1)
      do 1 i = 1, it
      na = n - 1
      lim = st + na*step
      if(halve) go to 17
      w(st) = 2.0*w(st)
      w(lim) = 2.0*w(lim)
17    lim = step + lim
c
c     fold real section
c
 2    p = st
      jump = na*step
      q = p + jump
      nb = na/2
      do 3 j = 1, nb
      u = w(p)
      w(p) = u + w(q)
      w(q) = u - w(q)
      p = p + step
3     q = q - step
      w(p) = 2.0*w(p)
c
c     test scan complete
c
      st1 = st + step + jump
      if(st1.eq.lim) go to 4
c
c     re-order real group to complex form
c
      p = st1
      q = p + jump - step
      do 5 j = 1, nb
      u = w(p)
      w(p) = w(q)
      w(q) = - u
      p = p + step
5     q = q - step
c
c     adjust last real part
c
      w(p) = - w(p)
c
c     test scan complete
c
      st1 = st1 + jump
      if(st1.eq.lim) go to 4
c
c     fold first complex group
c
c     initialise and set first terms
c
      p = st1
      q = p + jump
      qi = q
      pi = q + jump
      u = 1.414213562373*w(q)
      w(q) = w(p) - u
      w(p) = w(p) + u
      if(nb.eq.1) go to 10
c
c     fold intermediate terms
c
      do 7 j = 2, nb
      p = p + step
      q = q - step
      pi = pi - step
      qi = qi + step
      u = 0.707106781187*(w(q) - w(qi))
      v = 0.707106781187*(w(q) + w(qi))
      w(q) = w(pi) - v
      w(qi) = w(p) - u
      w(p) = w(p) + u
 7    w(pi) = - w(pi) - v
c
c     calculate end terms
c
10    p = p + step
      pi = pi - step
      u = 0.707106781187*(w(p) - w(pi))
      w(pi) = w(p) - u
      w(p) = w(p) + u
c
c     test scan complete
c
      st1 = st1 + 2*jump
      if(st1.eq.lim) go to 4
c
c     fold remaining groups
c
      iw = 1
      rev = .true.
8     rev = .not.rev
      if(rev) go to 13
      iw = iw + 2
      c = fact(iw)
      s = fact(iw + 1)
      go to 14
13    u = c
      c = s
      s = u
c
c     pointers and first term
c
14    p = st1
      q = p + jump
      pi = q + jump
      qi = q
      u = w(q)/c
      w(q) = w(p) - u
      w(p) = w(p) + u
      if(nb.eq.1) go to 11
c
c     fold intermediate terms
c
      do 9 j = 2, nb
      p = p + step
      q = q - step
      pi = pi - step
      qi = qi + step
      u = c*w(q) - s*w(qi)
      v = s*w(q) + c*w(qi)
      w(q) = w(pi) - v
      w(pi) = - w(pi) - v
      w(qi) = w(p) - u
 9    w(p) = w(p) + u
c
c     calculate end terms
c
11    p = p + step
      pi = pi - step
      u = c*w(p) - s*w(pi)
      w(pi) = w(p) - u
      w(p) = w(p) + u
c
c     next group if not end of scan
c
      st1 = st1 + 2*jump
      if(st1.lt.lim) go to 8
c
c     jump back for next fold
c
 4    na = nb
      if(na.ne.1) go to 2
c
c     final reduction
c
6     st1 = st + 3*step
      iw = 1
      p = st + step
      u = w(p)
      w(p) = w(st) - u
      w(st) = w(st) + u
      na = (n - 1)/2
      nb = 2*step
      do 12 j = 2, na
      iw = iw + 1
      p = st1 + step
      u = w(p)/fact(iw)
      w(p) = w(st1) - u
      w(st1) = w(st1) + u
12    st1 = st1 + nb
      if(.not.halve) go to 1
      w(st) = 0.5*w(st)
      st1 = st + step
      w(st1) = 0.5*w(st1)
1     st = st + incr
15    st = st + group*n - it*incr
      return
      end
