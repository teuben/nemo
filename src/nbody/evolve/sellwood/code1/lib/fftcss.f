      subroutine fftcss(n, no, group, is, unity)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to perform a fourier cosine analysis along one axis of a
c        2 or 3 dimensional mesh
c
      include 'rjinclude.h'
c
c local variables
      integer i, incr, is, iw, j, jump, lim, n, na, nb, no
      real c, fcc, s, u, v
c
      integer group, st, st1, step, p, q, pi, qi
      logical unity, rev, leave
      leave = .false.
      go to 6
      entry fftchs(n, no, group, is, unity)
      leave = .true.
c
c     set initial values
c
 6    if(unity) then
        incr = n
        step = 1
      else
        incr = 1
        step = group
      end if
      fcc = 0.5
      na = n - 1
17    fcc = 0.5*fcc
      na = na/2
      if(na.ne.1) go to 17
      st = is
c
c     cycle over transforms
c
      do 1 i = 1, no
      lim = st + n*step
c
c     preliminary unfolding
c
      st1 = st + 3*step
      iw = 1
      p = st + step
      if(.not.leave) go to 15
      w(st) = 2.0*w(st)
      w(p) = 2.0*w(p)
15    u = w(st)
      w(st) = u + w(p)
      w(p) = u - w(p)
      p = p + step
      w(p) = 2.0*w(p)
      nb = 2*step
      do 12 j = 5, n, 2
      iw = iw + 1
      p = st1 + step
      u = fact(iw)*(w(st1) - w(p))
      w(st1) = w(st1) + w(p)
      w(p) = u
12    st1 = st1 + nb
      nb = 1
c
c     unfold real section
c
 2    na = 2*nb
      p = st
      jump = na*step
      q = p + jump
      do 3 j = 1, nb
      u = w(p)
      w(p) = u + w(q)
      w(q) = u - w(q)
      p = p + step
3     q = q - step
c
c     test scan complete
c
      st1 = st + step + jump
      if(st1.eq.lim) go to 4
c
c     re-order first complex group to real form
c
      p = st1
      q = p + jump - step
      do 5 j = 1, nb
      u = w(p)
      w(p) = - 2.0*w(q)
      w(q) = 2.0*u
      p = p + step
5     q = q - step
c
c     adjust last real part
c
      w(q) = - w(q)
c
c     test scan complete
c
      st1 = st1 + jump
      if(st1.eq.lim) go to 4
c
c     unfold second complex group to form new first
c
c     initialise and set first terms
c
      p = st1
      q = p + jump
      qi = q
      pi = q + jump
      u = 0.707106781187*(w(p) - w(q))
      w(p) = w(p) + w(q)
      w(q) = u
      if(nb.eq.1) go to 10
c
c     unfold intermediate terms
c
      do 7 j = 2, nb
      p = p + step
      q = q - step
      pi = pi - step
      qi = qi + step
      u = 0.707106781187*(w(p) - w(qi))
      v = 0.707106781187*(w(q) + w(pi))
      w(p) = w(p) + w(qi)
      w(pi) = w(q) - w(pi)
      w(q) = u - v
 7    w(qi) = - u - v
c
c     calculate middle term
c
10    p = p + step
      pi = pi - step
      u = 1.414213562373*(w(p) - w(pi))
      w(p) = w(p) + w(pi)
      w(pi) = w(p) - u
c
c     test scan complete
c
      st1 = st1 + 2*jump
      if(st1.eq.lim) go to 4
c
c     unfold remaining groups
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
      qi = q
      pi = q + jump
      u = c*(w(p) - w(q))
      w(p) = w(p) + w(q)
      w(q) = u
      if(nb.eq.1) go to 11
c
c     unfold intermediate terms
c
      do 9 j = 2, nb
      p = p + step
      q = q - step
      pi = pi - step
      qi = qi + step
      u = w(q) + w(pi)
      v = w(p) - w(qi)
      w(p) = w(p) + w(qi)
      w(pi) = w(q) - w(pi)
      w(q) = c*v - s*u
 9    w(qi) = - c*u - s*v
c
c     calculate middle terms
c
11    p = p + step
      pi = pi - step
      u = w(p) - w(pi)
      w(p) = w(p) + w(pi)
      w(pi) = (c*w(p) - u)/s
c
c     next group if not end of scan
c
      st1 = st1 + 2*jump
      if(st1.ne.lim) go to 8
c
c     jump back for next unfolding stage
c
 4    nb = na
      na = 2*na
      st1 = lim - step
      if(na.lt.n) go to 2
c
c     correct normalisation and end values
c
      p = st
      do 16 j = 1, n
      w(p) = fcc*w(p)
16    p = p + step
      if(leave) go to 1
      p = p - step
      w(p) = 0.5*w(p)
      w(st) = 0.5*w(st)
1     st = st + incr
      return
      end
