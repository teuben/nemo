      subroutine fftsna(n, no, group, is, unity)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to perform a fourier sine analysis along one axis of a
c        2 or 3 dimensional mesh
c
      include 'rjinclude.h'
c
c local variables
      integer i, lii, incr, ip, ipi, iq, iqi, is, it, iw, j, jump, lim,
     1        n, na, nb, no
      real c, c1, s, u
c
      integer group, st, st1, step, p, q, pi, qi
      logical unity, rev
c
c     set initial values
c
      if(unity) then
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
      na = n - 1
      lim = st + na*step
c
c     fold real section
c
 2    p = st + step
      jump = na*step
      q = st + jump - step
      nb = na/2
      if(nb.eq.1) go to 6
      do 3 j = 2, nb
      ip = p
      iq = q
      do 1 i = 1, it
      work1(i) = w(ip)
      w(ip) = work1(i) - w(iq)
      w(iq) = work1(i) + w(iq)
      ip = ip + incr
1     iq = iq + incr
      p = p + step
3     q = q - step
      iq = q
      do 5 i = 1, it
      w(iq) = 2.0*w(iq)
5     iq = iq + incr
c
c     no action required to convert real group to complex form
c
c     test scan complete
c
 6    st1 = st + 2*jump
      if(st1.ge.lim) go to 4
c
c     fold first complex group
c
c     initialise and set first terms
c
      pi = st1
      q = pi + jump
      qi = q
      p = q + jump
      ipi = pi
      iqi = qi
      do 16 i = 1, it
      work1(i) = 1.414213562373*w(iqi)
      w(iqi) = work1(i) - w(ipi)
      w(ipi) = work1(i) + w(ipi)
      ipi = ipi + incr
16    iqi = iqi + incr
      if(nb.eq.1) go to 10
c
c     fold intermediate terms
c
      do 7 j = 2, nb
      pi = pi + step
      qi = qi - step
      p = p - step
      q = q + step
      ip = p
      iq = q
      ipi = pi
      iqi = qi
      do 7 i = 1, it
      work1(i) = 0.707106781187*(w(iq) - w(iqi))
      work2(i) = 0.707106781187*(w(iq) + w(iqi))
      w(iq) = work2(i) - w(ipi)
      w(ipi) = work2(i) + w(ipi)
      w(iqi) = w(ip) - work1(i)
      w(ip) = w(ip) + work1(i)
      ip = ip + incr
      iq = iq + incr
      ipi = ipi + incr
7     iqi = iqi + incr
c
c     calculate end terms
c
10    p = p - step
      pi = pi + step
      ip = p
      ipi = pi
      do 17 i = 1, it
      work1(i) = 0.707106781187*(w(ipi) + w(ip))
      w(ip) = work1(i) - w(ipi)
      w(ipi) = work1(i) + w(ipi)
      ip = ip + incr
17    ipi = ipi + incr
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
14    pi = st1
      qi = pi + jump
      q = qi
      p = q + jump
      ipi = pi
      iqi = qi
      c1 = 1.0/c
      do 18 i = 1, it
      work1(i) = c1*w(iqi)
      w(iqi) = work1(i) - w(ipi)
      w(ipi) = work1(i) + w(ipi)
      ipi = ipi + incr
18    iqi = iqi + incr
      if(nb.eq.1) go to 11
c
c     fold intermediate terms
c
      do 9 j = 2, nb
      pi = pi + step
      qi = qi - step
      q = q + step
      p = p - step
      ip = p
      iq = q
      ipi = pi
      iqi = qi
      do 9 i = 1, it
      work1(i) = c*w(iq) - s*w(iqi)
      work2(i) = s*w(iq) + c*w(iqi)
      w(iq) = work2(i) - w(ipi)
      w(ipi) = work2(i) + w(ipi)
      w(iqi) = w(ip) - work1(i)
      w(ip) = w(ip) + work1(i)
      ip = ip + incr
      iq = iq + incr
      ipi = ipi + incr
9     iqi = iqi + incr
c
c     calculate end terms
c
11    pi = pi + step
      p = p - step
      ip = p
      ipi = pi
      do 19 i = 1, it
      work1(i) = s*w(ip) + c*w(ipi)
      w(ip) = work1(i) - w(ipi)
      w(ipi) = work1(i) + w(ipi)
      ip = ip + incr
19    ipi = ipi + incr
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
      st1 = st + 2*step
      iw = 1
      na = (n - 1)/2
      nb = 2*step
      do 12 j = 2, na
      iw = iw + 1
      p = st1 + step
      ip = p
      iq = st1
      c1 = 1.0/fact(iw)
      do 20 i = 1, it
      work1(i) = c1*w(ip)
      w(ip) = work1(i) - w(iq)
      w(iq) = work1(i) + w(iq)
      ip = ip + incr
20    iq = iq + incr
12    st1 = st1 + nb
      p = st + step
      ip = p
      do 21 i = 1, it
      w(ip) = 2.0*w(ip)
21    ip = ip + incr
15    st = st + group*n
      return
      end
