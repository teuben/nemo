      subroutine fftsns(n, no, group, is, unity)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to perform a fourier sine synthesis along one axis of a
c        2 or 3 dimensional mesh
c
      include 'rjinclude.h'
c
c local variables
      integer i, iil, incr, ip, iq, ipi, iqi, is, it, iw, j, jump, lim,
     1        n, na, nb, no
      real c, c1, fcc, s, u
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
c     calculate normalising factor
c
      fcc = 0.5
      na = n - 1
17    fcc = 0.5*fcc
      na = na/2
      if(na.ne.1) go to 17
c
c     cycle over transforms
c
      do 16 iil = 1, no, group
      it = min0(group, no - iil + 1)
      lim = st + (n - 1)*step
c
c     preliminary reduction
c
      p = st + step
      ip = p
      do 1 i = 1, it
      w(ip) = 2.0*w(ip)
1     ip = ip + incr
      st1 = st + 2*step
      iw = 1
      nb = 2*step
      do 12 j = 5, n, 2
      iw = iw + 1
      p = st1 + step
      c1 = fact(iw)
      ip = p
      iq = st1
      do 18 i = 1, it
      work1(i) = w(iq)
      w(iq) = work1(i) - w(ip)
      w(ip) = c1*(work1(i) + w(ip))
      ip = ip + incr
18    iq = iq + incr
12    st1 = st1 + nb
      na = 1
c
c     unfold real section
c
 2    p = st + step
      nb = na
      na = 2*nb
      jump = na*step
      q = st + jump - step
      if(nb.eq.1) go to 6
      do 3 j = 2, nb
      ip = p
      iq = q
      do 19 i = 1, it
      work1(i) = w(ip) + w(iq)
      w(iq) = w(iq) - w(ip)
      w(ip) = work1(i)
      ip = ip + incr
19    iq = iq + incr
      p = p + step
3     q = q - step
c
c     test scan complete
c
 6    st1 = st + jump
      if(st1.eq.lim) go to 4
c
c     convert complex group to real form
c
      p = st1 + jump - step
      do 15 j = st1, p, step
      ip = j
      do 15 i = 1, it
      w(ip) = 2.0*w(ip)
15    ip = ip + incr
c
c     test scan complete
c
      st1 = st1 + jump
      if(st1.ge.lim) go to 4
c
c     unfold first complex group
c
c     initialise and set first terms
c
      pi = st1
      q = pi + jump
      qi = q
      p = q + jump
      ipi = pi
      iqi = qi
      do 20 i = 1, it
      work1(i) = 0.707106781187*(w(ipi) + w(iqi))
      w(ipi) = w(ipi) - w(iqi)
      w(iqi) = work1(i)
      ipi = ipi + incr
20    iqi = iqi + incr
      if(nb.eq.1) go to 10
c
c     unfold intermediate terms
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
      work1(i) = 0.707106781187*(w(iq) + w(ipi))
      work2(i) = 0.707106781187*(w(ip) - w(iqi))
      w(ip) = w(ip) + w(iqi)
      w(ipi) = w(ipi) - w(iq)
      w(iq) = work1(i) + work2(i)
      w(iqi) = work1(i) - work2(i)
      ipi = ipi + incr
      iqi = iqi + incr
      ip = ip + incr
7     iq = iq + incr
c
c     calculate middle terms
c
10    p = p - step
      pi = pi + step
      ip = p
      ipi = pi
      do 21 i = 1, it
      work1(i) = 1.414213562373*(w(ip) + w(ipi))
      w(ipi) = w(ipi) - w(ip)
      w(ip) = work1(i) - w(ipi)
      ip = ip + incr
21    ipi = ipi + incr
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
14    pi = st1
      qi = pi + jump
      q = qi
      p = q + jump
      ipi = pi
      iqi = qi
      do 22 i = 1, it
      work1(i) = c*(w(ipi) + w(iqi))
      w(ipi) = w(ipi) - w(iqi)
      w(iqi) = work1(i)
      ipi =ipi + incr
22    iqi = iqi + incr
      if(nb.eq.1) go to 11
c
c     unfold intermediate terms
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
      work1(i) = w(ip) - w(iqi)
      work2(i) = w(iq) + w(ipi)
      w(ip) = w(ip) + w(iqi)
      w(ipi) = w(ipi) - w(iq)
      w(iq) = c*work1(i) + s*work2(i)
      w(iqi) = c*work2(i) - s*work1(i)
      ip = ip + incr
      iq = iq + incr
      ipi = ipi + incr
9     iqi = iqi + incr
c
c     calculate middle terms
c
11    pi = pi + step
      p = p - step
      c1 = 1.0/s
      ip = p
      ipi = pi
      do 24 i = 1, it
      work1(i) = w(ip) + w(ipi)
      w(ipi) = w(ipi) - w(ip)
      w(ip) = c1*(work1(i) - c*w(ipi))
      ip = ip + incr
24    ipi = ipi + incr
c
c     next group if not end of scan
c
      st1 = st1 + 2*jump
      if(st1.lt.lim) go to 8
c
c     jump back for mext fold
c
 4    continue
      if((2*na).lt.n) go to 2
c
c     normalise
c
      if(unity) go to 23
      st1 = st + step
      q = lim -  1
      do 5 j = st1, q
 5    w(j) = fcc*w(j)
      go to 16
23    q = st + n*group - 1
      do 25 j = st, q
25    w(j) = fcc*w(j)
      fcc = 1.0/fcc
      ip = st + group*incr - incr
      iq = st + n - 1
      do 26 j = st, ip, incr
      w(j) = fcc*w(j)
      w(iq) = fcc*w(iq)
26    iq = iq + incr
16    st = st + group*n
      return
      end
