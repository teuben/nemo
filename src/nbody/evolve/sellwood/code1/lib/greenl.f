      subroutine greenl(lim1, analyt)
c routine written by Richard James for 3-D Poisson solver
      integer lim1
      logical analyt
c
c     purpose
c
c        to border the greens function over a small
c        mesh and extend it to a large mesh, using an analytic form
c        for the function at large distances.
c
c     parameters
c
c     lim1
c
c     analyt  -  (type logical)  specifies use of an analytical form
c                for the greens function in the near zone.
c
      include 'rjinclude.h'
c
c external
      real grenf3
c
c local variables
      integer i, ibufl, in, ir, irecl, ireclo, it, j, jd, jt, j1, j2,
     1        j3, k, kt, n
      real fac, x, lw1, lw2, lw3
c
      integer p, q, r, t, lr1, lr2, recl, rows, lq2, lq3, recs(4),
     1        buf(4)
      logical t1, t2, par, npar
c
c     set buffer length at load time
c
      data ibufl/1995/
c
c preserve value of factor
c
      fac = factor
c
c     jump if greens function prescribed analytically
c
      if(analyt) go to 41
c
c     send to file
c
      r = ((n1 - 1)*n2*n3)/2 + scmbuf - 1
      q = n2*((n1 + 1)/2)
      call scfilu(s3, 'tempft03', 100)
      rewind s3
      do 9 i = 1, q
      p = r + 1
      r = r + n3
      write(s3) (w(j), j = p, r)
9     continue
      rewind s3
      it = n1
      jt = n2
      kt = n3
      rewind s3
c
c     find specification of large mesh
c
      call psnset(lev3, lim1, 1)
c
c     rotate mesh dimensions to compensate for rotations between sine
c     analyses.
c
      q = n3
      n3 = n2
      n2 = n1
      n1 = q
      q = nom(3)
      nom(3) = nom(2)
      nom(2) = nom(1)
      nom(1) = q
      x = wt3
      wt3 = wt2
      wt2 = wt1
      wt1 = x
      write( s2, 200 )n1, n2, n3
 200  format( / 'Mesh rotated for final Green function' /
     +       ' new mesh dimensions are', 3i10 )
c
c     set up greens function on sequential file for large mesh
c
      i = skip + 2*n1*(n2 + n3) + 2*n2*n3
      call scfilu(s6, 'tempft06', i/512 + 5)
      rewind s6
      q = (jt + 1)/2
      it = (it + 1)/2
      recl = n3
      jt = q
      q = (kt + 1)/2
      r = q + n3 - 1
      write(s6) n1, n2, n3, wt1, wt2, wt3, recl
      q = q + scmbuf - 1
      r = r + scmbuf - 1
      kt = kt + scmbuf - 1
c
c     cycle over planes
c
      do 14 i = 1, n1
      t1 = i.le.it
c
c     step over unwanted lines
c
      if(.not.t1) go to 35
      do 34 j = 2, jt
      read(s3) (w(k), k = scmbuf, kt)
34    continue
c
c     cycle over lines
c
35    do 17 j = 1, n2
      t2 = t1.and.(j.le.jt)
      if(t2) read(s3) (w(k), k = scmbuf, kt)
      t = q
      if(t2) t = kt + 1
c
c     fill in distant values
c
      if(r.lt.t) go to 18
      do 16 k = t, r
16    w(k) = factor*grenf3(i - 1, j - 1, k - q)
c
c     send to file
c
18    write(s6) (w(k), k = q, r)
17    continue
c
c     skip any excess lines
c
      if((jt.le.n2).or.(.not.t1)) go to 14
      t = n2 + 1
      do 19 j = t, jt
      read(s3) (w(k), k = scmbuf, kt)
19    continue
14    continue
c
c     read into store
c
41    rewind s6
      read(s6) it, jt, kt, lw1, lw2, lw3, p
      r = n1*n2
      q = scmbuf - 1
      do 25 i = 1, r
      p = q + 1
      q = q + n3
      read(s6) (w(k), k = p, q)
25    continue
c
      factor = 1.0/float((n1 - 1)*(n2 - 1)*(n3 - 1))
c
c    calculate smoothing factor.
c
      smooth = w(scmbuf) - w(scmbuf + 1)
c
c     preserve boundaries
c
      call bdpres(.false., 0, 0)
c
c     scan over transform directions
c
      do 23 in = 1, 3
      par = in.eq.3
      npar = .not.par
c
c     cosine analysis
c
      rows = nom(in)
      n = n1
      if(in.eq.2) n = n2
      if(par) n = n3
      ir = n3
      if(in.ne.2) ir = rows
      call fftcha(n, rows, ir, scmbuf, par)
      if(in.ne.1) go to 23
c
c     normalise
c
      k = rows*n + scmbuf - 1
      do 43 j = scmbuf, k
43    w(j) = factor*w(j)
23    continue
c
c     set up control areas for boundary writes
c
      j1 = n2 - 2
      j2 = n3 - 2
      do 7 i = 1, 3
      recs(i) = 2*j1*j2
      buf(i) = min0((ibufl/(2*j2))*2*j2, recs(i))
      j1 = n1 - 2
      if(i.eq.2) j2 = n2 - 2
7     continue
      recs(4) = 0
      buf(4) = 0
c
c     send to file
c
      rewind s3
      k = n2*n3
      p = n1*n2*n3 + scmbuf - 1
      irecl = buf(1)
      write(s3) n1, n2, n3, wt1, wt2, wt3, irecl
      do 36 i = scmbuf, p, k
      q = i + k - 1
      write(s3) (w(j), j = i, q)
36    continue
c
c     transform boundary planes
c
      ireclo = irecl
      do 33 i = 1, 3
      lr1 = n2
      lr2 = n3
      if(i.gt.1) lr1 = n1
      if(i.gt.2) lr2 = n2
      lq2 = lr1*lr2 + 1
      call bdpres(.true., 1, 2*i - 1)
      call bdpres(.true., lq2, 2*i)
c
c     normalise
c
      lq3 = 2*lq2 - 2
      factor = 1.0/float((lr1 - 1)*(lr2 - 1))
      do 45 j = 1, lq3
45    w(j) = factor*w(j)
      call dtrfm1(1, 1, lr1, lr2, .false., .true.)
      call dtrfm3(1, 1, lr1, lr2, .false., .true.)
      call dtrfm1(lq2, lq2, lr1, lr2, .false., .true.)
      call dtrfm3(lq2, lq2, lr1, lr2, .false., .true.)
c
c     clear unwanted values
c
      j1 = 2*lr2 + 1
      j2 = j1 + lq2 - 1
      do 6 q = 3, lr1
      w(j1) = 0.0
      w(j1 + 1) = 0.0
      w(j2) = 0.0
      w(j2 + 1) = 0.0
      j1 = j1 + lr2
6     j2 = j2 + lr2
      j3 = lq2 + lr1*lr2
      w(j3) = 0.0
c
c     write to file.  the write is offset one element as the
c     cosine transform elements being handled here are offset
c     relative to the sine cosine elements with which they are
c     used in the boundary convolution routines.
c
      j1 = 2*lr2 + 2
      j2 = j1 + lr1*lr2
      jd = j2 - j1
      do 1 q = j1, lq2
      x         = 0.5*(w(q) + w(q + jd))
      w(q + jd) = 0.5*(w(q) - w(q + jd))
 1    w(q)      = x
      write(s3) (w(q), q = j1, lq2)
      write(s3) (w(q), q = j2, j3)
33    continue
      gren = .false.
      permit = .true.
c
c restore value of factor
c
      factor = fac
c
      return
      end
