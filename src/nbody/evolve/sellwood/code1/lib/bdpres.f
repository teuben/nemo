      subroutine bdpres(in, place, record)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to preserve and recover boundary planes for use by  bdpot3
c
c     parameters
c
c        in      -  (type logical) indicates if a boundary plane needs
c                   to be recovered from /lev3/ or file.  the value
c                   .false.  indicates that all boundaries should be
c                   preserved, together with part of the mesh if
c                   necessary
c
c        place   -  (type integer) is the number of the first location
c                   in /scm/ to which the boundary plane should be
c                   returned.  this parameter is not used when  .not.in
c
c        record  -  (type integer) identifies the boundary plane pair
c                   required when recovering data.  for value 7 the
c                   subroutine recovers preserved mesh values.  the
c                   parameter is not used in the case (.not.in).
c
      include 'rjinclude.h'
c
c local variables
      integer i, iil, j, jj, k, kk, kl, kp, l, m
c
      logical in
      integer place, record
c
c     jump to recovery section if required
c
      if(in) go to 1
c
c     preserve all boundary planes and (if necessary) part or all of
c     the mesh
c
c     write out 1 - boundaries
c
      k = n2*n3
      j = scmbuf + (n1 - 1)*n2*n3
      if(fileb) go to 5
      l = pl3 + k - 1
      kk = scmbuf
      do 3 jj = pl3, l
      lev3(jj) = w(kk)
3     kk = kk + 1
      kp = l + 1
      l = l + k
      kk = j
      do 4 jj = kp, l
      lev3(jj) = w(kk)
4     kk = kk + 1
      l = l + 1
      go to 6
5     l = scmbuf + k - 1
      write(s7) (w(m), m = scmbuf, l)
      l = j + k - 1
      write(s7) (w(m), m = j, l)
c
c     write out 2 - boundaries
c
6     j = scmbuf
      do 7 i = 1, 2
      if(fileb) go to 10
      do 12 iil = 1, n1
      k = j + n3 - 1
      do 13 jj = j, k
      lev3(l) = w(jj)
13    l = l + 1
12    j = j + n2*n3
      go to 7
10    m = 1
      do 8 iil = 1, n1
      k = j + n3 - 1
      do 9 jj = j, k
      w(m) = w(jj)
9     m = m + 1
8     j = j + n2*n3
      j = n1*n3
      write(s7) (w(k), k = 1, j)
 7    j = scmbuf + (n2 - 1)*n3
c
c     write out 3 - boundaries
c
      m = scmbuf
      do 11 i = 1, 2
      k = n1*n2*n3 + m - 1
      if(fileb) go to 16
      do 15 kk = m, k, n3
      lev3(l) = w(kk)
15    l = l + 1
      go to 11
16    write(s7) (w(kk), kk = m, k, n3)
      l = l + k
11    m = scmbuf + n3 - 1
      if(place.eq.0) return
c
c     preserve part of mesh to create work space if necessary
c
      if(.not.meshpr) return
      if(filep) go to 21
      kk = pl4 + pp - 1
      jj = scmbuf
      do 17 k = pl4, kk
      lev3(k) = w(jj)
17    jj = jj + 1
21    write(s7) (w(m), m = scmbuf, scmbuf + pp - 1)
      return
c
c     find position and size of record
c
1     k = n2*n3
      l = k*min0(record - 1, 2) + pl3
      if(record.lt.3) go to 22
      k = n1*n3
      l = k*min0(record - 3, 2) + l
      if(record.lt.5) go to 22
      k = n1*n2
      l = k*min0(record - 5, 2) + l
      if(record.lt.7) go to 22
c
c     restore mesh values
c
      if(.not.meshpr) return
      if(filep) go to 20
      iil = scmbuf
      l = pll4 - 1
      do 19 k = pp, l
      w(iil) = lev3(k)
19    iil = iil + 1
      return
20    read(s7) (w(k), k = scmbuf, pp)
      return
c
c     bring in boundary plane
c
22    if(fileb) go to 23
      kk = l + k - 1
      jj = place
      do 18 kl = l, kk
      w(jj) = lev3(kl)
18    jj = jj + 1
      return
23    k = place + k - 1
      read(s7) (w(m), m = place, k)
      return
      end
