      subroutine vbdtab(bdy, horiz, vert, itrnsp)
c routine written by Richard James for 3-D Poisson solver
      integer  bdy( 6 ), horiz( 6 ), itrnsp( 3 ), vert( 257 )
c
c     purpose
c
c        to calculate the sc, cs and cc transforms of the sums
c        and differences of potentials on the 3 sets of parallel
c        boundaries.  the data for the vertical boundaries is
c        packed into a set of vectors, one for each horizontal
c        plane of the mesh.
c
c     parameters:
c
c     bdy     -  an integer array with 6 elements, pointing
c                on entry to the vectors holding the boundary
c                potentials.  these vectors are replaced by the sums
c                and differences of potentials for parallel pairs.
c
c     horiz   -  an integer array with 6 elements, pointing to the
c                vectors to receive the horizontal boundary transforms.
c
c     vert    -  an integer array with one element per horizontal
c                plane, pointing to the vectors to receive the
c                vertical boundary transforms.  these vectors are
c                assumed to be of length
c                       6*(sum of horizontal dimensions).
c                transforms for the 2-boundaries precede those for
c                the 3-boundaries.
c
c     itrnsp  -  an integer descriptor array of 3 elements, pointing
c                to the vectors holding the transposition indices
c                for boundary planes.
c
c     notes:
c
c     1)   the mesh is assumed to have experienced two rotations
c          during the fourier transformation stage preceding this
c          routine, and so has  n3  planes, each containing  n1  rows
c          of  n2  elements.
c
c     2)   the order adopted for sets of transforms (segmented or
c          otherwise) is cc+, cc-, cs+, cs-, sc+, sc-.
c
c     3)   for efficiency reasons transposition indices are
c          pre-calculated by  vbdpot(xf).
c
      include 'rjinclude.h'
c
c local variables
      integer i, idir, iminus, iplus, is, iseg, isign, istck0, isv, j,
     1        j0, k, llen, len1, len2, lenrow, lr, lsegs, l1, l12,
     2        l13, l2, l23, l3, nrem, nrows, nsegs, num
      real u
c
c      integer bdy(6), horiz(6), vert(257), itrnsp(3), 
      integer cc(2), cs(2),
     1  sc(2), ccv(2, 2:3, 20), csv(2, 2:3, 20), scv(2, 2:3, 20),
     2  indxa(2:3), ind
      logical rem, plus
      integer l(3), ll(3)
c
c     set up array dimensions and check working area size
c
      l1 = n3
      l2 = n1
      l3 = n2
      l23 = l2*l3
      l13 = l1*l3
      l12 = l1*l2
      l(1) = l1
      l(2) = l2
      l(3) = l3
      ll(1) = l23
      ll(2) = l13
      ll(3) = l12
      i = 2*max0(l23, l13, l12)
      if(i.gt.lnwork) then
        write( s2, 200 )i, lnwork
 200    format( 'Insufficient space in /work/ for vbdtab(xg)' /
     1          'space required =', i10, 5x, 'space available =', i10 )
        call crash( 'vbdtab', 'xg 1')
      end if
c
c     calculate number of segments required for vertical boundaries
c
      lenrow = 6*(l2 + l3)
      nrows = min0(lnwork/lenrow, l1)
      nsegs = l1/nrows
      nrem = l1 - nrows*nsegs
      rem = nrem.ne.0
      if((nsegs.gt.20).or.(rem.and.(nsegs.eq.20))) then
        write( s2, 201 )nsegs, lnwork, lenrow, nrem
 201    format( 'Too many vertical segments in vbdtab(xg)' /
     1           ' nsegs =', i10, 5x,
     2           'lnwork =', i10, 5x, 'lenrow =', i10, 5x,
     3           'nrem =', i10 )
        call crash( 'vbdtab', 'xg 2')
      end if
c
c     set pointers for results.
c
      cc(1) = horiz(1)
      cc(2) = horiz(2)
      cs(1) = horiz(3)
      cs(2) = horiz(4)
      sc(1) = horiz(5)
      sc(2) = horiz(6)
c
c     preserve stack pointer.
c
      istck0 = istack
c
c     cycle over required pointers for vertical boundaries.
c
      len1 = nrows*l3
      len2 = nrem*l3
      do 3 idir = 2, 3
      do 2 isign = 1, 2
      do 1 iseg = 1, nsegs
      ccv(isign, idir, iseg) = istack
      csv(isign, idir, iseg) = istack + len1
      scv(isign, idir, iseg) = istack + 2*len1
 1    istack = istack + 3*len1
      if(rem) then
        ccv(isign, idir, nsegs + 1) = istack
        csv(isign, idir, nsegs + 1) = istack + len2
        scv(isign, idir, nsegs + 1) = istack + 2*len2
        istack = istack + 3*len2
      end if
 2    continue
      len1 = nrows*l2
 3    len2 = nrem*l2
c
c     set up pointers to scatter index vectors.
c
      indxa(2) = istack
      istack   = istack + nrows*l3
      if(l2.eq.l3) then
        indxa(3) = indxa(2)
      else
        indxa(3) = istack
        istack   = istack + nrows*l2
      end if
      if(istack.gt.maxstk) then
        maxstk = istack
        mxstid = 'vbdtab'
      end if
c
c     check for stack overflow.
c
      if(istack.gt.lstack) then
        write( s2, 202 )istck0, istack, lstack
 202    format( 'Stack overflow in vbdtab(xg)'/
     1          ' old pointer =', i10, ' new pointer =', i10,
     2          ' limit =', i10 )
        call crash( 'vbdtab', 'xg 3')
        stop
      end if
c
c     cycle over pairs of parallel boundaries
c
      len1 = l2
      len2 = l3
      do 8 idir = 1, 3
      if(idir.eq.3) len2 = l2
c
c     sums and differences of potentials
c
      llen = ll(idir)
      iminus = idir + idir
      iplus = iminus - 1
      do 20 i = 1, llen
      work(i) = wstack(bdy(iplus) + i) + wstack(bdy(iminus) + i)
20    vect(i) = wstack(bdy(iplus) + i) - wstack(bdy(iminus) + i)
      do 21 i = 1, llen
21    wstack(bdy(iplus) + i) = work(i)
      do 22 i = 1, llen
22    wstack(bdy(iminus) + i) = vect(i)
c
c     cycle over signs
c
      do 7 isign = 1, 2
      plus = isign.eq.1
c
c     cs transforms
c
      llen = ll(idir)
      if(.not.plus) then
        do 23 i = 1, llen
23      work(i) = wstack(bdy(iminus) + i)
      end if
      call vfthil(len1, len2, work, .true.)
      if(idir.eq.1) then
c
c       store directly in result area for horizontal boundaries
c
        do 24 i = 1, llen
24      wstack(cs(isign) + i) = work(i)
      else
c
c       cycle over segments for vertical boundaries
c
        is = 0
        llen = nrows*l(5 - idir)
        do 4 iseg = 1, nsegs
        do 25 i = 1, llen
25      wstack(csv(isign, idir, iseg) + i) = work(is + i)
 4      is = is + llen
        if(rem) then
          llen = nrem*l(5 - idir)
          do 26 i = 1, llen
26        wstack(csv(isign, idir, nsegs + 1) + i) = work(is + i)
        end if
      end if
c
c     sc transforms
c
      j0 = 0
      if(idir.gt.1) j0 = l(5 - idir)
      j = ll(idir) + j0
      llen = ll(idir)
      if(plus) then
        do 27 i = 1, llen
27      work(j + i) = wstack(bdy(iplus) + i)
      else
        do 28 i = 1, llen
28      work(j + i) = wstack(bdy(iminus) + i)
      end if
c
c     transpose ss transform for transformation along rows
c
      do 29 i = 1, llen
29    vect(i) = work(j + jstack(itrnsp(idir) + i))
      do 30 i = 1, llen
30    work(j0 + i) = vect(i)
c
c     transform and re-transpose
c
      call vfthil(len2, len1, work(j0 + 1), .true.)
      do 31 i = 1, llen
31    vect(i) = work(j0 + i)
      do 32 i = 1, llen
32    work(j + jstack(itrnsp(idir) + i)) = vect(i)
c
c     store transform
c
      if(idir.eq.1) then
c
c       direct to result area if horizontal plane
c
        do 33 i = 1, l23
33      wstack(sc(isign) + i) = work(j + i)
      else
c
c       cycle over segments for vertical boundaries
c
        is = j - l(5 - idir)
        llen = nrows*l(5 - idir)
        do 5 iseg = 1, nsegs
        do 34 i = 1, llen
34      wstack(scv(isign, idir, iseg) + i) = work(is + i)
 5      is = is + llen
        if(rem) then
          llen = nrem*l(5 - idir)
          do 46 i = 1, llen
46        wstack(scv(isign, idir, nsegs + 1) + i) = work(is + i)
        end if
      end if
c
c     cc transform
c
      call vfthil(len1, len2, work(j + 1), .true.)
      if(idir.eq.1) then
c
c       direct to result area if horizontal plane
c
        do 35 i = 1, l23
35      wstack(cc(isign) + i) = work(j + i)
      else
c
c       cycle over segments for vertical boundaries
c
        is = j
        llen = nrows*l(5 - idir)
        do 6 iseg = 1, nsegs
        do 36 i = 1, llen
36      wstack(ccv(isign, idir, iseg) + i) = work(is + i)
 6      is = is + llen
        if(rem) then
          llen = nrem*l(5 - idir)
          do 37 i = 1, llen
37        wstack(ccv(isign, idir, nsegs + 1) + i) = work(is + i)
        end if
      end if
c
c     end of loop over signs
c
 7    continue
c
c     end of loop over pairs of boundaries
c
 8    len1 = l1
c
c     assemble vertical boundary data
c
c     construct scatter index vectors
c
      k = 3
      if(l2.eq.l3) k = 2
      do 9 idir = 2, k
      llen = l(5 - idir)
      num = nrows*llen
      u = 1.0/float(llen)
      do 38 i = 1, num
38    jstack(indxa(idir) + i) = i
     1                       + (lenrow - llen)*int(u*(float(i) - 0.5))
 9    continue
c
c     cycle over vertical segments for interleaving
c
      lsegs = nsegs
      if(rem) lsegs = lsegs + 1
      isv = 0
      do 11 iseg = 1, lsegs
c
c     interleave transforms in /work/
c
      is = 0
      if(iseg.le.nsegs) then
        len1 = nrows*l3
      else
        len1 = nrem*l3
      end if
      do 10 idir = 2, 3
      ind = indxa(idir)
      llen = l(5 - idir)
      do 39 i = 1, len1
39    work(is + jstack(ind + i)) = wstack(ccv(1, idir, iseg) + i)
      is = is + llen
      do 40 i = 1, len1
40    work(is + jstack(ind + i)) = wstack(ccv(2, idir, iseg) + i)
      is = is + llen
      do 41 i = 1, len1
41    work(is + jstack(ind + i)) = wstack(csv(1, idir, iseg) + i)
      is = is + llen
      do 42 i = 1, len1
42    work(is + jstack(ind + i)) = wstack(csv(2, idir, iseg) + i)
      is = is + llen
      do 43 i = 1, len1
43    work(is + jstack(ind + i)) = wstack(scv(1, idir, iseg) + i)
      is = is + llen
      do 44 i = 1, len1
44    work(is + jstack(ind + i)) = wstack(scv(2, idir, iseg) + i)
c
c     reset vector lengths for second pass.
c
      if(iseg.le.nsegs) then
        len1 = nrows*l2
      else
        len1 = nrem*l2
      end if
10    is = is + llen
c
c     write out rows to result area
c
      is = 0
      lr = nrows
      if(iseg.gt.nsegs) lr = nrem
      do 11 j = 1, lr
      isv = isv + 1
      do 45 i = 1, lenrow
45    wstack(vert(isv) + i) = work(is + i)
11    is = is + lenrow
c
c     finish
c
c
c     reset stack pointer.
c
      istack = istck0
      return
      end
