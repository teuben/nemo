      subroutine vbdpot(bdy)
c routine written by Richard James for 3-D Poisson solver
      integer bdy( 6 )
c
c     purpose
c
c        to calculate the contribution on the mesh boundary of the
c        boundary charges (including corrective charges).
c
c     parameter
c
c     bdy  -  an integer array with 6 elements, pointing on
c             entry to the vectors holding the charges in the 6
c             boundary planes.  charges are supplied as sine-sine
c             transforms, and the routine overwrites them by the
c             potentials.  these are also held as sine-sine transforms.
c
      include 'rjinclude.h'
c
c local variables
      integer i, ib, idir, iplane, ipm, ipp, is, istck0, isw, is0, j, k,
     1        llen, len1, len12, len2, ll0, ll1, lr, l1, l12, l13, l2,
     2        l2m2, l2p, l23, l3, l3m2, l3p
      real u
c
      integer horiz(6), ccp1, ccm1, csp1, csm1, scp1,
     1 scm1, cc2, cs2, sc2, cc3, cs3, sc3, isprd2, isprd3, rcc, rcs,
     2 rsc, itrnsp(3), vert(257), result(6), dcc, dcs, dsc
      integer l(3), ll(3)
c
c     set up mesh dimensions
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
      l2m2 = l2 - 2
      l3m2 = l3 - 2
c
c     check adequacy of vert
c
      if(l1.gt.257) then
       write( s2, 200 )l
 200   format( 'There are only 257 pointers in array vert',
     1         ' in vbdpot(xf)' / 'rotated mesh dimensions are',
     2         3i10 / ' run stopped' )
        call crash( 'vbdpot', 'xf 1')
      end if
c
c     check mesh is not too small
c
      if((l2.lt.9).or.(l3.lt.9)) then
        write( s2, 201 )l
 201    format( 'Mesh is too small for vbdpot(xf)' /
     1          ' dimensions of rotated mesh are', 3i10 /
     2          ' run stopped by vbdpot' )
        call crash( 'vbdpot', 'xf 2')
      end if
c
c     allocate space for transposition index vectors.  a separate
c     pointer is provided for each direction, and in the case of
c     symmetry 2 or 3 pointers may point to the same vector.
c
c     cycle over directions
c
      len1 = l2
      len2 = l3
      istck0 = istack
      do 1 idir = 1, 3
c
c     reset second dimension if necessary.
c
      if(idir.eq.3) len2 = l2
c
c     decide if vector needs calculating
c
      if((idir.eq.2).and.(len1.eq.l2)) then
        itrnsp(2) = itrnsp(1)
        go to 1
      end if
      if(idir.eq.3) then
        if((len1.eq.l2).and.(len2.eq.l3)) then
          itrnsp(3) = itrnsp(1)
          go to 1
        else if(len2.eq.l3) then
          itrnsp(3) = itrnsp(2)
          go to 1
        end if
      end if
c
c     allocate space on stack for index vector.
c
      itrnsp(idir) = istack
      istack = istack + len1*len2
c
c     reset first dimension and terminate loop over directions.
c
 1    len1 = l1
c
c     allocate space for distribution index vectors.
c
      isprd2 = istack
      isprd3 = isprd2 + l23
      istack = isprd3 + l23
c
c     allocate memory for boundary transforms and set values in the
c     alternative pointers to these transforms.
c
      do 2 i = 1, 5, 2
      horiz(i)      = istack
      result(i)     = horiz(i) + l23
      horiz(i + 1)  = result(i) + l23
      result(i + 1) = horiz(i)
 2    istack        = horiz(i + 1) + l23
      ccp1 = horiz(1)
      ccm1 = horiz(2)
      csp1 = horiz(3)
      csm1 = horiz(4)
      scp1 = horiz(5)
      scm1 = horiz(6)
c
c     allocate memory for vertical transforms.
c
      llen = 6*(l2 + l3)
      do 3 i = 1, l1
 3    vert(i) = istack + llen*(i - 1)
      istack = vert(l1) + l23
      if(istack.gt.maxstk) then
        maxstk = istack
        mxstid = 'vbdpot'
      end if
c
c     check stack space.
c
      if(istack.gt.lstack) then
        write( s2, 202 )istck0, istack, lstack
 202    format( 'Stack overflow in vbdpot(xf)' /
     1          ' old stack pointer =', i10, ' new pointer =',
     2          i10, ' limit =', i10 )
        call crash( 'vbdpot', 'xf 3')
        stop
      end if
c
c     check vectorisation work space.
c
      i = max0(l23, l13, l12)
      if(i.gt.lnvect) then
        write( s2, 204 )l23, l13, l12, lnvect
 204    format( 'Common area /vect/ too small in vbdpot(xf)' /
     1          ' boundary planes have sizes', 3i10 /
     2          ' limit =', i10 )
        call crash( 'vbdpot', 'xf 4')
        stop
      end if
c
c     assign pointers to areas in /work/.
c
      l3p = l3 + l3
      l2p = l2 + l2
      cc2 = 0
      cs2 = cc2 + l3p
      sc2 = cs2 + l3p
      cc3 = sc2 + l3p
      cs3 = cc3 + l2p
      sc3 = cs3 + l2p
c
c     calculate transposition index vectors.
c
      len1 = l2
      len2 = l3
      do 13 idir = 1, 3
      if(idir.eq.3) then
        len2 = l2
        if(len2.eq.l3) go to 13
      end if
      len12 = len1*len2
      if((idir.eq.2).and.(len1.eq.l2)) go to 13
c
c     calculate index vector for current direction.
c
      u = 1.0/float(len1)
      do 20 k = 1, len12
      j = int(u*(float(k) - 0.5))
20    jstack(itrnsp(idir) + k) = k*len2 - len2 + 1 - j*(len12 - 1)
c
c     end of loop over transpostion index vectors.
c
13    len1 = l1
c
c     calculate distribution index vector for 2-boundaries.  values
c     referring to odd 2-indices are incremented by l3 to pick up
c     differences instead of sums.
c
      u = 1.0/float(l3)
      do 21 k = 1, l23
21    jstack(isprd2 + k) = k - l3*int(u*(float(k) - 0.5))
      do 22 k = ((l2 + 1)/2)*l3 + 1, l23
22    jstack(isprd2 + k) = jstack(isprd2 + k) + l3
c
c     calculate distribution index vector for 3-boundaries.  values
c     referring to odd 3-indices are incremented by l2 to pick up
c     differences instead of sums.
c
      do 23 k = 1, (l3 + 1)/2
23    ivect(k) = 0
      do 24 k = (l3 + 3)/2, l3
24    ivect(k) = l2
      is = isprd3 - l3
      do 25 j = 1, l2
      is = is + l3
      do 25 k = 1, l3
25    jstack(is + k) = ivect(k) + j
c
c     tabulate boundary data
c
      call vbdtab(bdy, horiz, vert, itrnsp)
c
c     set initial pointers for horizontal data and results.
c
      dcc = ccp1
      dcs = csp1
      dsc = scp1
      rcc = result(1)
      rcs = result(3)
      rsc = result(5)
c
c     cycle over horizontal planes
c
      isw = 1
      do 11 iplane = 1, l1
c
c     reassign data and result pointers if first odd-index
c
      if((iplane + iplane).eq.(l1 + 3)) then
        dcc = ccm1
        dcs = csm1
        dsc = scm1
        rcc = result(2)
        rcs = result(4)
        rsc = result(6)
        isw = iplane
      end if
c
c     recover data
c
      do 26 k = 1, llen
26    work(k) = wstack(vert(iplane) + k)
      call grget(grein, l23)
      grein(l23 + 1) = 0.0
c
c     ccc terms
c
      do 27 k = 1, l23
27    vect(k) = work(cc2 + jstack(isprd2 + k))
     1        + work(cc3 + jstack(isprd3 + k))
     2        + wstack(dcc + k)
      do 29 k = 1, l23
29    work(llen + k) = vect(k)*grein(k)
c
c     accumulate horizontal terms
c
      if(iplane.eq.isw) then
        do 30 k = 1, l23
30      wstack(rcc + k) = work(llen + k)
      else
        do 31 k = 1, l23
31      wstack(rcc + k) = work(llen + k) + wstack(rcc + k)
      end if
c
c     initial pointers for vertical boundary accumulation.  initially
c     the mesh is divided into 4 areas, 2 to be added for the + terms
c     and 2 for the - terms.
c
      lr = l3
      ll0 = ((l2 - 1)/4)*l3
      ll1 = ll0
      ipp = l23 + llen
      is0 = 0
c
c     cycle over vertical boundary pairs
c
      do 5 idir = 2, 3
c
c     prepare accumulation on current boundaries and perform first
c     summation stage.
c
      ipm = ipp + ll1
      is = llen + lr
      do 32 k = 1, ll1
32    vect(k) = work(is + k) + work(is + ll1 + k)
      is = is + ll1 + ll1
      do 33 k = 1, ll1
33    vect(ll1 + k) = work(is + k) + work(is + ll1 + k)
      do 34 k = 1, lr
34    vect(k) = vect(k) + work(llen + k)
c
c     compress to 2 lines of values for each of the + and the - terms
c
 4    ll1 = ll1/2
      if(ll1.ge.(lr + lr)) then
        do 35 k = 1, ll1
35      work(ipp + k) = vect(k)         + vect(ll1 + k)
        do 36 k = 1, ll1
36      work(ipm + k) = vect(2*ll1 + k) + vect(3*ll1 + k)
        ll1 = ll1/2
        if(ll1.ge.(lr + lr)) then
          do 85 k = 1, ll1
85        vect(k)       = work(ipp + k) + work(ipp + ll1 + k)
          do 86 k = 1, ll1
86        vect(ll1 + k) = work(ipm + k) + work(ipm + ll1 + k)
          go to 4
        else
          do 87 k = 1, 2*ll1
87        vect(k)         = work(ipp + k)
          do 48 k = 1, 2*ll1
48        vect(2*ll1 + k) = work(ipm + k)
        end if
      end if
c
c     final sums overwrite data at beginning of area
c
      do 37 k = 1, lr
37    work(is0 + k)      = vect(k)        + vect(lr + k)
      do 38 k = 1, lr
38    work(is0 + lr + k) = vect(2*lr + k) + vect(3*lr + k)
c
c     adjust pointers, transpose matrix and terminate loop
c
      if(idir.eq.2) then
        lr = l2
        ll1 = ((l3 - 1)/4)*l2
        is0 = 6*l3
        do 39 k = 1, l23
39      vect(k) = work(llen + jstack(itrnsp(1) + k))
        do 40 k = 1, l23
40      work(llen + k) = vect(k)
      end if
 5    continue
c
c     ccs terms - note that the greens function origin is shifted
c     to allow for the offset in the storage of sine transforms
c
      do 41 k = 1, l23
41    vect(k) = (work(cs2 + jstack(isprd2 + k)) + wstack(dcs + k))
      do 50 k = 1, l23
50    work(llen + k) = vect(k)*grein(k + 1)
c
c     accumulate horizontal terms
c
      if(iplane.eq.isw) then
        do 42 k = 1, l23
42      wstack(rcs + k) = work(llen + k)
      else
        do 43 k = 1, l23
43      wstack(rcs + k) = work(llen + k) + wstack(rcs + k)
      end if
c
c     prepare accumulation on 2-boundaries
c
      ipp = llen + l3
      ll1 = ll0
      ipm = ipp + ll1 + ll1
      do 44 k = 1, l3
44    vect(k) = work(ipp + k) + work(llen + k)
      do 45 k = 1, l3
45    work(ipp + k) = vect(k)
c
c     accumulate even and odd terms
c
 6    do 46 k = 1, ll1
46    vect(k)       = work(ipp + k) + work(ipp + ll1 + k)
      do 47 k = 1, ll1
47    vect(ll1 + k) = work(ipm + k) + work(ipm + ll1 + k)
      ll1 = ll1/2
      if(ll1.gt.l3) then
        do 88 k = 1, ll1
88      work(ipp + k) = vect(k)         + vect(ll1 + k)
        do 89 k = 1, ll1
89      work(ipm + k) = vect(2*ll1 + k) + vect(3*ll1 + k)
        ll1 = ll1/2
        if(ll1.gt.l3) then
          go to 6
        else
          do 90 k = 1, 2*l3
90        vect(k)        = work(ipp + k)
          do 91 k = 1, 2*l3
91        vect(2*l3 + k) = work(ipm + k)
        end if
      end if
c
c     final sums overwrite data for cs transforms on 2-boundaries.
c     values on the 3-edges are not required.
c
      is = l3 + l3 + 1
      do 49 k = 1, l3m2
49    work(is + k) = vect(k + 1)        + vect(k + l3 + 1)
      is = is + l3
      do 51 k = 1, l3m2
51    work(is + k) = vect(k + 2*l3 + 1) + vect(k + 3*l3 + 1)
c
c     csc terms
c
      do 52 k = 1, l23
52    vect(k) = work(cs3 + jstack(isprd3 + k)) + wstack(dsc + k)
      do 84 k = 1, l23
84    work(llen + k) = vect(k)
      do 53 k = 1, l2m2*l3
53    work(llen + l3 + k) = work(llen + l3 + k)*grein(2*l3 + k)
c
c     accumulate horizontal terms
c
      if(iplane.eq.isw) then
        do 54 k = 1, l23
54      wstack(rsc + k) = work(llen + k)
      else
        do 55 k = 1, l23
55      wstack(rsc + k) = work(llen + k) + wstack(rsc + k)
      end if
c
c     transpose matrix for 3-boundary summation and prepare
c     accumulation on 3-boundaries
c
      do 56 k = 1, l23
56    vect(k) = work(llen + jstack(itrnsp(1) + k))
      do 57 k = 1, l23
57    work(llen + k) = vect(k)
      ipp = llen + l2
      do 58 k = 1, l2
58    work(ipp + k) = work(ipp + k) + vect(k)
      ll1 = ((l3 - 1)/4)*l2
      ipm = ipp + ll1 + ll1
c
c     accumulate even and odd terms
c
 7    do 59 k = 1, ll1
59    vect(k)       = work(ipp + k) + work(ipp + ll1 + k)
      do 60 k = 1, ll1
60    vect(ll1 + k) = work(ipm + k) + work(ipm + ll1 + k)
      ll1 = ll1/2
      if(ll1.gt.l2) then
        do 92 k = 1, ll1
92      work(ipp + k) = vect(k)         + vect(ll1 + k)
        do 93 k = 1, ll1
93      work(ipm + k) = vect(2*ll1 + k) + vect(3*ll1 + k)
        ll1 = ll1/2
        if(ll1.gt.l2) then
          go to 7
        else
          do 61 k = 1, 2*ll1
61        vect(k)         = work(ipp + k)
          do 63 k = 1, 2*ll1
63        vect(2*ll1 + k) = work(ipm + k)
        end if
      end if
c
c     final sums overwrite data for cs terms on 3-boundaries.
c     values on the 2-edges are not required.
c
      is = 6*l3 + l2 + l2 + 1
      do 62 k = 1, l2m2
62    work(is + k) = vect(k + 1)        + vect(k + l2 + 1)
      is = is + l2
      do 64 k = 1, l2m2
64    work(is + k) = vect(k + 2*l2 + 1) + vect(k + 3*l2 + 1)
c
c     scc terms if appropriate
c
      if(iplane.le.2) go to 10
      do 65 k = 1, l23
65    vect(k) = work(sc2 + jstack(isprd2 + k))
     1        + work(sc3 + jstack(isprd3 + k))
      do 66 k = 1, l23
66    work(llen + k) = vect(k)*grein(k)
c
c     initial pointers for vertical boundary accumulation.
c
      lr = l3
      ll1 = ll0
      ipp = l23 + llen
      is0 = 4*l3
c
c     cycle over vertical boundary pairs
c
      do 9 idir = 2, 3
c
c     prepare accumulation on current boundaries
c
      ipm = ipp + ll1
      is = llen + lr
      do 67 k = 1, ll1
67    vect(k) = work(is + k) + work(is + ll1 + k)
      do 68 k = 1, lr
68    vect(k) = vect(k) + work(llen + k)
      do 69 k = 1, ll1
69    work(ipp + k) = vect(k)
      is = is + ll1 + ll1
      do 70 k = 1, ll1
70    vect(k) = work(is + k) + work(is + ll1 + k)
      do 71 k = 1, ll1
71    work(ipm + k) = vect(k)
c
c     compress to 2 lines of values for + terms and for - terms
c
        ll1 = ll1/2
        if(ll1.ge.(lr + lr)) then
 8      do 73 k = 1, ll1
73      vect(k)       = work(ipp + k) + work(ipp + ll1 + k)
        do 74 k = 1, ll1
74      vect(ll1 + k) = work(ipm + k) + work(ipm + ll1 + k)
        ll1 = ll1/2
        if(ll1.ge.(lr + lr)) then
          do 75 k = 1, ll1
75        work(ipp + k) = vect(k)         + vect(ll1 + k)
          do 76 k = 1, ll1
76        work(ipm + k) = vect(2*ll1 + k) + vect(3*ll1 + k)
          ll1 = ll1/2
          if(ll1.ge.(lr + lr)) then
            go to 8
          else
            do 72 k = 1, 2*ll1
72          vect(k)         = work(ipp + k)
            do 83 k = 1, 2*ll1
83          vect(2*ll1 + k) = work(ipm + k)
          end if
        end if
      end if
c
c     final sums overwrite data for sc transforms
c
      do 77 k = 1, lr
77    work(is0 + k)      = vect(k)        + vect(lr + k)
      do 94 k = 1, lr
94    work(is0 + lr + k) = vect(2*lr + k) + vect(3*lr + k)
c
c     adjust pointers and transpose matrix if looping again
c
      if(idir.eq.2) then
        lr = l2
        ll1 = ((l3 - 1)/4)*l2
        is0 = 6*l3 + 4*l2
        do 78 k = 1, l23
78      vect(k) = work(llen + jstack(itrnsp(1) + k))
        do 79 k = 1, l23
79      work(llen + k) = vect(k)
      end if
 9    continue
c
c     replace results in stack working areas and terminate loop over
c     planes
c
10    do 80 k = 1, llen
80    wstack(vert(iplane) + k) = work(k)
11    continue
c
c     scan over pairs of parallel boundaries for ss contributions.
c     the 'vertical' edges of the array pointed to by bdy( ) are
c     multiplied by 0's inserted into the relevant greens function
c     when it is first calculated.  the horizontal edges (not
c     calculated here) are not required and the values left in them
c     are removed by vbdass.
c
      len1 = l2
      len2 = l3
      do 12 idir = 1, 3
      if(idir.eq.3) len2 = l2
      ib = idir + idir - 1
c
c     get averages of greens functions and form averages of potentials
c     note that the sine transform offset is compensated in the
c     greens function in this part of the calculation.  the densities
c     have been summed and differenced by vbdtab(xg).
c
      len12 = (len1 - 2)*len2
      call grget(grein(len2 + 1), len12)
      do 81 k = len2 + 1, ll(idir) - len2
81    wstack(bdy(ib) + k) = wstack(bdy(ib) + k)*grein(k)
c
c     get 0.5*(differences of greens functions) and form halved
c     potential differences
c
      call grget(grein(len2 + 1), len12)
      do 82 k = len2 + 1, ll(idir) - len2
82    wstack(bdy(ib + 1) + k) = wstack(bdy(ib + 1) + k)*grein(k)
c
c     end of loop over directions
c
12    len1 = l1
c
c     assemble data.
c
      call vbdass(result, vert, itrnsp, bdy)
c
c     release stack and finish.
c
      istack = istck0
      return
      end
