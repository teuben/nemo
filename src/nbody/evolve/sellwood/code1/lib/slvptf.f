      subroutine slvptf(bdy)
c routine written by Richard James for 3-D Poisson solver
      integer bdy( 6 )
c
c     purpose
c
c        to solve the poisson equation in transform space.  when not
c        calculating a greens function the routine calculates
c        boundary screening charges.
c
c     parameter:
c
c     bdy   -  an integer array pointing to vectors to receive the
c              screening charges.
c
c     notes:
c
c     1)  the routine operates on the mesh after the two mesh rotations
c         required in the fourier transform stage.  the mesh therefore
c         contains n3 planes, each with n1 rows of n2 elements.
c
c     2)  the logical scalar  permit  in /marks/ enables the
c         accumulation of screening charges .in this case, if
c         .not.permit, no information is written to the boundary
c         result areas pointed to by bdy.
c
      include 'rjinclude.h'
c
c local variables
      integer i, ip, iplane, ip23, is, iswpln, is1, it, j, k, l, ll,
     1        lll, ll12, ll13, lm2, lpl, lw, l1, l12, l13, l2, l2m2,
     2        l23, l3, l3m2
      real fn1, fn2, fn3, f1, f2, f3, u, v, wa1, wa2, wa3, wb1, wb2,
     1     wb3, x
c
      integer lp1, lp2, lp3, iad(4)
c      integer bdy(6), 
      integer writ23, writ13, writ12, rfac23, rfac13, rfac12,
     1        cfac23, cfac13, cfac12, ar1, ar2, temp, w23a, w23b,
     2        w13a, w13b, w12a, w12b, indx, ass23, e1(4), e2(4),
     3        e3(4), edge1, edge2, edge2w, edge3, ww13a, ww13b,
     4        ww12a, ww12b, fac1, fac2, fac3, sfac1, sfac3,
     5        a12, a13
      logical form1, form2, form3, equal
c
c     return if greens function calculation.
c
      if(gren) return
c
c     set up controls
c
      form1 = formul.eq.1
      form2 = formul.eq.2
      form3 = formul.eq.3
c
c     calculate plane sizes and associated quantities
c
      l1 = n3
      l2 = n1
      l3 = n2
      lm2 = l1 - 2
      l2m2 = l2 - 2
      l3m2 = l3 - 2
      l23 = l2m2*l3
      l13 = lm2*l3
      l12 = lm2*l2
      ll = l2*l3
      lll = l1*ll
      ll13 = l1*l3
      ll12 = l1*l2
c
c     weights and scaling factors
c
      wa1 = wt3
      wa2 = wt1
      wa3 = wt2
      wb1 = w3
      wb2 = w1
      wb3 = w2
      fn1 = 1.0/float(l1 - 1)
      fn2 = 1.0/float(l2 - 1)
      fn3 = 1.0/float(l3 - 1)
c
c     check plane sizes
c
      i = 2*max0(ll, ll13, ll12)
      if(i.gt.lnwork) then
        write( s2, 200 )ll, ll13, ll12, i, lnwork
 200    format( ' area  /work/  too small in vbdpot(xe)' /
     1          ' plane sizes are', 3i10, 5x, 'length required =',
     2          i10 / ' actual length =', i10 )
        call crash( 'slvptf', 'xe 1')
      end if
c
c     check size of vectorisation work area.
c
      if(max0(ll, ll13, ll12).gt.lnvect) then
        write( s2, 201 )ll, ll13, ll12, lnvect
 201    format( ' area /vect/ smaller than a boundary plane' /
     1          ' plane sizes are', 3i10, '   limit =', i10 )
        call crash( 'slvptf', 'xe 2')
        stop
      end if
c
c     set up pointers to new stack.
c
c     index vector for access to 13-planes of mesh.
c
      indx   = istack
c
c     pointers for temporary arrays of twiddle factors and sines.
c
      fac1   = indx  + l13
      sfac1  = fac1  + lm2
      l      = sfac1 + lm2
      if(l1.eq.l2) then
        fac2 = fac1
      else
        fac2 = l
        l     = l + l2m2
      end if
      if(l1.eq.l3) then
        fac3   = fac1
        sfac3  = sfac1
      else if(l2.eq.l3) then
        fac3   = fac2
      else
        fac3   = l
        l      = l + l3m2
      end if
      if(l1.ne.l3) then
        sfac3  = l
        l      = l + l3m2
      end if
c
c     pointers for write masks.  these masks consist of unity where
c     writing is permitted and zero where it is not.
c
      lpl = max0(l23, l13)
      writ23   = l
      writ13   = l
      l        = l + lpl
      if(l2.eq.l3) then
        writ12 = writ13
      else
        writ12 = l
        l      = l + l12
      end if
c
c     pointers for row factor arrays.
c
      rfac23   = l
      rfac13   = rfac23
      l        = l + lpl
      if(l2.eq.l3) then
        rfac12 = rfac13
      else
        rfac12 = l
        l      = l + l12
      end if
c
c     pointers for column factor arrays.
c
      cfac23   = l
      l        = l + lpl
      cfac13 = cfac23
      if(l2.eq.l3) then
        cfac12 = cfac13
      else
        cfac12 = l
        l      = l + l12
      end if
c
c     pointers for factor contributions for planes.
c
      ar1      = l
      ar2      = l + l23
      temp     = ar2 + l23
      l        = temp + ll
c
c     pointers to factor arrays for 13- and 12-planes.
c
      if(form1) then
c
c       case of formula 1.
c
        if((l1.eq.l2).and.(abs(wa1 - wa2).lt.1.0e-6)) then
          a13  = ar1
        else
          a13  = l
          l    = l + l13
        end if
        if((l2.eq.l3).and.(abs(wa2 - wa3).lt.1.0e-6)) then
          a12  = a13
        else
          a12  = l
          l    = l + l12
        end if
      else if(form2) then
c
c       case of formula 2.
c
        equal  = (abs(wa1 - wa2).lt.1.0e-6).and.
     1           (abs(wa1 - wa3).lt.1.0e-6)
        if(equal.and.(l1.eq.l2)) then
          a13  = ar1
        else
          a13  = l
          l    = l + l13
        end if
        if(equal.and.(l2.eq.l3)) then
          a12  = a13
        else
          a12  = l
          l    = l + l12
        end if
      else
c
c       case of formula 3.
c
        if(l1.eq.l2) then
          a13  = ar1
        else
          a13  = l
          l    = l + l13
        end if
        if(l2.eq.l3) then
          a12  = a13
        else
          a12  = l
          l    = l + l12
        end if
      end if
c
c     assign stack space for edges.
c
      edge1  = l
      edge2  = edge1  + lm2
      edge2w = edge2  + l2m2
      edge3  = edge2w + l2m2
      l      = edge3  + l3m2
      do 6 i = 1, 4
      e1(i)  = l
      e2(i)  = e1(i) + lm2
      e3(i)  = e2(i) + l2m2
 6    l      = e3(i) + l3m2
      if(permit) then
c
c       working arrays for screening charge accumulation.  note that
c       the area pointed to by  w13b  must immediately follow that
c       pointed to by  w13a  (and similarly for  w12b, w12a), as these
c       areas are accessed together when gathering for the edge-
c       adjacent terms.
c
        w23a   = l
        w23b   = l     + l23
        w13a   = w23b  + l23
        w13b   = w13a  + l13
        w12a   = w13b  + l13
        w12b   = w12a  + l12
        ww13a  = w12b  + l12
        ww13b  = ww13a + l13
        ww12a  = ww13b + l13
        ww12b  = ww12a + l12
        l      = ww12b + l12
      end if
      if(l.gt.maxstk) then
        maxstk = l
        mxstid = 'slvptf'
      end if
c
c     check stack space.
c
      if(l.ge.lstack) then
        write( s2, 202 )istack, l, lstack
 202    format( ' stack overflow in vslvpt(xe)' /
     1          ' stack pointer =', i10, '   last pointer =', i10,
     2          ' limit =', i10 )
        call crash( 'slvptf', 'xe 3')
        stop
      end if
c
c     construct index vectors for access to 13-planes of mesh
c
      u = 1.0/float(l3)
      k = ll - l3
      do 51 l = 1, l13
51    jstack(indx + l) = k*int(u*(float(l) - 0.5)) + l
c
c     collect twiddle and weight factors
c
      do 52 l = 1, lm2
52    wstack(fac1 + l)  = fact(ctrlc + 2*l - 2)
      do 53 l = 1, lm2
53    wstack(sfac1 + l) = fact(ctrlc + 2*l - 1)
      if(l1.ne.l2) then
        do 54 l = 1, l2m2
54      wstack(fac2 + l) = fact(ctrlc + 2*l - 2)
      end if
      if((l1.ne.l3).and.(l2.ne.l3)) then
        do 55 l = 1, l3m2
55      wstack(fac3 + l) = fact(ctrlc + 2*l - 2)
      end if
      if(l1.ne.l3) then
        do 56 l = 1, l3m2
56      wstack(sfac3 + l) = fact(ctrlc + 2*l - 1)
      end if
c
c     write masks for 23- and 13-planes
c
      do 57 l = 1, lpl
57    wstack(writ23 + l) = 1.0
      do 58 l = 1, lpl, l3
58    wstack(writ23 + l) = 0.0
      do 59 l = l3, lpl, l3
59    wstack(writ23 + l) = 0.0
c
c     write mask for 12-planes
c
      if(l2.ne.l3) then
        do 60 l = 1, l12
60      wstack(writ12 + l) = 1.0
        do 61 l = 1, l12, l2
61      wstack(writ12 + l) = 0.0
        do 62 l = l2, l12, l2
62      wstack(writ12 + l) = 0.0
      end if
c
c     row factor arrays for 23- and 13-planes
c
      do 63 l = 2, l3 - 1
63    vect(l) = wstack(fac3 + l - 1)
      vect(1) =  0.0
      vect(l3) = 0.0
      do 64 k = 0, lpl - 1, l3
      do 64 l = 1, l3
64    wstack(rfac23 + k + l) = vect(l)
c
c     row factor array for 12-planes if not as for 13- or 23-planes.
c
      if(l2.ne.l3) then
        do 65 l = 2, l2 - 1
65      vect(l) = wstack(fac2 + l - 1)
        vect(1)  = 0.0
        vect(l2) = 0.0
        do 66 k = 0, l12 - 1, l2
        do 66 l = 1, l2
66      wstack(rfac12 + k + l) = vect(l)
      end if
c
c     column factor arrays for 23- and 13-planes
c
      is = ctrlc
      do 4 i = 2, lpl - l3 + 2, l3
      do 67 l = 0, l3m2 - 1
67    wstack(cfac23 + i + l) = fact(is)
 4    is = is + 2
      do 50 l = 1, lpl - 1, l3
50    wstack(cfac23 + l) = 0.0
      do 277 l = l3, lpl, l3
277   wstack(cfac23 + l) = 0.0
c
c     column factor array for 12-planes
c
      if(l2.ne.l3) then
        is = ctrlc
        do 5 i = 2, l12 - l2 + 2, l2
        do 68 l = 0, l2m2 - 1
68      wstack(cfac12 + i + l) = fact(is)
 5      is = is + 2
      end if
      do 278 l = 1, l12, l2
278   wstack(cfac12 + l) = 0.0
      do 279 l = l2, l12, l2
279   wstack(cfac12 + l) = 0.0
      if(form1) then
c
c     arrays of factor contributions for planes.  the array pointed to
c     by  ar1  has a zero border, and that pointed to by  ar2  has
c     unity on the border.  arrays pointed to by  a12, a13  have a
c     zero border.
c
c     formula 1 arrays of factor contributions for planes
c
        do 69 l = 1, l23
69      vect(l) = wa1*wstack(writ23 + l)
        do 273 l = 1, l23
273     wstack(ar1 + l) = vect(l)
        do 70 l = 1, l23
70      vect(l) = (wa2*wstack(cfac23 + l) + wa3*wstack(rfac23 + l)
     1          - 1.0)*wstack(writ23 + l)
        do 71 l = 1, l23
71      wstack(ar2 + l) = vect(l) + 1.0
        if((l1.ne.l2).or.(abs(wa1 - wa2).ge.1.0e-6)) then
          do 72 l = 1, l13
72        vect(l) = wa2*wstack(writ13 + l)
          do 280 l = 1, l13
280       wstack(a13 + l) = vect(l)
        end if
        if((l2.ne.l3).or.(abs(wa2 - wa3).ge.1.0e-6)) then
          do 73 l = 1, l12
73        vect(l) = wa3*wstack(writ12 + l)
          do 281 l = 1, l12
281       wstack(a12 + l) = vect(l)
        end if
      else if(form2) then
c
c     formula 2 arrays of factor contributions for plane
c
        do 74 l = 1, l23
74      vect(l) = (wa1 - wb2*wstack(rfac23 + l)
     1                 - wb3*wstack(cfac23 + l))*wstack(writ23 + l)
        do 75 l = 1, l23
75      wstack(ar1 + l) = vect(l)
        do 76 l = 1, l23
76      vect(l) = ((wa3 - wb1*wstack(cfac23 + l))*wstack(rfac23 + l)
     1          -   1.0 +  wa2*wstack(cfac23 + l))*wstack(writ23 + l)
        do 77 l = 1, l23
77      wstack(ar2 + l) = vect(l) + 1.0
        if((.not.equal).or.(l1.ne.l2)) then
          do 78 l = 1, l13
78        vect(l) = (wa2 - wb3*wstack(rfac13 + l)
     1                   - wb1*wstack(cfac13 + l))*wstack(writ13 + l)
          do 79 l = 1, l13
79        wstack(a13 + l) = vect(l)
        end if
        if((.not.equal).or.(l2.ne.l3)) then
          do 80 l = 1, l12
80        vect(l) = (wa3 - wb1*wstack(rfac12 + l)
     1                   - wb2*wstack(cfac12 + l))*wstack(writ12 + l)
          do 81 l = 1, l12
81        wstack(a12 + l) = vect(l)
        end if
      else
c
c     formula 3 arrays of factor contributions for plane
c
        x = 0.033333333333333*wa1/wb3
        do 82 l = 1, l23
82      vect(l) = wstack(rfac23 + l)*wstack(cfac23 + l)
        do 83 l = 1, l23
83      wstack(temp + l) = vect(l)
        do 84 l = 1, l23
84      vect(l) = (wa2*wstack(cfac23 + l) + wa3*wstack(rfac23 + l)
     1          -  wb1*wstack(temp + l) - 1.0)*wstack(writ23 + l)
        do 85 l = 1, l23
85      wstack(ar2 + l) = vect(l) + 1.0
        do 86 l = 1, l23
86      vect(l) = (wa1 + wb3*(x*wstack(temp + l) - wstack(cfac23 + l))
     1          -  wb2*wstack(rfac23 + l))*wstack(writ23 + l)
        do 87 l = 1, l23
87      wstack(ar1 + l) = vect(l)
        if(l1.ne.l2) then
          do 88 l = 1, l13
88        vect(l) = (wa1 + wb1*(x*wstack(rfac13 + l) - 1.0)
     1            *  wstack(cfac13 + l) - wb1*wstack(rfac13 + l))
     2            *  wstack(writ13 + l)
          do 89 l = 1, l13
89        wstack(a13 + l) = vect(l)
        end if
        if(l2.ne.l3) then
          do 90 l = 1, l12
90        vect(l) = (wa1 + wb1*(x*wstack(rfac12 + l) - 1.0)
     1            *  wstack(cfac12 + l) - wb1*wstack(rfac12 + l))
     2            *  wstack(writ12 + l)
          do 91 l = 1, l12
91        wstack(a12 + l) = vect(l)
        end if
      end if
c
c     clear vectors for edges
c
      do 92 l = e1(1) + 1, e3(4) + l3m2
92    wstack(l) = 0.0
c
c     potential calculation - initial pointers and work area
c     assignments
c
      iplane = 1
      lp1 = ctrlc - 2
      is1 = scmbuf
      ip23 = 2
      iswpln = (l1 + 1)/2
      if(permit) ass23 = w23b
c
c     cycle over planes, stepping pointers and collecting factors.
c
      do 7 i = 1, lm2
      is1 = is1 + ll
      is = is1 + l3
      lp1 = lp1 + 2
      iplane = iplane + 1
      u = fact(lp1)
      v = fn1*fact(lp1 + 1)
c
c     adjust assignment and control if first odd-numbered plane
c
      if((iplane.eq.iswpln).and.permit) then
        ass23 = w23a
        ip23 = iswpln
      end if
c
c     calculate potentials and (if required) accumulate horizontal
c     boundary adjacent terms
c
      do 94 l = 1, l23
94    w(l + is - 1) = w(l + is - 1)/
     1                (u*wstack(ar1 + l) + wstack(ar2 + l))
      if(permit) then
c
c       border terms can be ignored here as their values are eventually
c       multiplied by the array with pointer  ar1.
c
        if(iplane.eq.ip23) then
          do 95 l = 1, l23
95        wstack(ass23 + l) = v*w(l + is - 1)
        else
          do 96 l = 1, l23
96        wstack(ass23 + l) = v*w(l + is - 1) + wstack(ass23 + l)
        end if
      end if
 7    continue
c
c     omit screening charge accumulation if not required.
c
      if(.not.permit) go to 10
c
c     boundary adjacent potentials for 13-boundaries
c
c     cycle over 13-planes
c
      lp2 = ctrlc - 1
      is = scmbuf + ll + l3
      do 8 i = 2, l2 - 1
      lp2 = lp2 + 2
      f2 = fn2*fact(lp2)
      do 97 l = 1, l13
97    work(l) = w(jstack(indx + l) + is - 1)
      if((i + i).lt.l2) then
        if(i.eq.2) then
          do 98 l = 1, l13
98        wstack(ww13b + l) = f2*work(l)
        else
          do 99 l = 1, l13
99        wstack(ww13b + l) = f2*work(l) + wstack(ww13b + l)
        end if
      else
        if((i + i).eq.(l2 + 1)) then
          do 101 l = 1, l13
101       wstack(ww13a + l) = f2*work(l)
        else
          do 102 l = 1, l13
102       wstack(ww13a + l) = f2*work(l) + wstack(ww13a + l)
        end if
      end if
 8    is = is + l3
c
c     store in main work areas
c
      do 104 l = 1, l13
104   vect(l) = (wstack(ww13a + l) + wstack(ww13b + l))
     1        * wstack(writ13 + l)
      do 105 l = 1, l13
105   wstack(w13a + l) = vect(l)
      do 106 l = 1, l13
106   vect(l) = (wstack(ww13a + l) - wstack(ww13b + l))
     1        * wstack(writ13 + l)
      do 107 l = 1, l13
107   wstack(w13b + l) = vect(l)
c
c     boundary adjacent potentials for 12-planes
c
c     cycle over planes
c
      lp3 = ctrlc - 1
      is = scmbuf + ll
      do 9 i = 2, l3 - 1
      is = is + 1
      lp3 = lp3 + 2
      f3 = fn3*fact(lp3)
      do 108 l = 1, l12
108   work(l) = w(is + l3*l - l3)
       if((i + i).lt.l3) then
        if(i.eq.2) then
          do 109 l = 1, l12
109       wstack(ww12b + l) = f3*work(l)
        else
          do 110 l = 1, l12
110       wstack(ww12b + l) = f3*work(l) + wstack(ww12b + l)
        end if
      else
        if((i + i).eq.(l3 + 1)) then
          do 112 l = 1, l12
112       wstack(ww12a + l) = f3*work(l)
        else
          do 113 l = 1, l12
113       wstack(ww12a + l) = f3*work(l) + wstack(ww12a + l)
        end if
      end if
 9    continue
c
c     accumulate into main work areas
c
      do 115 l = 1, l12
115   work(l) = (wstack(ww12a + l) + wstack(ww12b + l))
     1        * wstack(writ12 + l)
      do 116 l = 1, l12
116   wstack(w12a + l) = work(l)
      do 117 l = 1, l12
117   work(l) = (wstack(ww12a + l) - wstack(ww12b + l))
     1        * wstack(writ12 + l)
      do 118 l = 1, l12
118   wstack(w12b + l) = work(l)
c
c     calculate edge adjacent potentials if required
c
10    if(form1) go to 15
c
c     1-edge terms from 13-planes.  as pairs of boundary-adjacent
c     charges are juxtaposed in memory we may access both members
c     of a pair in a single double length gather, using the first
c     pointer only.
c
      lp3 = ctrlc + 1
      ip = 2
      lw = lm2 + lm2
      do 11 i = 2, l3 - 1
      if((i + i).eq.(l3 + 1)) ip = 1
      do 119 l = 1, lw
119   work(l) = wstack(w13a + i + l3*l - l3)
      f3 = fn3*fact(lp3)
      do 120 l = 1, lm2
120   wstack(e1(ip) + l) = f3*work(l) + wstack(e1(ip) + l)
      do 121 l = 1, lm2
121   wstack(e1(ip + 2) + l) = f3*work(l + lm2)
     1                        + wstack(e1(ip + 2) + l)
11    lp3 = lp3 + 2
      do 122 l = 1, lm2
      vect(l)       = wstack(e1(1) + l) + wstack(e1(2) + l)
122   vect(l + lm2) = wstack(e1(1) + l) - wstack(e1(2) + l)
      do 123 l = 1, lm2
123   wstack(e1(1) + l) = vect(l)
      do 124 l = 1, lm2
124   wstack(e1(2) + l) = vect(l + lm2)
      do 125 l = 1, lm2
      vect(l)       = wstack(e1(3) + l) + wstack(e1(4) + l)
125   vect(l + lm2) = wstack(e1(3) + l) - wstack(e1(4) + l)
      do 126 l = 1, lm2
126   wstack(e1(3) + l) = vect(l)
      do 127 l = 1, lm2
127   wstack(e1(4) + l) = vect(l + lm2)
c
c     3-edge terms from 13-planes
c
      ip = 3
      is = w13a + 1
      lp1 = ctrlc + 1
      do 12 i = 2, l1 - 1
      f1 = fn1*fact(lp1)
      if((i + i).eq.(l1 + 1)) ip = 1
      do 128 l = 1, l3m2
      vect(l)  = f1*wstack(l + is) + wstack(e3(ip) + l)
128   vect(l + l3m2) = f1*wstack(l + is + l13) + wstack(e3(ip + 1) + l)
      do 161 l = 1, l3m2
161   wstack(e3(ip) + l) = vect(l)
      do 129 l = 1, l3m2
129   wstack(e3(ip + 1) + l) = vect(l + l3m2)
      lp1 = lp1 + 2
12    is = is + l3
      do 130 l = 1, l3m2
      vect(l)        = wstack(e3(1) + l) + wstack(e3(3) + l)
      vect(l + l3m2) = wstack(e3(1) + l) - wstack(e3(3) + l)
      work(l)        = wstack(e3(2) + l) + wstack(e3(4) + l)
130   work(l + l3m2) = wstack(e3(2) + l) - wstack(e3(4) + l)
      do 131 l = 1, l3m2
131   wstack(e3(1) + l) = vect(l)
      do 132 l = 1, l3m2
132   wstack(e3(3) + l) = vect(l + l3m2)
      do 133 l = 1, l3m2
133   wstack(e3(2) + l) = work(l)
      do 134 l = 1, l3m2
134   wstack(e3(4) + l) = work(l + l3m2)
c
c     2-edge terms from 12-planes.  as before, we exploit the
c     juxtaposition of the areas referred to by  w12a, w12b.
c
      ip = 3
      is = w12a + 1
      lp1 = ctrlc + 1
      do 13 i = 2, l1 - 1
      f1 = fn1*fact(lp1)
      if((i + i).eq.(l1 + 1)) ip = 1
      it = is + l12
      do 135 l = 1, l2m2
      vect(l)        = f1*wstack(is + l) + wstack(e2(ip) + l)
135   vect(l + l2m2) = f1*wstack(it + l) + wstack(e2(ip + 1) + l)
      do 136 l = 1, l2m2
136   wstack(e2(ip) + l)     = vect(l)
      do 137 l = 1, l2m2
137   wstack(e2(ip + 1) + l) = vect(l + l2m2)
      lp1 = lp1 + 2
13    is = is + l2
      do 138 l = 1, l2m2
      vect(l)        = wstack(e2(1) + l) + wstack(e2(3) + l)
      vect(l + l2m2) = wstack(e2(1) + l) - wstack(e2(3) + l)
      work(l)        = wstack(e2(2) + l) + wstack(e2(4) + l)
138   work(l + l2m2) = wstack(e2(2) + l) - wstack(e2(4) + l)
      do 139 l = 1, l2m2
139   wstack(e2(1) + l) = vect(l)
      do 140 l = 1, l2m2
140   wstack(e2(3) + l) = vect(l + l2m2)
      do 141 l = 1, l2m2
141   wstack(e2(2) + l) = work(l)
      do 142 l = 1, l2m2
142   wstack(e2(4) + l) = work(l + l2m2)
c
c     corner adjacent potentials (if required) from 1-edges
c
      if(form3) then
        do 143 l = 1, (l1 - 3)/2
143     work(l) = - fn1
        do 144 l = (l1 - 1)/2, lm2
144     work(l) =   fn1
        do 145 l = 1, lm2
145     work(l) =   work(l)*wstack(sfac1 + l)
        do 14 i = 1, 4
        x = 0.0
        do 146 l = 1, lm2
        x = x + wstack(e1(i) + l)*wstack(sfac1 + l)
146     corner(i + 4) = corner(i + 4) + wstack(e1(i) + l)*work(l)
14      corner(i)     = fn1*x
      end if
c
c     assemble 23-boundaries
c
15    do 147 l = 1, ll
147   work(l)      = w(scmbuf - 1 + l)
      do 148 l = 1, ll
148   work(ll + l) = w(scmbuf + lll - ll - 1 + l)
      do 149 l = 1, l23
149   work(l3 + l)      = work(l3 + l) + wstack(ar1 + l)
     1                  *(wstack(w23a + l) + wstack(w23b + l))
      do 150 l = 1, l23
150   work(ll + l3 + l) = work(ll + l3 + l) + wstack(ar1 + l)
     1                  *(wstack(w23a + l) - wstack(w23b + l))
      if(form1) go to 19
c
c     3-edges
c
      if(form3) then
        u = 0.2*wb3
        do 151 l = 1, l3m2
151     vect(l) = u*wstack(fac3 + l)
        do 152 l = 1, l3m2
152     wstack(edge3 + l) = wb3 - vect(l)
      end if
      iad(1) = 2
      iad(3) = ll + 2
      iad(2) = iad(3) - l3
      iad(4) = iad(2) + ll
      do 16 i = 1, 4
      j = iad(i) - 1
      if(form2) then
        do 153 l = 1, l3m2
153     work(j + l) = wb3*wstack(e3(i) + l) + work(j + l)
      else
        do 154 l = 1, l3m2
154     work(j + l) = wstack(edge3 + l)*wstack(e3(i) + l)
     1              + work(j + l)
      end if
16    continue
c
c     2-edges
c
      if(form3) then
        u = 0.2*wb2
        do 155 l = 1, l2m2
155     vect(l) = u*wstack(fac2 + l)
        do 156 l = 1, l2m2
156     wstack(edge2w + l) = wb2 - vect(l)
      end if
      iad(1) = l3 + 1
      iad(2) = l3 + l3
      iad(3) = iad(1) + ll
      iad(4) = iad(2) + ll
      do 17 i = 1, 4
      j = iad(i) - l3
      do 157 l = 1, l2m2
157   vect(l) = work(j + l*l3)
      if(form2) then
        do 158 l = 1, l2m2
158     wstack(e2(i) + l) = wb2*wstack(e2(i) + l) + vect(l)
      else
        do 159 l = 1, l2m2
159     vect(l) = wstack(edge2w + l)*wstack(e2(i) + l) + vect(l)
        do 160 l = 1, l2m2
160     wstack(e2(i) + l) = vect(l)
      end if
      do 17 l = 1, l2m2
17    work(j + l*l3) = wstack(e2(i) + l)
c
c     corners if appropriate
c
      if(form3) then
        iad(1) = 1
        iad(2) = l3
        iad(3) = ll - l3 + 1
        iad(4) = ll
        u = 0.2*wb1
        do 18 i = 1, 4
        j = iad(i)
        work(j) = u*corner(i) + work(j)
18      work(j + ll) = u*corner(i + 4) + work(j + ll)
      end if
c
c     transfer assembled planes to destinations
c
19    do 162 l = 1, ll
162   wstack(bdy(1) + l) = work(l)
      do 163 l = 1, ll
163   wstack(bdy(2) + l) = work(ll + l)
c
c     assemble 13-boundaries
c
      i = ll + scmbuf
      do 164 l = 1, l13
164   work(l3 + l)        = w(jstack(indx + l) + i - 1)
      i = i + ll - l3
      do 165 l = 1, l13
165   work(ll13 + l3 + l) = w(jstack(indx + l) + i - 1)
      do 166 l = 1, l13
166   work(l3 + l)        = work(l3 + l)
     1                    + wstack(a13 + l)*wstack(w13a + l)
      do 167 l = 1, l13
167   work(ll13 + l3 + l) = work(ll13 + l3 + l)
     1                    + wstack(a13 + l)*wstack(w13b + l)
c
c     1-edge values if appropriate
c
      if(form1) go to 22
      iad(1) = l3 + 1
      iad(2) = l3 + l3
      iad(3) = iad(1) + ll13
      iad(4) = iad(2) + ll13
      if(form3) then
        u = 0.2*wb1
        do 168 l = 1, lm2
168     vect(l) = u*wstack(fac1 + l)
        do 169 l = 1, lm2
169     wstack(edge1 + l) = wb1 - vect(l)
      end if
      do 21 i = 1, 4
      j = iad(i) - l3
      do 170 l = 1, lm2
170   vect(l) = work(j + l*l3)
      if(form2) then
        do 171 l = 1, lm2
171     wstack(e1(i) + l) = wb1*wstack(e1(i) + l) + vect(l)
      else
        do 172 l = 1, lm2
172     vect(l) = wstack(edge1 + l)*wstack(e1(i) + l) + vect(l)
        do 173 l = 1, lm2
173     wstack(e1(i) + l) = vect(l)
      end if
      do 276 l = 1, lm2
276   work(j + l*l3) = wstack(e1(i) + l)
21    continue
c
c     clear 3-edges - these charges are included in the 23-planes
c
22    do 174 l = 1, l3
174   work(l) = 0.0
      do 175 l = 1, 2*l3
175   work(ll13 - l3 + l) = 0.0
      do 176 l = 1, l3
176   work(2*ll13 - l3 + l) = 0.0
c
c     transfer to destination
c
      do 177 l = 1, ll13
177   wstack(bdy(3) + l) = work(l)
      do 178 l = 1, ll13
178   wstack(bdy(4) + l) = work(ll13 + l)
c
c     assemble 12-boundaries
c
      is1 = ll12 + l2
      is = scmbuf + ll - l3
      do 179 l = 1, l12
179   work(l2 + l)       = w(is + l3*l)
      is = is + l3 - 1
      do 180 l = 1, l12
180   work(is1 + l) = w(is + l3*l)
      do 181 l = 1, l12
181   work(l2 + l) = work(l2 + l) + wstack(a12 + l)*wstack(w12a + l)
      do 182 l = 1, l12
182   work(is1 + l) = work(is1 + l) + wstack(a12 + l)*wstack(w12b + l)
c
c     clear edges - charges are in 23- or 13-planes
c
      do 183 l = 1, l2
183   work(l)             = 0.0
      do 184 l = 1, 2*l2
184   work(ll12 - l2 + l) = 0.0
      do 185 l = 1, l2
185   work(2*ll12 - l2 + l) = 0.0
      do 186 l = 1, l1 + lm2
186   work(l2*l + 1)      = 0.0
      do 187 l = 1, l1 + lm2
187   work(l2*l + l2)     = 0.0
c
c     return to destination and finish
c
      do 188 l = 1, ll12
188   wstack(bdy(5) + l) = work(l)
      do 189 l = 1, ll12
189   wstack(bdy(6) + l) = work(ll12 + l)
      return
      end
