      subroutine slvptt(bdy)
c routine written by Richard James for 3-D Poisson solver
      integer bdy( 6 )
c
c     purpose
c
c        to solve the poisson equation in transform space.  when not
c        calculating a greens function the routine supplements the
c        potentials by terms arising from the boundary potentials
c        and copies the boundary potentials back to the mesh.
c
c     parameter
c
c     bdy   -  an integer array pointing to vectors holding the
c              boundary values.
c
c     note:
c
c         the routine operates on the mesh after the two mesh rotations
c         required in the fourier transform stage.  the mesh therefore
c         contains n3 planes, each with n1 rows of n2 elements.
c
      include 'rjinclude.h'
c
c local variables
      integer i, iplane, is, isp, is0, is1, is23, j, jj, j0, j1, j2, j3,
     1      j4, k, l, ll, lll, ll12, ll13, lm2, lpl, lsw, l1, l12, l13,
     2      l2, l2m2, l23, l3, l3m2
      real fn1, fn2, fn3, u, v, wa1, wa2, wa3, wb1, wb2, wb3, x
c
      integer lp1, lp2, lp3, iad(4)
c      integer bdy(6)
      integer writ23, writ13, writ12, rfac23, rfac13, rfac12,
     1        cfac23, cfac13, cfac12, ar1, ar2, temp, e1(4), e2(4),
     2        e3(4), edge1, edge2, edge2w, edge3, indx,
     3        bdyc(2), fac1, fac2, fac3, sfac1, sfac3,
     4        a12, a13, chg13a, chg13b, chg12a, chg12b
      logical form1, form2, form3, equal
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
      fn1 = 2.0
      fn2 = 2.0
      fn3 = 2.0
c
c     check plane sizes
c
      i = 2*max0(ll, ll13, ll12)
      if(i.gt.lnwork) then
        write( s2, 900 )ll, ll13, ll12, i, lnwork
 900    format( ' area  /work/  too small in vbdpot(xe)' /
     1          ' plane sizes are', 3i10, 5x, 'length required =',
     2          i10 / ' actual length =', i10 )
        call crash( 'slvptt', 'xe 1')
      end if
c
c     check size of vectorisation work area.
c
      if(max0(ll, ll13, ll12).gt.lnvect) then
        write( s2, 901 )ll, ll13, ll12, lnvect
 901    format( ' area /vect/ smaller than a boundary plane' /
     1          ' plane sizes are', 3i10, ' limit =', i10 )
        call crash( 'slvptt', 'xe 2')
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
c
c     arrays for equivalent charges.
c
      bdyc(1)  = l
      bdyc(2)  = l + l23
c
c     arrays for 13 charge planes.
c
      chg13a   = bdyc(2) + l23
      chg13b   = chg13a  + l13
      chg12a   = chg13b  + l13
      chg12b   = chg12a  + l12
      l        = chg12b  + l12
      if(l.gt.maxstk) then
        maxstk = l
        mxstid = 'slvptt'
      end if
c
c     check stack space.
c
      if(l.ge.lstack) then
        write( s2, 902 )istack, l, lstack
 902    format( ' stack overflow in vslvpt(xe)' /
     1          ' stack pointer =', i10, ' last pointer =', i10,
     2          ' limit =', i10 )
        call crash( 'slvptt', 'xe 3')
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
c     supplementation of the potential by boundary terms
c
c     first transfer the potentials for the lower horizontal
c     boundary to the mesh
c
      do 190 l = 1, ll
190   w(scmbuf - 1 + l) = wstack(bdy(1) + l)
c
c     assemble equivalent charge layers for 23-planes.  these contain
c     all corner contributions and all the 2- and 3-edge terms.
c
      do 191 l = 1, ll
191   work(l)      = wstack(bdy(1) + l)
      do 192 l = 1, ll
192   work(ll + l) = wstack(bdy(2) + l)
      do 282 l = 1, l23
282   vect(l) = wstack(ar1 + l) + 1.0 - wstack(writ23 + l)
      do 193 l = 1, l23
193   work(l3 + l)      = vect(l)*work(l3 + l)
      do 194 l = 1, l23
194   work(ll + l3 + l) = vect(l)*work(ll + l3 + l)
c
c     omit edge calculation for simplest formula
c
      if(form1) go to 29
c
c     deal with corners if necessary.  corner terms are spread
c     along the 3-edges.
c
      j0 = ll + ll
      if(form3) then
        iad(1) = 1
        iad(2) = l3
        iad(3) = ll - l3 + 1
        iad(4) = ll
        u = 0.2*wb1*fn3
        do 25 i = 1, 4
        j = iad(i)
        corner(i) = work(j)
25      corner(i + 4) = work(j + ll)
        k = (l3 - 3)/2
        do 195 l = 1, k
195     work(j0 + l)     = - u
        do 196 l = 1, k + 1
196     work(j0 + k + l) =   u
        do 197 l = 1, l3m2
197     work(j0 + l) = wstack(sfac3 + l)*work(j0 + l)
        v = 0.2*wb3
        do 198 l = 1, l3m2
198     vect(l) = v*wstack(fac3 + l)
        do 199 l = 1, l3m2
199     wstack(edge3 + l) = wb3 - vect(l)
      else
        do 200 l = 1, l3m2
200     wstack(edge3 + l) = wb3
      end if
c
c     3-edge equivalent charges
c
      iad(1) = 2
      iad(3) = ll + 2
      iad(2) = iad(3) - l3
      iad(4) = iad(2) + ll
      do 26 i = 1, 4
      jj = i + i
      j = iad(i)
      do 201 l = 1, l3m2
201   vect(l) = wstack(edge3 + l)*work(j - 1 + l)
      do 202 l = 1, l3m2
202   wstack(e3(i) + l) = vect(l)
      if(form3) then
        do 203 l = 1, l3m2
203     vect(l) = (u*corner(jj - 1))*wstack(sfac3 + l)
        do 204 l = 1, l3m2
204     vect(l) = vect(l) + corner(jj)*work(j0 + l)
        do 114 l = 1, l3m2
114     wstack(e3(i) + l) = vect(l) + wstack(e3(i) + l)
      end if
26    continue
c
c     sum and difference in 2-direction
c
      do 205 l = 1, l3m2
205   work(j0 + l) = wstack(e3(1) + l) + wstack(e3(2) + l)
      do 206 l = 1, l3m2
206   vect(l)      = wstack(e3(1) + l) - wstack(e3(2) + l)
      do 207 l = 1, l3m2
207   wstack(e3(1) + l) = vect(l)
      do 208 l = 1, l3m2
208   wstack(e3(2) + l) = work(j0 + l)
      do 209 l = 1, l3m2
209   work(j0 + l) = wstack(e3(3) + l) + wstack(e3(4) + l)
      do 210 l = 1, l3m2
210   vect(l)      = wstack(e3(3) + l) - wstack(e3(4) + l)
      do 211 l = 1, l3m2
211   wstack(e3(3) + l) = vect(l)
      do 212 l = 1, l3m2
212   wstack(e3(4) + l) = work(j0 + l)
c
c     2-edge equivalent charges
c
      iad(1) = l3 + 1
      iad(2) = l3 + l3
      iad(3) = ll + iad(1)
      iad(4) = ll + iad(2)
      if(form3) then
        u = 0.2*wb2
        do 213 l = 1, l2m2
213     vect(l) = u*wstack(fac2 + l)
        do 214 l = 1, l2m2
214     wstack(edge2 + l) = wb2 - vect(l)
      else
        do 215 l = 1, l2m2
215     wstack(edge2 + l) = wb2
      end if
      do 27 i = 1, 4
      j = iad(i) - l3
      do 216 l = 1, l2m2
216   vect(l) = work(j + l3*l)*wstack(edge2 + l)
      do 217 l = 1, l2m2
217   wstack(e2(i) + l) = vect(l)
27    continue
c
c     sum and difference in 3-direction
c
      j1 = j0 + l3m2
      j2 = j1 + l2m2
      j3 = j2 + l2m2
      j4 = j3 + l2m2
      do 218 l = 1, l2m2
218   work(j1 + l) = fn3*(wstack(e2(1) + l) - wstack(e2(2) + l))
      do 219 l = 1, l2m2
219   work(j2 + l) = fn3*(wstack(e2(1) + l) + wstack(e2(2) + l))
      do 220 l = 1, l2m2
220   work(j3 + l) = fn3*(wstack(e2(3) + l) - wstack(e2(4) + l))
      do 221 l = 1, l2m2
221   work(j4 + l) = fn3*(wstack(e2(3) + l) + wstack(e2(4) + l))
c
c     spread edge contributions over mesh
c
      lsw = (l3 - 3)/2
      lp2 = ctrlc + 1
      is = l3 + 2
      do 28 j = 2, l2 - 1
      is1 = is + ll
      u = fn2*fact(lp2)
      lp2 = lp2 + 2
      if((j + j).lt.l2) then
        do 222 l = 1, l3m2
222     work(is  - 1 + l) = work(is  - 1 + l) + u*wstack(e3(1) + l)
        do 223 l = 1, l3m2
223     work(is1 - 1 + l) = work(is1 - 1 + l) + u*wstack(e3(3) + l)
      else
        do 224 l = 1, l3m2
224     work(is  - 1 + l) = work(is  - 1 + l) + u*wstack(e3(2) + l)
        do 225 l = 1, l3m2
225     work(is1 - 1 + l) = work(is1 - 1 + l) + u*wstack(e3(4) + l)
      end if
c
c     contribution vector for 2-edge terms
c
      do 226 l = 1, lsw
226   vect(l) = work(j1 + j - 1)
      do 227 l = lsw + 1, l3m2
227   vect(l) = work(j2 + j - 1)
      do 228 l = 1, lsw
228   vect(l3m2 + l) = work(j3 + j - 1)
      do 229 l = lsw + 1, l3m2
229   vect(l3m2 + l) = work(j4 + j - 1)
      do 230 l = 1, l3m2
230   work(is  - 1 + l) = work(is  - 1 + l) + wstack(sfac3 + l)
     1                                      * vect(l)
      do 231 l = 1, l3m2
231   work(is1 - 1 + l) = work(is1 - 1 + l) + wstack(sfac3 + l)
     1                                      * vect(l3m2 + l)
28    is = is + l3
c
c     preserve equivalent charges on horizontal boundaries
c
29    do 232 l = 1, l23
232   wstack(bdyc(1) + l) = work(l3 + l) - work(ll + l3 + l)
      do 233 l = 1, l23
233   wstack(bdyc(2) + l) = work(l3 + l) + work(ll + l3 + l)
c
c     equivalent charges on 13-boundaries.  only values selected by
c     the mask with pointer  writ13  are changed.
c
      do 234 l = 1, ll13
234   work(l)        = wstack(bdy(3) + l)
      do 235 l = 1, ll13
235   work(ll13 + l) = wstack(bdy(4) + l)
      do 93 l = 1, l13
93    vect(l) = wstack(a13 + l) + 1.0 - wstack(writ13 + l)
      do 236 l = 1, l13
236   work(l3 + l)        = vect(l)*work(l3 + l)
      do 237 l = 1, l13
237   work(ll13 + l3 + l) = vect(l)*work(ll13 + l3 + l)
c
c     1-edge terms - this is the last set to be included
c
c     omit code if simplest formula
c
      if(form1) go to 32
c
c     collect edge values
c
      iad(1) = 1
      iad(2) = l3
      iad(3) = ll13 + iad(1)
      iad(4) = ll13 + iad(2)
      if(form3) then
        u = 0.2*wb1
        do 238 l = 1, lm2
238     vect(l) = wb1 - u*wstack(fac1 + l)
      else
        do 239 l = 1, lm2
239     vect(l) = wb1
      end if
      do 30 i = 1, 4
      j = iad(i)
      do 240 l = 1, lm2
240   wstack(e1(i) + l) = vect(l)*work(j + l3*l)
30    continue
c
c     sums and differences in 3-direction
c
      j1 = ll13 + ll13
      j2 = j1 + lm2
      j3 = j2 + lm2
      j4 = j3 + lm2
      do 241 l = 1, lm2
241   work(j1 + l) = fn3*(wstack(e1(1) + l) - wstack(e1(2) + l))
      do 242 l = 1, lm2
242   work(j2 + l) = fn3*(wstack(e1(1) + l) + wstack(e1(2) + l))
      do 243 l = 1, lm2
243   work(j3 + l) = fn3*(wstack(e1(3) + l) - wstack(e1(4) + l))
      do 244 l = 1, lm2
244   work(j4 + l) = fn3*(wstack(e1(3) + l) + wstack(e1(4) + l))
c
c     cycle over lines
c
      is0 = l3 + 2
      is1 = is0 + ll13
      do 31 i = 2, l1 - 1
      do 245 l = 1, lsw
245   vect(l)        = work(j1 + i - 1)
      do 246 l = lsw + 1, l3m2
246   vect(l)        = work(j2 + i - 1)
      do 247 l = 1, lsw
247   vect(l3m2 + l) = work(j3 + i - 1)
      do 248 l = lsw + 1, l3m2
248   vect(l3m2 + l) = work(j4 + i - 1)
      do 249 l = 1, l3m2
249   work(is0 - 1 + l) = work(is0 - 1 + l) + wstack(sfac3 + l)
     1                                      * vect(l)
      do 250 l = 1, l3m2
250   work(is1 - 1 + l) = work(is1 - 1 + l) + wstack(sfac3 + l)
     1                                      * vect(l3m2 + l)
      is0 = is0 + l3
31    is1 = is1 + l3
c
c     preserve equivalent charges for 13-boundaries
c
32    is0 = l3
      is1 = is0 + ll13
      do 251 l = 1, l13
251   wstack(chg13a + l) = work(is0 + l) - work(is1 + l)
      do 252 l = 1, l13
252   wstack(chg13b + l) = work(is0 + l) + work(is1 + l)
c
c     equivalent charges on 12-boundaries.  edge terms do not enter
c     here.  only values selected by the mask with pointer  writ12
c     are multiplied by array values with pointer  a12.
c
      do 253 l = 1, ll12
253   work(l)        = wstack(bdy(5) + l) - wstack(bdy(6) + l)
      do 254 l = 1, ll12
254   work(ll12 + l) = wstack(bdy(5) + l) + wstack(bdy(6) + l)
      do 100 l = 1, l12
100   vect(l) = wstack(a12 + l) + 1.0 - wstack(writ12 + l)
      do 103 l = 1, l12
103   work(l2 + l)        = vect(l)*work(l2 + l)
      do 111 l = 1, l12
111   work(ll12 + l2 + l) = vect(l)*work(ll12 + l2 + l)
c
c     store equivalent charges
c
      is0 = l2
      is1 = ll12 + is0
      do 255 l = 1, l12
255   wstack(chg12a + l) = work(is0 + l)
      do 256 l = 1, l12
256   wstack(chg12b + l) = work(is1 + l)
c
c     now ready to down 3ment and smooth interior potentials
c
      isp = scmbuf + ll + l3
      is23 = 1
      iplane = 2
      lp1 = ctrlc - 2
c
c     transfer potentials on vertical boundaries to mesh
c
      i = isp - 2*l3
      do 257 l = 1, l12
257   w(i + l3*l) = wstack(bdy(5) + l2 + l)
      i = i + l3 - 1
      do 258 l = 1, l12
258   w(i + l3*l) = wstack(bdy(6) + l2 + l)
      i = isp - l3 - 1
      do 259 l = 1, l13
259   w(i + jstack(indx + l)) = wstack(bdy(3) + l3 + l)
      i = i + ll - l3
      do 260 l = 1, l13
260   w(i + jstack(indx + l)) = wstack(bdy(4) + l3 + l)
c
c     cycle over 23-planes, restoring original charges.  if not
c     calculating a greens function, decrease the charges to give
c     correct smoothing when the potentials are calculated, and
c     also supplement from the 23-boundaries at this stage.
c
      do 35 i = 1, lm2
      lp1 = lp1 + 2
      if((iplane + iplane).eq.(l1 + 1)) is23 = 2
      iplane = iplane + 1
      u = fact(lp1)
      v = fn1*fact(lp1 + 1)
      if(gren) then
        do 261 l = 1, l23
261     w(isp - 1 + l) = w(isp - 1 + l) + v*wstack(writ23 + l)
     1                                     *wstack(bdyc(is23) + l)
      else
        do 262 l = 1, l23
262     vect(l) = u*wstack(ar1 + l) + wstack(ar2 + l)
        do 263 l = 1, l23
263     vect(l) = vect(l)*(1.0 - smooth*vect(l))*w(isp - 1 + l)
     1          + v*wstack(bdyc(is23) + l)
        do 264 l = 1, l23
264     w(isp - 1 + l) = w(isp - 1 + l)*(1.0 - wstack(writ23 + l))
     1                 + vect(l)*wstack(writ23 + l)
      end if
35    isp = isp + ll
c
c     supplement from 13-boundaries
c
      lp2 = ctrlc - 1
      isp = isp - lm2*ll
      is = isp
      do 36 i = 2, l2 - 1
      do 265 l = 1, l13
265   work(l) = w(is - 1 + jstack(indx + l))
      lp2 = lp2 + 2
      u = fn2*fact(lp2)
      if((i + i).lt.l2) then
        do 266 l = 1, l13
266     work(l) = work(l) + u*wstack(writ13 + l)*wstack(chg13a + l)
      else
        do 267 l = 1, l13
267     work(l) = work(l) + u*wstack(writ13 + l)*wstack(chg13b + l)
      end if
      do 268 l = 1, l13
268   w(is - 1 + jstack(indx + l)) = work(l)
36    is = is + l3
c
c     supplement from 12-boundaries
c
      lp3 = ctrlc - 1
      is = isp - l3 + 1
      do 37 i = 2, l3 - 1
      do 269 l = 1, l12
269   work(l) = w(is - l3 + l3*l)
      lp3 = lp3 + 2
      u = fn3*fact(lp3)
      if((i + i).lt.l3) then
        do 270 l = 1, l12
270     work(l) = work(l) + u*wstack(writ12 + l)*wstack(chg12a + l)
      else
        do 271 l = 1, l12
271     work(l) = work(l) + u*wstack(writ12 + l)*wstack(chg12b + l)
      end if
      do 272 l = 1, l12
272   w(is - l3 + l3*l) = work(l)
37    is = is + 1
c
c     find new potentials
c
      lp1 = lp1 - 2*lm2
      do 38 i = 1, lm2
      lp1 = lp1 + 2
      u = fact(lp1)
      do 274 l = 1, l23
274   w(isp - 1 + l) = w(isp - 1 + l)/
     1                 (u*wstack(ar1 + l) + wstack(ar2 + l))
38    isp = isp + ll
c
c     transfer potential for top horizontal boundary and finish
c
      do 275 l = 1, ll
275   w(scmbuf + lll - ll - 1 + l) = wstack(bdy(2) + l)
      return
      end
