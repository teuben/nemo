      subroutine vbdgrn(bdy)
c routine written by Richard James for 3-D Poisson solver
      integer bdy( 6 )
c
c     purpose
c
c        to calculate the greens function on the mesh boundary
c        from the analytical form given by  grenf3
c
c    parameter:
c
c     bdy - an integer array holding pointers to vectors to receive
c           the results.
c
c    note:   no attempt has been made to vectorise this routine.
c
      include 'rjinclude.h'
c
c external
      real grenf3
c
c local variables
      integer i, j, k
      real x, y
c
      integer p, lp1, lp2, lp3, lq1, lq2, lq3
c, bdy(6)
      logical yes
c
c     check mesh symmetry
c
      if((n1.ne.n2).or.(n1.ne.n3)) then
        write( s2, 201 )n1, n2, n3
 201    format( 'A symmetrical mesh is required for vbdgrn(xi)' /
     1          ' the actual dimensions are', 3i10 )
        stop
      end if
c
c     initial correction factor and initial points
c
      factor = 0.0
      lp1 = (n1 + 1)/2
      lp2 = (n2 + 1)/2
      lp3 = (n3 + 1)/2
c
c     values on 1 - boundaries
c
      p = 1
      do 2 j = 1, n2
      yes = j.le.lp2
      y = 8.0*wt1
      if(j.eq.lp2) y = 4.0*wt1
      do 2 k = 1, n3
      if(k.eq.lp3) y = 0.5*y
      x = grenf3(lp1 - 1, j - lp2, k - lp3)
      if(yes.and.(k.le.lp3))
     1 factor = y*(x - grenf3(lp1, j - lp2, k - lp3)) + factor
      w(p) = x
2     p = p + 1
c
c     adjust scaling factor for 3 pairs of boundaries
c
      factor = 3.0*factor
      if(formul.eq.1) go to 9
c
c     edge contribution
c
      y = - 0.666666666667*(wt2 + wt3)
      y = 3.0*y
      lq1 = 1 - lp1
      lq2 = 1 - lp2
      lq3 = 1 - lp3
      do 12 i = 1, lp1
      if(i.eq.lp1) y = 0.5*y
12    factor = y*(grenf3(i - lp1, lq2, lq3)
     1 + grenf3(i - lp1, lq2 - 1, lq3 - 1) - grenf3(i - lp1,lq2,lq3-1)
     2 - grenf3(i - lp1, lq2 - 1, lq3)) + factor
      if(formul.eq.2) go to 9
c
c     corner contribution
c
      lq1 = - lq1
      lq2 = - lq2
      lq3 = - lq3
      factor = - 0.0444444444444*(grenf3(lp1, lp2, lp3)
     1 - grenf3(lq1, lp2, lp3) - grenf3(lp1, lq2, lp3)
     2 - grenf3(lp1, lp2, lq3) + grenf3(lp1, lq2, lq3)
     3 + grenf3(lq1, lp2, lq3) + grenf3(lq1, lq2, lp3)
     4 - grenf3(lq1, lq2, lq3)) + factor
9     factor = 1.0/factor
      write(s2, 200) factor
200   format('0normalising factor for greens function =', f20.10)
c
c     sine transform boundaries
c
      p = 1
      lq1 = n2
      lq2 = n3
      lq3 = p + lq1*lq2 - 1
      do 4 j = p, lq3
4     w(j) = factor*w(j)
      call dtrfm1(p, p, lq1, lq2, .false., .false.)
      call dtrfm3(p, p, lq1, lq2, .false., .false.)
      lq3 = bd(i + 1)
13    p = p + lq1*lq2
c
c     transfer to result area in stack.
c
      do 6 i = 1, 6
      do 6 j = 1, n2*n3
 6    wstack(bdy(i) + j) = w(j)
      return
      end
