      subroutine markrs
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to calculate pointer values for the potential solver
c
      include 'rjinclude.h'
      integer n(3)
      equivalence (n(1), n1)
c
c     initialise
c
      opt3 = .false.
      copy = .false.
      b(1) = 1
      b(2) = 1
      b(3) = 1
c
c     set up three dimensional case
c
      incrl = n(2)
c
c     limits on records transferred
c
      nom(1) = n(2)*n(3)
      nom(2) = n(3)*n(1)
      nom(3) = n(1)*n(2)
c
c     wafer normals
c
      norm(1) = 2
      norm(2) = 1
      norm(3) = 1
c
c     limits of scans within prism
c
      limit(1) = n(3)
      limit(2) = n(3)
      limit(3) = n(2)*n(1)
c
c     markers for edge of mesh
c
      mark(1) = b(3)
      mark(2) = b(3)
      mark(3) = b(2)
c
c     excesses at edge
c
      edge(1) = 1
      edge(2) = edge(1)
      edge(3) = 1
c
c     plane advance increment and termination
c
      pl(1, 1) = n(2)
      pl(1, 2) = 1
      pl(1, 3) = n(2)
      pl(1, 4) = n(2)
      pl(2, 1) = n(1)*n(2)
      pl(2, 2) = n(2)
      pl(2, 3) = n(2)*n(1)
      pl(2, 4) = n(2)*n(1)
      pl(3, 1) = n(1)*n(2)
      pl(3, 2) = n(2)
      pl(3, 3) = n(2)*n(1)
      pl(3, 4) = n(2)*(n(1) - 1)
c
c     suppress check for last few planes
c
      edpl(1, 1) = 2
      edpl(2, 1) = 2
      edpl(3, 1) = 2
c
c     adjust terminators and limits
c
      limit(1) = n(3)
      limit(2) = n(3)
      limit(3) = n(2)
      pl(1, 1) = n(2)
      pl(2, 1) = n(1)*n(2)
      pl(3, 1) = pl(2, 1)
c
c     check parameters for last few planes
c
      edpl(1, 1) = (b(2) - 1)*b(3)
      edpl(2, 1) = (b(1) - 1)*b(3)
      edpl(3, 1) = (b(1) - 1)*b(2)
      edpl(1, 2) = 1
      edpl(2, 2) = n(2)
      edpl(3, 2) = edpl(2, 2)
c
c     set up factors for formula 2
c
      w1 = 0.083333333333*(wt2 + wt3)
      w2 = 0.083333333333*(wt3 + wt1)
      w3 = 0.083333333333*(wt1 + wt2)
c
c     initialise charge factors
c
      cf(1) = wt1
      cf(2) = wt2
      cf(3) = wt3
      cf(4) = w1
      cf(5) = w2
      cf(6) = w3
      if(formul.le.2) go to 3
c
c     check geometry
c
      if(((wt1 - wt2)**2 + (wt1 - wt3)**2 + (wt2 - wt3)**2).lt.1.0e-12)
     1 go to 3
      write(s2, 201) formul, wt1, wt2, wt3
201   format('0cell geometry not consistent with formula', i10
     1/'0parameters are', 3f20.7)
      stop
c
c     fault print if insufficient store
c
3     if((nom(1).ge.1).and.(nom(2).ge.1).and.(nom(3).ge.1)) return
      write(s2, 200) scmbuf, lnscm, n
200   format('0insufficient store'/'0work space from', i10, 5x
     1, 'to', i10/ '0mesh parameters are', 3i10)
      stop
      end
