      real function grenf3(i, j, k)
c routine written by Richard James for 3-D Poisson solver
      integer i, j, k
c
c     purpose
c
c        to evaluate the Newtonian approximation to the
c        greens function
c
      include 'rjinclude.h'
c
c     evaluate Newtonian expression for distant potential
c
      real x
c
      x = float(i)**2/wt1 + float(j)**2/wt2 + float(k)**2/wt3
      grenf3 = 1.0/sqrt(x)
      return
      end
