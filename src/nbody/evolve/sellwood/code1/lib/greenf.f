      subroutine greenf(analyt)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to calculate a source function as preparation for finding
c        the greens function  over a small mesh.
c
      include 'rjinclude.h'
c
c local variables
      integer i
c
      logical analyt
      integer p
c
c     report entry, jump if greens function not prescribed analytically
c
      write(s2, 200)
200   format('0calculation of greens function')
      if(.not.analyt) go to 42
      call grenfn
      return
c
c     clear mesh, apart from a unit charge at the centre
c
42    p = scmbuf + n1*n2*n3 - 1
      do 1 i = scmbuf, p
1     w(i) = 0.0
      p = (((n1 - 1)/2)*n2 + (n2 - 1)/2)*n3 + (n3 - 1)/2 + scmbuf
      w(p) = 1.0
      return
      end
