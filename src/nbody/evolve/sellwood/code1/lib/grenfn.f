      subroutine grenfn
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        this routine is a dummy version leading to an abort if
c        called.  it is intended to be replaced by a routine to
c        evaluate an analytical form for the near field greens
c        function if this is ever required.
c
      call crash( 'grenfn', 'uf 1' )
      stop
      end
