c
c	Example of a Fortran program
c
c       On *Unix* to be compiled and linked as:
c       
c       f77 -g -C -u -o mainf mainf.f
c     
c       The flags '-g -C -u' are all defensive programming flags
c       and are optional
c
      PROGRAM example
c
      INTEGER   NMAX
      PARAMETER (NMAX=10)
c
      REAL    a
      INTEGER i

      a = 1.0
      DO 100 i=1,NMAX
         a = a + a
  100 CONTINUE
   
      WRITE(*,'(A,G20.10)') 'The sum is ',a

      END
