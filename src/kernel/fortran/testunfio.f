      PROGRAM unfio
c
c     write some unformatted fortran I/O
c
      INTEGER i,n,iarr(10)
      REAL rarr(10)
      DOUBLE PRECISION darr(10)

      n=10
      DO i=1,n
         iarr(i) = i
         rarr(i) = REAL(i)
         darr(i) = DBLE(i)
      ENDDO

      OPEN(9,FILE='tmp.o',FORM='UNFORMATTED')
      WRITE(9) n
      WRITE(9) iarr
      WRITE(9) rarr
      WRITE(9) darr
      CLOSE(9)

      WRITE(*,*) 'A file tmp.o has been created'
      WRITE(*,*) 'To test, run:'
      WRITE(*,*) '   unfio tmp.o'

      END
