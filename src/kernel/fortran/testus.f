C
C  test to see how symbol names get mangled by fortran compilers
C  Not recommended
C
      SUBROUTINE a_b
      WRITE(*,*) 'Running a_b'
      END
C
      SUBROUTINE ab
      WRITE(*,*) 'Running ab'
      END
