C
C:      Test program for NEMO's footran interface
C+
C   x=4.0\n             number to test
C   VERSION=1.1\n       15-mar-01 PJT
C-

        SUBROUTINE NEMOMAIN  
        INCLUDE 'getparam.inc'
        DOUBLE PRECISION x, funcf

        x = getdparam('x')
         
        WRITE(*,*) 'funcf(x) = ',funcf(x)

        RETURN
        END


