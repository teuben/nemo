      SUBROUTINE FCLOUD(I,FI,FD)
*
*
*       Force & first derivative from interstellar clouds.
*       --------------------------------------------------
*
      INCLUDE 'common4.h'
      COMMON/CLOUDS/  XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &                CLM(MCL),CLMDOT(MCL),VCL,SIGMA,RB2,RB3,PCL2,
     &                TCL,STEPCL,TBIG,DTBIG,NCL,NEWCL
      REAL*8  FI(3),FD(3),A(6)
*
*
*       Sum over all clouds.
      DO 20 ICL = 1,NCL
          DO 5 K = 1,3
              A(K) = XCL(K,ICL) - X(K,I)
              A(K+3) = XDOTCL(K,ICL) - XDOT(K,I)
    5     CONTINUE
*
          RIJ2 = A(1)**2 + A(2)**2 + A(3)**2 + RCL2(ICL)
          A5 = BODYCL(ICL)/(RIJ2*SQRT(RIJ2))
          A6 = BODYCL(ICL)*RB3
          A7 = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))/RIJ2
*
*       Add cloud force & first derivative and subtract average field.
          DO 10 K = 1,3
              FI(K) = FI(K) + A(K)*A5 + (X(K,I) - RDENS(K))*A6
              FD(K) = FD(K) + (A(K+3) - A(K)*A7)*A5 + XDOT(K,I)*A6
   10     CONTINUE
   20 CONTINUE
*
      RETURN
*
      END
