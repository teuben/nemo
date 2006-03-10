      SUBROUTINE G6ACC(NP,GPINDX,XGP,VGP,GPACC,GPJERK,GPPOT)
*
*
*       GRAPE force emulator.
*       ---------------------
*
      INCLUDE 'common4.h'
      REAL*8  A(9),F1(3),F1DOT(3),XGP(3,NMAX),VGP(3,NMAX),
     &        GPACC(3,48),GPJERK(3,48),GPPOT(48)
      INTEGER GPINDX(48)
*
*
      DO 40 L = 1,NP
          I = GPINDX(L)
*       Initialize force & first derivative of body #I.
          DO 10 K = 1,3
              GPACC(K,L) = 0.0D0
              GPJERK(K,L) = 0.0D0
   10     CONTINUE
          GPPOT(L) = 0.0
*
*       Obtain force & first derivative by summing over all particles.
          DO 30 J = IFIRST,NTOT
              IF (J.EQ.I) GO TO 30
*
              DO 15 K = 1,3
                  A(K) = X(K,J) - XGP(K,L)
                  A(K+3) = XDOT(K,J) - VGP(K,L)
   15         CONTINUE
*
              RIJ2 = A(1)**2 + A(2)**2 + A(3)**2
              A(7) = 1.0/RIJ2
              A(8) = BODY(J)*A(7)*SQRT(A(7))
              A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
              DO 20 K = 1,3
                  F1(K) = A(K)*A(8)
                  F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
   20         CONTINUE
*
              DO 25 K = 1,3
                  GPACC(K,L) = GPACC(K,L) + F1(K)
                  GPJERK(K,L) = GPJERK(K,L) + F1DOT(K)
   25         CONTINUE
              GPPOT(L) = GPPOT(L) - BODY(J)/SQRT(RIJ2)
   30     CONTINUE
*     WRITE (6,35) I, L, (GPACC(K,L),K=1,3)
*  35 FORMAT (' G6ACC:   I L F  ',2I5,1P,3E10.2)
   40 CONTINUE
*
      RETURN
*
      END
