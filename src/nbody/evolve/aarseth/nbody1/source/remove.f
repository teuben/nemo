      SUBROUTINE REMOVE(I)
*
*
*       Removal of particle.
*       --------------------
*
      INCLUDE 'common1.h'
      REAL*4  A(6)
*
*
*       Correct forces & first derivatives.
      DO 10 J = 1,N
          IF (J.EQ.I) GO TO 10
          RIJ2 = 0.0
          RIJDOT = 0.0
*
          DO 1 K = 1,3
              A(K) = X(K,I) - X(K,J)
              A(K+3) = XDOT(K,I) - XDOT(K,J)
              RIJ2 = RIJ2 + A(K)**2
              RIJDOT = RIJDOT + A(K)*A(K+3)
    1     CONTINUE
*
          RIJ2 = RIJ2 + EPS2
          A5 = BODY(I)/(RIJ2*SQRT(RIJ2))
          A6 = 3.0*RIJDOT/RIJ2
*
*       Subtract force & first derivative (first order) due to escaper.
          DT = TIME - T0(J)
          DO 5 K = 1,3
              DF = A(K)*A5
              DD = (A(K+3) - A6*A(K))*A5
              F(K,J) = F(K,J) - 0.5*DF
              FDOT(K,J) = FDOT(K,J) - ONE6*DD
              D1(K,J) = D1(K,J) - DD
*       Remove contribution to X & XDOT for consistency (small effect).
              X(K,J) = X(K,J) - (ONE6*DD*DT + 0.5*DF)*DT**2
              XDOT(K,J) = XDOT(K,J) - (0.5*DD*DT + DF)*DT
    5     CONTINUE
   10 CONTINUE
*
*       Skip table updates for last body.
      IF (I.GT.N) GO TO 30
*
*       Move up all COMMON variables.
      DO 20 J = I,N
          J1 = J + 1
          DO 15 K = 1,3
              X(K,J) = X(K,J1)
              X0(K,J) = X0(K,J1)
              X0DOT(K,J) = X0DOT(K,J1)
              XDOT(K,J) = XDOT(K,J1)
              F(K,J) = F(K,J1)
              FDOT(K,J) = FDOT(K,J1)
              D1(K,J) = D1(K,J1)
              D2(K,J) = D2(K,J1)
              D3(K,J) = D3(K,J1)
   15     CONTINUE
*
          BODY(J) = BODY(J1)
          NAME(J) = NAME(J1)
          STEP(J) = STEP(J1)
          T0(J) = T0(J1)
          T1(J) = T1(J1)
          T2(J) = T2(J1)
          T3(J) = T3(J1)
   20 CONTINUE
*
   30 RETURN
*
      END
