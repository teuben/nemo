      SUBROUTINE DFSP(I,I1,FIRR,FD)
*
*
*       Force on single particle due to KS pair.
*       ----------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  FIRR(3),FD(3),DX(3),DV(3)
*
*
*       Add contributions from each component or just the c.m.
      K = I1
    1 dr2 = 0.0
      drdv = 0.0
      DO 5 L = 1,3
          dx(L) = X(L,K) - X(L,I)
          dv(L) = XDOT(L,K) - XDOT(L,I)
          dr2 = dr2 + dx(L)**2
          drdv = drdv + dx(L)*dv(L)
    5 CONTINUE
*
      dr2i = 1.0/dr2
      dr3i = BODY(K)*dr2i*SQRT(dr2i)
      drdv = 3.0*drdv*dr2i
*
      DO 10 L = 1,3
          FIRR(L) = FIRR(L) + dx(L)*dr3i
          FD(L) = FD(L) + (dv(L) - dx(L)*drdv)*dr3i
   10 CONTINUE
*
*       Treat the second component similarly (unless I1 > N).
      IF (I1.LE.N) THEN
          IF (K.EQ.I1) THEN
              K = I1 + 1
              GO TO 1
          END IF
      END IF
*
      RETURN
*
      END
