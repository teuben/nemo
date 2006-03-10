      SUBROUTINE DFSP2(I,J1,FIRR,FD)
*
*
*       Force correction on single particle due to KS pair.
*       ---------------------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  FIRR(3),FD(3),DX(3),DV(3)
*
*
*       Define pair index and obtain c.m. distance. 
      JPAIR = KVEC(J1)
      J = N + JPAIR
      dr2 = 0.0
      DO 2 L = 1,3
          dx(L) = X(L,J) - X(L,I)
          dr2 = dr2 + dx(L)**2
    2 CONTINUE
*
*       Skip correction if c.m. approximation is satisfied.
      IF (dr2.GT.CMSEP2*R(JPAIR)**2) GO TO 20
*
*       Ensure prediction of neighbour for new polynomial.
      IF (IPHASE.NE.0) THEN
          CALL XVPRED(J,-2)
      END IF
*
*       Subtract the c.m. terms included on GRAPE.
      drdv = 0.0
      DO 3 L = 1,3
          dv(L) = XDOT(L,J) - XDOT(L,I)
          drdv = drdv + dx(L)*dv(L)
    3 CONTINUE
*
      dr2i = 1.0/dr2
      dr3i = BODY(J)*dr2i*SQRT(dr2i)
      drdv = 3.0*drdv*dr2i
*
      DO 4 L = 1,3
          FIRR(L) = FIRR(L) - dx(L)*dr3i
          FD(L) = FD(L) - (dv(L) - dx(L)*drdv)*dr3i
    4 CONTINUE
      PHI(I) = PHI(I) + BODY(J)*SQRT(dr2i)
*
*       Add contributions from each component to yield net correction.
      K = J1
    5 dr2 = 0.0
      drdv = 0.0
      DO 6 L = 1,3
          dx(L) = X(L,K) - X(L,I)
          dv(L) = XDOT(L,K) - XDOT(L,I)
          dr2 = dr2 + dx(L)**2
          drdv = drdv + dx(L)*dv(L)
    6 CONTINUE
*
      dr2i = 1.0/dr2
      dr3i = BODY(K)*dr2i*SQRT(dr2i)
      drdv = 3.0*drdv*dr2i
*
      DO 10 L = 1,3
          FIRR(L) = FIRR(L) + dx(L)*dr3i
          FD(L) = FD(L) + (dv(L) - dx(L)*drdv)*dr3i
   10 CONTINUE
      PHI(I) = PHI(I) - BODY(K)*SQRT(dr2i)
*
*       Treat the second component similarly.
      IF (K.EQ.J1) THEN
          K = J1 + 1
          GO TO 5
      END IF
*
   20 RETURN
*
      END
