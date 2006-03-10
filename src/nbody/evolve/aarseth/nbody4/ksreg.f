      SUBROUTINE KSREG
*
*
*       New KS regularization.
*       ----------------------
*
      INCLUDE 'common4.h'
      REAL*8  SAVE(21)
      EXTERNAL RENAME
*
*
*       Subtract dominant F & FDOT terms for routine fpolyi.F.
      CALL FFDOT2(ICOMP,JCOMP)
      CALL FFDOT2(JCOMP,ICOMP)
*
*       Save basic variables for components unless in correct location.
      DO 10 KCOMP = 1,2
*       Treat the first & second component in turn.
          IF (KCOMP.EQ.1) THEN
              I = ICOMP
          ELSE
              I = JCOMP
          END IF
          J = 2*NPAIRS + KCOMP
          IF (I.EQ.J) GO TO 10
*
          DO 2 K = 1,3
              SAVE(K) = X(K,I)
              SAVE(K+3) = X0DOT(K,I)
              SAVE(K+14) = F(K,I)
              SAVE(K+17) = FDOT(K,I)
    2     CONTINUE
          SAVE(7) = BODY(I)
          SAVE(8) = RADIUS(I)
          SAVE(9) = TEV(I)
          SAVE(10) = BODY0(I)
          SAVE(11) = EPOCH(I)
          SAVE(12) = TEV0(I)
          SAVE(13) = SPIN(I)
          SAVE(14) = ZLMSTY(I)
          SAVE(21) = PHI(I)
          NAMEI = NAME(I)
          KSI = KSTAR(I)
*
*       Exchange first & second single particle with ICOMP & JCOMP.
          DO 4 K = 1,3
              X(K,I) = X(K,J)
              X0(K,I) = X0(K,J)
              X0DOT(K,I) = X0DOT(K,J)
              XDOT(K,I) = XDOT(K,J)
              F(K,I) = F(K,J)
              FDOT(K,I) = FDOT(K,J)
              D2(K,I) = D2(K,J)
              D3(K,I) = D3(K,J)
              X(K,J) = SAVE(K)
              X0DOT(K,J) = SAVE(K+3)
              F(K,J) = SAVE(K+14)
              FDOT(K,J) = SAVE(K+17)
    4     CONTINUE
*
          BODY(I) = BODY(J)
          RADIUS(I) = RADIUS(J)
          ZLMSTY(I) = ZLMSTY(J)
          SPIN(I) = SPIN(J)
          TEV(I) = TEV(J)
          TEV0(I) = TEV0(J)
          BODY0(I) = BODY0(J)
          EPOCH(I) = EPOCH(J)
          NAME(I) = NAME(J)
          KSTAR(I) = KSTAR(J)
          STEP(I) = STEP(J)
          T0(I) = T0(J)
          PHI(I) = PHI(J)
          BODY(J) = SAVE(7)
          RADIUS(J) = SAVE(8)
          ZLMSTY(J) = SAVE(14)
          SPIN(J) = SAVE(13)
          TEV(J) = SAVE(9)
          TEV0(J) = SAVE(12)
          BODY0(J) = SAVE(10)
          EPOCH(J) = SAVE(11)
          PHI(J) = SAVE(21)
          NAME(J) = NAMEI
          KSTAR(J) = KSI
   10 CONTINUE
*
*       Increase pair index, total number & single particle index.
      NPAIRS = NPAIRS + 1
      NTOT = N + NPAIRS
      IFIRST = 2*NPAIRS + 1
*
*       Update all relevant COMMON list arrays.
      CALL RENAME
*
*       Set components in JLIST(1&2) for updating NLIST in FPOLY.
      JLIST(1) = ICOMP
      JLIST(2) = JCOMP
*
*       Initialize the regularized solution.
      CALL KSINIT
*
*       Check updating of global index & neighbour list for chain c.m.
      IF (NCH.GT.0) THEN
          CALL CHFIND
      END IF
*
      RETURN
*
      END
