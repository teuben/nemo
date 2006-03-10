      SUBROUTINE SEARCH(I,IKS)
*
*
*       Close encounter search.
*       -----------------------
*
      INCLUDE 'common4.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Increase counter for regularization attempts.
      NKSTRY = NKSTRY + 1
*
*       Find dominant neighbour by examining all STEP(J) =< 4*STEP(I).
      DTCR = MIN(4.0*STEP(I),SMIN)
      ITRY = 0
    1 FMAX = 0.0
      NCLOSE = 0
      JCOMP = 0
*
      NNB1 = LSHORT(1) + 1
      DO 6 L = 2,NNB1
          J = LSHORT(L)
          IF (STEP(J).GT.DTCR) GO TO 6
          A1 = X(1,J) - X(1,I)
          A2 = X(2,J) - X(2,I)
          A3 = X(3,J) - X(3,I)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*       Note index of every particle with small step inside 2*RMIN.
          IF (RIJ2.LT.RMIN22) THEN
              IF (J.EQ.I) GO TO 6
              NCLOSE = NCLOSE + 1
              JLIST(NCLOSE) = J
              FIJ = (BODY(I) + BODY(J))/RIJ2
              IF (FIJ.GT.FMAX) THEN
                  FMAX = FIJ
*       Save square distance and global index of dominant body.
                  RJMIN2 = RIJ2
                  JCOMP = J
              END IF
          END IF
    6 CONTINUE
*
*       See whether dominant component is a single particle inside RMIN.
      IF (JCOMP.LT.IFIRST.OR.JCOMP.GT.N) GO TO 10
      IF (RJMIN2.GT.RMIN2) GO TO 10
*
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
*       Only select approaching particles (include nearly circular case).
      RIJMIN = SQRT(RJMIN2)
      IF (RDOT.GT.0.02*SQRT((BODY(I) + BODY(JCOMP))*RIJMIN)) GO TO 10
*
*       Evaluate vectorial perturbation due to close bodies.
      IF (NCLOSE.GT.1) THEN
          CALL FPERT(I,JCOMP,NCLOSE,PERT)
      ELSE
          PERT = 0.0
      END IF
*
*       Check for second search (case of small STEP(I) after termination).
      IF (NCLOSE.EQ.1.AND.LSHORT(1).GT.2) THEN
          ITRY = ITRY + 1
          DTCR = SMIN
          IF (ITRY.LE.1) GO TO 1
      END IF
*
*       Accept #I & JCOMP if the relative motion is dominant (GI < 0.25).
      GI = PERT*RJMIN2/(BODY(I) + BODY(JCOMP))
      IF (GI.GT.0.25) THEN
*         IF (KZ(4).GT.0.AND.TIME-TLASTT.GT.4.44*TCR/FLOAT(N))
*    &                                             CALL EVOLVE(JCOMP,0)
          GO TO 10
      END IF
*
*       Exclude any c.m. body of compact subsystem (I < N: TRIPLE or QUAD).
      DO 8 ISUB = 1,NSUB
          NAMEI = NAMES(1,ISUB)
          IF (NAMEI.EQ.NAME(I).OR.NAMEI.EQ.NAME(JCOMP)) GO TO 10
    8 CONTINUE
*
*       Also check possible c.m. body of chain regularization (NAME = 0).
      IF (NCH.GT.0) THEN
          IF (NAME(I).EQ.0.OR.NAME(JCOMP).EQ.0) GO TO 10
      END IF
*
*       Save index and increase indicator to denote new regularization.
      ICOMP = I
      IKS = IKS + 1
*
   10 RETURN
*
      END
