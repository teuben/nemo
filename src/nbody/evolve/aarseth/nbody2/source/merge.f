      SUBROUTINE MERGE
*
*
*       Merging of colliding bodies.
*       ----------------------------
*
      INCLUDE 'common2.h'
*
*
*       Predict current coordinates & velocities of second body.
      CALL XVPRED(JCOMP,JCOMP)
*
*       Ensure ICOMP < JCOMP for merger procedure.
      I = JCOMP
      JCOMP = MAX(ICOMP,JCOMP)
      ICOMP = MIN(ICOMP,I) 
*
*       Obtain the potential energy with respect to both components.
      POT1 = 0.0
      I = ICOMP
      NNB1 = LIST(1,ICOMP) + 1
    4 DO 8 L = 2,NNB1
          J = LIST(L,ICOMP)
          IF (J.EQ.JCOMP) GO TO 8
          RIJ2 = EPS2
          DO 6 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
    6     CONTINUE
          POT1 = POT1 + BODY(I)*BODY(J)/SQRT(RIJ2)
    8 CONTINUE
*
      IF (I.EQ.ICOMP) THEN
          I = JCOMP
          GO TO 4
      END IF
*
*       Form the two-body binding energy.
      RIJ2 = 0.0
      VIJ2 = 0.0
      DO 10 K = 1,3
          RIJ2 = RIJ2 + (X(K,ICOMP) - X(K,JCOMP))**2
          VIJ2 = VIJ2 + (XDOT(K,ICOMP) - XDOT(K,JCOMP))**2
   10 CONTINUE
*
      BCM = BODY(ICOMP) + BODY(JCOMP)
      A1 = 2.0/SQRT(RIJ2 + EPS2) - VIJ2/BCM
      EPAIR = -0.5*BODY(ICOMP)*BODY(JCOMP)*A1
*
*       Set new variables for the merged body in ICOMP.
      DO 20 K = 1,3
          X(K,ICOMP) = (BODY(ICOMP)*X(K,ICOMP) +
     &                                       BODY(JCOMP)*X(K,JCOMP))/BCM
          XDOT(K,ICOMP) = (BODY(ICOMP)*XDOT(K,ICOMP) +
     &                                    BODY(JCOMP)*XDOT(K,JCOMP))/BCM
   20 CONTINUE
*
*       Set new mass and check maximum mass.
      BODYI = BODY(ICOMP)
      BODY(ICOMP) = BCM
      BODY1 = MAX(BCM,BODY1)
*
*       Obtain the potential energy with respect to new c.m.
      POT2 = 0.0
      I = ICOMP
      DO 40 L = 2,NNB1
          J = LIST(L,I)
          IF (J.EQ.JCOMP) GO TO 40 
          RIL2 = EPS2
          DO 30 K = 1,3
              RIL2 = RIL2 + (X(K,I) - X(K,J))**2
   30     CONTINUE
          POT2 = POT2 + BODY(I)*BODY(J)/SQRT(RIL2)
*       Reduce the neighbour steps because of discontinuity.
          DT = TIME - T0(J)
          STEP(J) = MAX(0.5*STEP(J),DT)
   40 CONTINUE
*
*       Correct the total energy for internal energy and tidal effect.
      BE(3) = BE(3) - EPAIR - (POT2 - POT1)
*
      IF (KZ(12).GT.1) THEN
          WRITE (6,50)  NAME(ICOMP), NAME(JCOMP), BODYI, BODY(JCOMP),
     &                  SQRT(RIJ2), EPAIR, BE(3), (POT2 - POT1)
   50     FORMAT (5X,' MERGER    NAME =',2I6,'  MASS =',2F8.4,'  RIJ =',
     &                         1P,E9.1,'  EPAIR =',E9.1,'  E =',0P,F9.5,
     &                                                   '  DP =',F10.5)
      END IF
*
*       Reduce particle number and update all COMMON arrays.
      N = N - 1
      CALL REMOVE(JCOMP,2)
*
*       Obtain new neighbour list if merged body has < 2 neighbours.
      IF (LIST(1,ICOMP).LT.2) THEN
          RSI = 2.0*RS(ICOMP)
          CALL NBLIST(ICOMP,RSI)
      END IF
*
*       Predict current coordinates & velocities of the neighbours to F2DOT.
      DO 60 L = 2,NNB1
          J = LIST(L,ICOMP)
          DT = TIME - T0(J)
          A3 = (T0(J) - T1(J)) + (T0(J) - T2(J))
          A4 = (T0(J) - T0R(J)) + (T0(J) - T1R(J)) + (T0(J) - T2R(J))
*
          DO 55 K = 1,3
              F2DOTK = D3R(K,J)*A4 + D2R(K,J) + D3(K,J)*A3 + D2(K,J)
              X(K,J) = (((ONE12*F2DOTK*DT + FDOT(K,J))*DT + F(K,J))*DT +
     &                                          X0DOT(K,J))*DT + X0(K,J)
              XDOT(K,J) = ((ONE3*F2DOTK*DT + 3.0*FDOT(K,J))*DT +
     &                                       2.0*F(K,J))*DT + X0DOT(K,J)
   55     CONTINUE
   60 CONTINUE
*
*       Obtain F & FDOT followed by F2DOT & F3DOT and set time-steps.
      CALL FPOLY1(ICOMP,ICOMP)
      CALL FPOLY2(ICOMP,ICOMP)
      CALL STEPS(ICOMP,ICOMP)
*
*       Restore the current neighbour velocities for the predictor.
      DO 70 L = 2,NNB1
          J = LIST(L,ICOMP)
          DO 65 K = 1,3
              XDOT(K,J) = X0DOT(K,J)
   65     CONTINUE
   70 CONTINUE
*
*       Set phase indicator = -1 for new NLIST in routine INTGRT.
      IPHASE = -1
      NMERG = NMERG + 1
*
      RETURN
*
      END
