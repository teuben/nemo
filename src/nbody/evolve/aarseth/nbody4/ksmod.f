      SUBROUTINE KSMOD(IPAIR,KMOD)
*
*
*       Modified KS motion.
*       -------------------
*
      INCLUDE 'common4.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
*
*
*       Set current modification level.
      IMOD = KSLOW(IPAIR)
*
*       Check transition type (KMOD <= 1 denotes standard restart).
      IF (KMOD.GT.1) THEN
*       Determine provisional index for KS slow-down.
          DO 5 K = 2,10
              ISBIN = K - 1
              IF (ISLOW(K).GT.KMOD) GO TO 10
    5     CONTINUE
*       Restrict increase to two levels.
   10     ISBIN = MIN(ISBIN,IMOD+2)
*       See whether standard solution is called for.
          IF (ISBIN.EQ.1) GO TO 30
      ELSE
*       Include termination (IMOD > 1 & KMOD <= 1).
          ISBIN = 1
          GO TO 30
      END IF
*
*       Estimate time interval to reach largest permitted perturbation.
      GMX = RANGE*GMIN
      CALL TPERT(IPAIR,GMX,DT)
*
*       Evaluate the unmodified Kepler period.
      ICM = N + IPAIR
      SEMI = -0.5*BODY(ICM)/H(IPAIR)
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
*
*       Reduce level if modification factor is too large.
      DO 20 K = 2,10
          IF (TK*FLOAT(ISLOW(ISBIN)).LT.DT.OR.ISBIN.EQ.1) GO TO 30
          ISBIN = ISBIN - 1
   20 CONTINUE
*
*       Exit if level is unchanged.
   30 IF (ISBIN.EQ.IMOD) GO TO 100
*
*       Ensure the Kepler period is defined for transition to standard case.
      IF (ISBIN.EQ.1) THEN
          ICM = N + IPAIR
          SEMI = -0.5*BODY(ICM)/H(IPAIR)
          TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
      END IF
*
*       Define auxiliary quantities.
      ZETA = 1.0 - R(IPAIR)/SEMI
      PSI = TDOT2(IPAIR)/SQRT(BODY(ICM))
*
*       Determine the eccentric anomaly with respect to pericentre (-PI,PI).
      THETA = ATAN2(PSI/SQRT(SEMI),ZETA)
*
*       Obtain apocentre time interval from Kepler's equation and the period.
      DT = SEMI*SQRT(SEMI/BODY(ICM))*(THETA - PSI/SQRT(SEMI))
      DT = 0.5D0*TK + DT
      DT = -DT
*
*       Evaluate regularized apocentre time (Stiefel & Scheifele, p. 85).
*     DTU = 2.0D0*(H(IPAIR)*DT - 0.5D0*TDOT2(IPAIR))/BODY(ICM)
*
*       Determine the regularized step by Newton-Raphson iteration (DT < 0).
      I1 = 2*IPAIR - 1
      DTU = DT/R(IPAIR)
      ITER = 0
*       Note: explicit relation agrees with iterated value (bug fix 9/99).
   40 Y0 = DT - ((ONE6*TDOT3(IPAIR)*DTU +
     &                             0.5*TDOT2(IPAIR))*DTU + R(IPAIR))*DTU
      YPR = -((0.5*TDOT3(IPAIR)*DTU + TDOT2(IPAIR))*DTU + R(IPAIR))
      DTU = DTU - Y0/YPR
      DT1 = ((ONE6*TDOT3(IPAIR)*DTU + 0.5*TDOT2(IPAIR))*DTU +
     &                                                     R(IPAIR))*DTU
      ITER = ITER + 1
      IF (ABS(DT - DT1).GT.1.0D-10*STEP(I1).AND.ITER.LT.5) GO TO 40
*
*       Reset reference energy and generate new Stumpff coefficients.
      H0(IPAIR) = H(IPAIR)
      Z = -0.5*H(IPAIR)*DTU**2
      CALL STUMPF(IPAIR,Z)
*
*       Integrate small step back to apocentre using temporary indicator # 0.
      DTAU(IPAIR) = DTU
      TIME = T0(I1) + DT
      IPHASE = -1
      CALL KSINT(I1)
      IPHASE = 0
*
*       Predict current coordinates & velocities for perturbers, ICM & IPAIR.
      NNB2 = LIST(1,I1) + 2
      LIST(NNB2,I1) = ICM
      DO 60 L = 2,NNB2
          J = LIST(L,I1)
          CALL XVPRED(J,0)
   60 CONTINUE
*
*       Set new KS level and increase restart counter (ISBIN > 1).
      KSLOW(IPAIR) = ISBIN
      IF (ISBIN.GT.1) NKSMOD = NKSMOD + 1
*
*       Check perturbers for small GAMMA (note KSLIST skip in KSINT).
      IF (GAMMA(IPAIR).LT.0.01*GMIN) THEN
          CALL KSLIST(IPAIR)
*       Warning: ISBIN up to 9 and GAMMA < 1D-12 seen without this check.
      END IF
*
*       Perform perturbed restart of KS motion if needed.
      IF (LIST(1,I1).GT.0) THEN
          CALL KSPOLY(IPAIR,ISBIN)
      END IF
*
*       Ensure that apocentre criterion will fail after next perturbed step.
      TDOT2(IPAIR) = 0.0D0
*
*       Set indicator = -1 to skip perturber selection in routine KSINT.
      KMOD = -1
*
  100 RETURN
*
      END
