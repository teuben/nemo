*
*             T R I P L E
*             ***********
*
*
*       Three-body regularization program.
*       ----------------------------------
*
*       Method of Aarseth & Zare, Celestial Mechanics 10, 185.
*       ......................................................
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*
      PROGRAM TRIPLE
*
      IMPLICIT  REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  Q(8),P(8),R,R1,R2,ENERGY,M(3),X(3,3),XDOT(3,3),
     &               RCOLL,ERROR,C11,C12,C19,C20,C24,C25,NSTEPS,NAME(3)
      COMMON/CLOSE/  RIJ(3,3),ICALL
      COMMON/BSSAVE/  EP(4),STEPC,TFAC,ITFAC,JC,NHALF2
      COMMON/ANGLES/  TWOPI,ALPHA,DEGREE
      COMMON/CALLS/   NFN
      COMMON/TBIN0/   RGRAV,EBIN0,ECC0,DB,I10,I20
      REAL*8  Y(17),RK(3),M0(3),X0(3,3),V0(3,3)
      REAL*8  XX(3,3),XMAX
*
*
*       COMMON variables
*       ****************
*
*       ------------------------------------------------------------------
*       C11     Inverse mass factor for DERQP (also C12,C19,C20,C24,C25).
*       ENERGY  Twice the initial total energy.
*       ERROR   Relative energy error at the end (ignore with T' = 1/U).
*       ICALL   Indicator for pericentre check (first function call only).
*       M       Particle mass.
*       NAME    Particle identity (initialized to 1,2,3).
*       NFN     Number of function calls.
*       NSTEPS  Number of DIFSY calls.
*       P       Regularized momenta.
*       Q       Regularized coordinates.
*       R       Distance between M(1) and M(2) (not regularized).
*       R1      Distance between M(1) and M(3).
*       R2      Distance between M(2) and M(3).
*       RCOLL   Minimum two-body separation (osculating pericentre).
*       RIJ     Minimum pairwise separations (osculating pericentre).
*       X       Particle coordinates (X(3,I) is Z-component).
*       XDOT    Velocity components (XDOT(3,I) is Z-component).
*       ------------------------------------------------------------------
*
*       Input Parameters
*       ------------------------------------------------------------------
*       TOL     Tolerance for Bulirsch-Stoer integrator.
*       TCRIT0  Maximum time (terminates on escape or TCRIT).
*       ECRIT   Energy error criterion (time reversal only).
*       IREV    Time reversal option (IREV = 0: standard).
*       IMOVE   Movie option (> 0).
*       XMAX    Maximum range (X & Y; movie option).
*       NWAIT   Delay loop (CALL WAIT; depends on CPU speed.
*       M,R,V   Mass, coordinates & velocities (routine DATA).
*       Template:
*       1.0D-12 100.0 1.0D-03 0 1 5.0
*       ------------------------------------------------------------------
*
*       Read and print main parameters (IREV > 0 for time reversal).
      READ (5,*)  TOL, TCRIT0, ECRIT, IREV, IMOVE, XMAX, NWAIT
      WRITE (6,999)  TOL, TCRIT0, ECRIT, IREV, XMAX
  999 FORMAT (//,5X,'TOL =',1P,E9.2,'  TCRIT =',E9.2,'  ECRIT =',E8.1,
     &              '  IREV =',0P,I3,'  XMAX =',F4.1,/)
*
      TWOPI = 8.0D0*ATAN(1.0D0)
      DEGREE = 360.0D0/TWOPI
      DELTA = 0.0
      ALPHA = 0.0
      PC = 0.0
      SEMI = 0.0
      E = 0.0
      IESC = 0
*
*       Specify the number of first-order equations for the integrator.
      N = 17
      IEXP = 0
      IF (IMOVE.GT.0) THEN
          CALL MOVIE0(XMAX)
      END IF
*
*       Read the basic configuration.
    1 CALL DATA
*
*       Save the tolerance, repeat counter & initial condition.
      EPS = TOL
      IERR = 0
      DO 3 I = 1,3
          M0(I) = M(I)
          DO 2 K = 1,3
              X0(K,I) = X(K,I)
              V0(K,I) = XDOT(K,I)
              XX(K,I) = X(K,I)
    2     CONTINUE
    3 CONTINUE
*
    5 PERIM = 100.0
      ZMI1 = 1000.0
      ZMI2 = 1000.0
      ZMI3 = 1000.0
      ZMIN = 1000.0
*
*       Set control indicators.
      IRUN = 0
      ITER = 0
*       Specify step increase factor and initialization index for DIFSY.
      STEPC = 1.1D0
      JC = -1
*
*       Initialize diagnostic variables & counters.
      R12MIN = 100.0
      RMIN = 100.0
      RCOLL = 100.0
      DO 10 J = 1,3
          DO 8 K = 1,3
              RIJ(J,K) = 100.0
    8     CONTINUE
   10 CONTINUE
      NSTEPS = 0
      NFN = 0
      NREG = 0
      ICALL = 0
*
*       Initialize local time, regularized time and termination time..
      TIME = 0.0D0
      TAU = 0.0D0
      Y(17) = 0.0
      TCRIT = TCRIT0
*
*       Obtain initial energy and transform to regularized variables.
      IF (IRUN.EQ.0) CALL TRANSF(1)
*
*       Define gravitational radius and pericentre check distance.
      RGRAV = (M(1)*M(2) + M(1)*M(3) + M(2)*M(3))/(0.5D0*ABS(ENERGY))
      RSTAR = 0.5*RGRAV
*
*       Form the two smallest distances (assume sensible reference body).
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
*
*       Set termination distance to maximum of separation and 5*RGRAV.
      RMAX0 = MAX(R1,R2,5.0*RGRAV)
*
*       Specify the crossing time (also meaningful if ENERGY > 0).
      TCR = (M(1) + M(2) + M(3))**2.5/ABS(ENERGY)**1.5
*
*       Set step for time transformation DT/DTAU = 1.0/U (ignore M1*M2/R).
*     TPR = R1*R2/(M(1)*M(3)*R2 + M(2)*M(3)*R1)
      TPR = R1*R2/SQRT(R1 + R2)
      DTAU = TCR*EPS**0.1/TPR
*
*       Try alternative expression for initial step (Seppo Mikkola).
      R0 = MAX(R1,R2)
      DTAU = 0.1*(EPS/1.0E-12)**0.1*R0*SQRT(R0/(M(1) + M(2) + M(3)))/TPR
*
      IF (IERR.EQ.0) THEN
          WRITE (6,15)  RGRAV, R1, R2, 0.5*ENERGY
   15     FORMAT (5X,'RGRAV =',F7.3,'  R1 =',F7.3,'  R2 =',F7.3,
     &               '  ENERGY =',F12.6,/)
          CALL BINARY(IESC)
      END IF

*       Initialize time constant & input array for the integrator.
      CONST = 0.0D0
      DO 20 K = 1,8
          CONST = CONST + Q(K)*P(K)
          Y(K) = Q(K)
          Y(K+8) = P(K)
   20 CONTINUE
*
*       Produce initial output.
      IF (IREV.EQ.0) CALL TRANSF(2)
*
*       Advance the equations of motion by Bulirsch-Stoer integrator.
   30 IF (IMOVE.GT.0) DTAU = 0.1*DTAU
      CALL DIFSY1(N,EPS,DTAU,TAU,Y)
*
      IF (DTAU.EQ.0.0) THEN
          WRITE(6,*)  ' STEPSIZE = 0!', char(7)
          STOP
      END IF

*       Copy regularized coordinates & momenta and obtain physical time.
      SUMQP = 0.0D0
      DO 40 K = 1,8
          Q(K) = Y(K)
          P(K) = Y(K+8)
          SUMQP = SUMQP + Q(K)*P(K)
*       Note that the momentum includes factor of 2 in AZ formulation.
   40 CONTINUE
*
*       Set explicit time (Baumgarte & Stiefel, 1974; Aarseth, 1976).
*     TEXP = (0.5D0*(SUMQP - CONST) - TAU)/ENERGY
      TIME = Y(17)
*
*       Update relative distances (NB! not quite latest value).
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
*
*       Check minimum two-body separations and increase step counter.
      RMIN = MIN(RMIN,R)
      RM = MIN(R1,R2)
      R12MIN = MIN(R12MIN,RM)
      RMAX = MAX(R1,R2,R)
      NSTEPS = NSTEPS + 1
*
*       Check minimum two-body separations.
      RK(1) = R1
      RK(2) = R2
      RK(3) = R
*       Consider pairs 1-2, 1-3 & 2-3 with identified names.
      DO 44 K = 1,3
          DO 42 L = K+1,3
              I = NAME(K)
              J = NAME(L)
*       Use cyclic loop index (3,1,2) for distances R, R1 & R2.
              KK = K - 1
              IF (KK.EQ.0) KK = 3
              RIJ(I,J) = MIN(RIJ(I,J),RK(KK))
              RIJ(J,I) = MIN(RIJ(J,I),RK(KK))
   42     CONTINUE
   44 CONTINUE
*
*       Switch on search indicator inside RSTAR (reset in DERQP).
      IF (RM.LT.RSTAR) THEN
          ICALL = 1
      END IF
*
*       Determine pericentre & impact parameter during forward motion.
      IF (IRUN.EQ.0.AND.R1 + R2 + R.LT.PERIM) THEN
          PERIM = R1 + R2 + R
          IF (PERIM.LT.3.0*RGRAV) THEN
              CALL IMPACT(PC)
          END IF
      END IF
*
*       Record minimum moment of inertia during forward integration.
      IF (IRUN.EQ.0) THEN
          ZMI = R1**2 + R2**2 + R**2
          ZMI1 = ZMI2
          ZMI2 = ZMI3
          ZMI3 = ZMI
          IF (ZMI3.GT.ZMI2.AND.ZMI2.LT.ZMI1) THEN
              ZMIN = MIN(ZMI2,ZMIN)
          END IF
      END IF
*
*       Include optional movie 
      IF (IMOVE.GT.0) THEN
          DO 45 L = 1,3
              IF (NAME(L).EQ.1) I1 = L
              IF (NAME(L).EQ.2) I2 = L
              IF (NAME(L).EQ.3) I3 = L
   45     CONTINUE
          CALL TRANSF(3)
          DO 48 I = 1,3
              DO 46 K = 1,3
                  XX(K,I) = X(K,I)
   46         CONTINUE
   48     CONTINUE
          CALL MOVIE(I1,I2,I3,XX,TIME)
          ZZ = 0.0
          DO 49 KK = 1,NWAIT
              CALL WAIT(KK,ZZ)
   49     CONTINUE
          IF (ZZ.LT.-1.0) STOP
      END IF
*
*       See whether switching of reference body is desirable.
      IF (R.GT.RM) GO TO 70
      IMIN = 1
*       Use a simple distance test to determine new reference body IMIN.
      IF (R2.LT.1.00001*R1) IMIN = 2
*
*       Transform to physical variables and rename the exchanged particles.
      CALL TRANSF(3)
*
      DO 50 K = 1,3
          TEMP1 = X(K,3)
          TEMP2 = XDOT(K,3)
          X(K,3) = X(K,IMIN)
          XDOT(K,3) = XDOT(K,IMIN)
          X(K,IMIN) = TEMP1
          XDOT(K,IMIN) = TEMP2
   50 CONTINUE
*
      TEMP1 = M(3)
      M(3) = M(IMIN)
      M(IMIN) = TEMP1
      NAME3 = NAME(3)
      NAME(3) = NAME(IMIN)
      NAME(IMIN) = NAME3
*
*       Transform back to regularized variables and initialize input array.
      CALL TRANSF(4)
      DO 60 K = 1,8
          Y(K) = Q(K)
          Y(K+8) = P(K)
   60 CONTINUE
*
*       Increase regularization counter at the end of switching procedure.
      NREG = NREG + 1
*
*       Update relative distances (NB! New R becomes largest old R1 or R2).
      R = MAX(R1,R2)
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
*
*       Check termination criteria (TIME > TCRIT or RMAX > RMAX0).
   70 IF (IRUN.GT.0.AND.(TIME.GT.TCRIT.OR.ITER.GT.0)) THEN
          IF (ABS(TIME - TCRIT).LT.1.0E-12) GO TO 120
*       Begin iteration for termination at TCRIT (NB! use correct TPR).
          DTAU = (TCRIT - TIME)/(R1*R2)
          DTAU = DTAU*SQRT(R1 + R2)
*       Note extra scaling of STEP to compensate for movie reduction.
          IF (IMOVE.GT.0) DTAU = 10.0*DTAU
*       ------------------------------------------------------------
*       Include alternative time transformation (N = 16; see DERQP).
*         U = M(1)*M(3)*R2 + M(2)*M(3)*R1 + M(1)*M(2)*R1*R2/R
*         DTAU = DTAU*U
*       ------------------------------------------------------------
          ITER = ITER + 1
          IF (ITER.LE.9) GO TO 30
          GO TO 90
      END IF
      IF ((RMAX.LT.RMAX0.AND.TIME.LT.TCRIT).OR.IRUN.GT.0) GO TO 30
*
*       Obtain final output after transforming to physical variables.
      CALL TRANSF(2)
*
*       Perform binary analysis.
      CALL BINARY(IESC)
*
*       Print diagnostics during forward run.
   90 IF (IRUN.EQ.0) THEN
          RCOLL = MIN(RCOLL,R12MIN)
          WRITE (6,95)  TIME, RMIN, R12MIN, RCOLL, RMAX, ERROR, NSTEPS,
     &                  NREG
   95     FORMAT  (/,'  TIME =',F7.3,'  MIN(R) =',1P,E8.1,
     &               '  MIN(R1,R2) =',E8.1,'  RCOLL =',E8.1,
     &               '  RMAX =',E8.1,'  DE/E =',E9.1,'  NSTEPS =',I5,
     &               '  NREG =',I3)
          WRITE (6,100)  RIJ(1,2), RIJ(1,3), RIJ(2,3), RCOLL
  100     FORMAT (' RIJ:  1-2 1-3 2-3 RCOLL ',1P,4E10.2)
      END IF
*
*       Check time reversal indicator.
      IF (IREV.GT.0) THEN
          IF (IRUN.EQ.0) THEN
              TCRIT = TIME
              TAU = 0.0D0
              Y(17) = 0.0D0
              NSTEPS = 0
              NREG = 0
              NFN0 = NFN
              NFN = 0
*       Initialize time constant & input array for the integrator.
              CONST = 0.0D0
              DO 110 K = 1,8
                  P(K) = -P(K)
                  Y(K) = Q(K)
                  Y(K+8) = P(K)
                  CONST = CONST + Q(K)*P(K)
  110         CONTINUE
              IRUN = IRUN + 1
          END IF
          GO TO 30
      END IF
*
*       Continue with the next experiment.
      IEXP = IEXP + 1
      IF (IMOVE.GT.0) CALL MOVIE1
      GO TO 1
*
*       Transform to physical variables and evaluate time reversal errors.
  120 CALL TRANSF(3)
      ERRX = 0.0
      ERRV = 0.0
      DO 130 I = 1,3
          L = NAME(I)
          DO 125 K = 1,3
              ERRX = ERRX + (X(K,I) - X0(K,L))**2
              ERRV = ERRV + (XDOT(K,I) + V0(K,L))**2
  125     CONTINUE
  130 CONTINUE
      ERRX = SQRT(ERRX/3.0)
      ERRV = SQRT(ERRV/3.0)
      RCOLL = MIN(RCOLL,R12MIN)
*       Obtain angle of escape in degrees.
      BETA = ATAN2(XDOT(2,IESC),XDOT(1,IESC))
      IF (BETA.LT.0.0) BETA = BETA + TWOPI
      BETA = DEGREE*BETA
*
      WRITE (10,135)  NSTEPS, NFN0, NFN, ERROR, ERRX, ERRV, EPS, TIME,
     &                SEMI, E
  135 FORMAT (' ',I4,2I7,1P,E9.1,3E9.1,0P,F7.2,2F7.4)
*
*       Print diagnostic (successful run on unit #7, unsuccessful on #8).
      IF (MAX(ERRX,ERRV).LT.ECRIT) THEN
          WRITE (7,140)  IEXP, TIME, DELTA, ALPHA, IESC, PC, BETA,
     &                   SEMI, E, RCOLL, PERIM, ZMIN, ERRX, ERRV, EPS
  140     FORMAT (I5,F7.1,F8.4,F8.2,I3,F7.3,F7.1,2F7.3,1P,6E10.2)
      ELSE
          WRITE (8,140)  IEXP, TIME, DELTA, ALPHA, IESC, PC, BETA,
     &                   SEMI, E, RCOLL, PERIM, ZMIN, ERRX, ERRV, EPS
*
*       Reduce tolerance by 100 and try again (until limit of 1.0E-15).
          IERR = IERR + 1
          IF (IERR.LT.7.AND.EPS.GT.0.99E-13) THEN
              EPS = 0.01*EPS
              DO 150 I = 1,3
                  M(I) = M0(I)
                  NAME(I) = I
                  DO 145 K = 1,3
                      X(K,I) = X0(K,I)
                      XDOT(K,I) = V0(K,I)
  145             CONTINUE
  150         CONTINUE
              GO TO 5
          END IF
      END IF
*
*       Continue with the next experiment.
      GO TO 1
*
      END
