      SUBROUTINE OUTPUT
*
*
*       Output and data save.
*       ---------------------
*
      INCLUDE 'common1.h'
      REAL*4  RDENS(3),A(8),XS(3,NMAX)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*
*       Predict X & XDOT for all particles (TIME > 0).
      IF (TIME.GT.0.0) THEN
          CALL XVPRED(1,N)
      END IF
*
*       Obtain the total kinetic & potential energy at current time.
      CALL ENERGY
*
*       Form approximate half-mass radius.
      RSCALE = 0.5*ZMASS**2/POT
*
*       Form virial theorem ratio (NB! VIR differs from POT if EPS2 > 0).
      Q = ZKIN/POT
*
*       Define the standard crossing time. (or code units)
      ETOT = ZKIN - POT
      TCR = ZMASS**2.5/(2.0*ABS(ETOT))**1.5
      IF(KZ(15).EQ.1) TCR=1
*
*       Update energies and form the relative error using MAX(ZKIN,POT).
      IF (TIME.EQ.0.0D0) THEN
          DE = 0.0
          BE(1) = ETOT
          BE(3) = ETOT
      ELSE
          BE(2) = BE(3)
          BE(3) = ETOT
          DE = (BE(3) - BE(2))/MAX(ZKIN,POT)
          ERRTOT = ERRTOT + DE
      END IF
*
*       Initialize c.m. integrals.
      DO 10 K = 1,3
          CMR(K) = 0.0
          CMRDOT(K) = 0.0
          RDENS(K) = 0.0
   10 CONTINUE
*
*       Obtain Z angular momentum & c.m. integrals.
      AZ = 0.0
      DO 30 I = 1,N
          AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
          DO 25 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
      RCM = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)/ZMASS
      VCM = SQRT(CMRDOT(1)**2 + CMRDOT(2)**2 + CMRDOT(3)**2)/ZMASS
*
*       Define tidal radius for isolated system (2*RTIDE used in ESCAPE).
      RTIDE = 10.0*RSCALE
*
*       Check print frequency indicator & optional model counter.
      NPRINT = NPRINT + 1
      IF (NPRINT.GT.NFIX.OR.TIME.LE.0.0D0) THEN
          NPRINT = 1
          IF (KZ(3).GT.0) MODEL = MODEL + 1
      END IF
*
      TC = TIME/TCR
      I6 = TSTAR*TIME
      WRITE (6,40)  TIME, Q, NSTEPI, DE, BE(3), TC
   40 FORMAT (//,' T =',F6.1,'  Q =',F5.2,'  STEPS =',I7,'  DE =',F10.6,
     &                                      '  E =',F10.6,'  TC =',F6.1)
      WRITE (6,50)  RSCALE, RCM, VCM, AZ, I6, NRUN
   50 FORMAT (/,' <R> =',F6.2,'  RCM =',F8.4,'  VCM =',F8.4,
     &                       '  AZ =',F10.5,'  T6 =',I4,'  NRUN =',I3,/)
*
*       Update next output time.
      TNEXT = TNEXT + DELTAT
*
*       Check optional output of single bodies & binaries.
      CALL BODIES
*
*       Perform automatic error control (restart with KZ(2) > 1).
      CALL CHECK(DE)
      IF (ABS(DE).GT.5.0*QE) GO TO 170
*
*       Check optional output of data bank (frequency NFIX).
      IF (KZ(3).EQ.0.OR.NPRINT.NE.1) GO TO 160
*
*       Write NK basic data on unit 3 (copy X(K,I) to *4 array if REAL*8).
      A(1) = TIME
      A(2) = RSCALE
      A(3) = RDENS(1)
      A(4) = RDENS(2)
      A(5) = RDENS(3)
      A(6) = BE(3)
      A(7) = RTIDE
      A(8) = 0.0
*
      DO 150 I = 1,N
          DO 145 K = 1,3
              XS(K,I) = X(K,I)
  145     CONTINUE
  150 CONTINUE
*
*       Split into WRITE (3) N & WRITE (3) (A(K), if disc instead of tape.
      IF (FIRST) THEN
          OPEN (UNIT=3,STATUS='NEW',FORM='UNFORMATTED',FILE='OUT3')
          FIRST = .FALSE.
      END IF
      NK = 8
      WRITE (3)  N, MODEL, NRUN, NK
      WRITE (3)  (A(K),K=1,NK), (BODY(J),J=1,N),
     &           ((XS(K,J),K=1,3),J=1,N), ((XDOT(K,J),K=1,3),J=1,N),
     &           (NAME(J),J=1,N)
*       Assume that the OS closes the file (no system crash or bug crash!).
*     CLOSE (UNIT=3)
*
*       See whether any escapers should be removed.
  160 IF (KZ(13).GT.0) THEN
          CALL ESCAPE
      END IF
*
*       Check adjustment of coordinates & velocities to c.m. condition.
      IF (KZ(14).GT.0) THEN
          CALL CMCORR
      END IF
*
*       Restore single precision velocities for coordinate predictions.
  170 DO 180 I = 1,N
          DO 175 K = 1,3
              XDOT(K,I) = X0DOT(K,I)
  175     CONTINUE
  180 CONTINUE
*
*       Obtain elapsed CPU time and update total since last output/restart.
      CALL CPUTIM(TCOMP)
      CPUTOT = CPUTOT + TCOMP - CPU0
      CPU0 = TCOMP
*
*       See whether COMMON should be saved after energy check.
      IF (KZ(2).GE.1) CALL MYDUMP(1,2)
*
*       Check termination time.
      IF (TIME.LT.TCRIT) RETURN
*       Save COMMON each time by CALL MYDUMP(1,1) before TIME < TCRIT test.
*
*       Terminate after optional COMMON save.
      WRITE (6,190)  TIME, CPUTOT, ERRTOT
  190 FORMAT (/,9X,'END RUN',3X,'TIME =',F8.2,'  CPUTOT =',F6.2,
     &                                                '  ERRTOT =',F9.5)
      IF (KZ(1).GT.0) CALL MYDUMP(1,1)
*
      STOP
*
      END
