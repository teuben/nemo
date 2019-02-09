      SUBROUTINE OUTPUT
*
*
*       Output and data save.
*       ---------------------
*
      INCLUDE 'common2.h'
      REAL*4  A(12),XS(3,NMAX)
      LOGICAL  FIRST
      DATA  FIRST /.TRUE./
      SAVE  FIRST
*
*
*       Predict X & XDOT for all particles (TIME > 0).
      IF (TIME.GT.0.0) THEN
          CALL XVPRED(1,N)
      END IF
*
*       Obtain the total energy at current time.
      CALL ENERGY
*
*       Set provisional half-mass radius (used for density contrast).
      IF (KZ(10).EQ.0) THEN
          RSCALE = 0.5*ZMASS**2/POT
      ELSE
          RSCALE = 2.0*RSMIN*(FLOAT(N)/100.0)**0.333
      END IF
*
*       Initialize c.m. integrals.
      DO 10 K = 1,3
          CMR(K) = 0.0
          CMRDOT(K) = 0.0
   10 CONTINUE
*
*          Perform time-step & neighbour statistics.
      DTI = 0.0
      DTRI = 0.0
      CNNB = 0.0
      CMAX = 0.0
      NNB = 0
      DO 20 I = 1,N
          DTI = DTI + 1.0/STEP(I)
          DTRI = DTRI + 1.0/STEPR(I)
          CNNB = CNNB + LIST(1,I)/STEP(I)
          NNB = NNB + LIST(1,I)
          RHO = LIST(1,I)/RS(I)**3
*       Find provisional density centre from maximum neighbour density.
          IF (RHO.GT.CMAX) THEN
              CMAX = RHO
              IMAX = I
          END IF
   20 CONTINUE
*
*          Estimate relative cost & effective neighbour number of AC scheme.
      COST = CNNB/(FLOAT(N)*DTRI)
      CNNB = CNNB/DTI
*
*       Set average neighbour number and update option 14 if non-zero.
      NNB = FLOAT(NNB)/FLOAT(N)
      IF (KZ(14).GT.0) KZ(14) = NNB
*
*       Find provisional density centre (adopt zero if density contrast < 5).
      A2 = 2.0*CMAX*RSCALE**3/FLOAT(N)
      DO 25 K = 1,3
          RDENS(K) = X(K,IMAX)
          IF (A2.LT.5.0.OR.KZ(10).GT.0) RDENS(K) = 0.0
   25 CONTINUE
*
*       Find density centre & core radius (Casertano & Hut, AP.J. 298, 80).
      IF (KZ(8).GT.0.AND.N.GE.20) THEN
          CALL CORE
      ELSE
          NC = N/2
          RC = RSCALE
          RC2 = RC**2
          RHOM = 0.0
          RHOD = 0.0
          VC = 0.0
      END IF
*
*       Check optional sorting of Lagrangian radii & half-mass radius.
      IF (KZ(7).GT.0) THEN
          CALL LAGR(RDENS)
      END IF
*
*       Scale CMAX, average & maximum core mass density by the mean value.
      CMAX = 2.0*CMAX*RSCALE**3/FLOAT(N)
      RHOD = 8.3773*RHOD*RSCALE**3/ZMASS
      RHOM = 8.3773*RHOM*RSCALE**3/ZMASS
*
*       Obtain Z angular momentum & c.m. integrals.
      AZ = 0.0
      DO 30 I = 1,N
          AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
          DO 28 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   28     CONTINUE
   30 CONTINUE
*
      RD = SQRT(RDENS(1)**2 + RDENS(2)**2 + RDENS(3)**2)
      RCM = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)/ZMASS
      VCM = SQRT(CMRDOT(1)**2 + CMRDOT(2)**2 + CMRDOT(3)**2)/ZMASS
*
*       Form virial theorem ratio (NB! VIR differs from POT if EPS2 > 0).
      Q = ZKIN/VIR
*
*       Define standard crossing time.
      ETOT = ZKIN - POT - ETIDE
      TCR = ZMASS**2.5/(2.0*ABS(ETOT))**1.5
      IF(KZ(19).EQ.1) TCR=1
*
*       Update energies and form the relative error using MAX(ZKIN,POT).
      IF (TIME.LE.0.0D0) THEN
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
*       Update density contrast factor for neighbour sphere modification.
      ALPHA = FLOAT(NNBMAX)*SQRT(0.08*RSCALE**3/FLOAT(N))
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
      WRITE (6,40)  TIME, Q, N, NNB, NSTEPI, NSTEPR, DE, BE(3), TC, I6
   40 FORMAT (//,' T =',F6.1,'  Q =',F5.2,'  N =',I5,'  <NB> =',I3,
     &           '  NSTEPS =',I11,I10,'  DE =',F10.6,'  E =', F9.5,
     &           '  TC =',F6.1,'  T6 =',I5)
*
      CALL CPUTIM(TCOMP)
      WRITE (6,50)
   50 FORMAT (/,'    <R>  RTIDE  RDENS   RC   NC    MC   RHOD   RHOM',
     &          '  CMAX  <Cn>  Ir/R  NRUN  MODEL    CPU',
     &          '   RCM    VCM       AZ     RSMIN   BODY1')
      WRITE (6,55)  RSCALE, RTIDE, RD, RC, NC, ZMC, RHOD, RHOM, CMAX,
     &              CNNB, COST, NRUN, MODEL, TCOMP, RCM, VCM, AZ, RSMIN,
     &              BODY1
   55 FORMAT (' #1',F5.2,F6.1,F7.2,F6.2,I4,F7.3,F6.0,F7.0,2F6.1,F6.2,I6,
     &                                 I7,F8.1,F7.3,F8.4,F9.5,F7.3,F8.4)
*
      WRITE (6,60)
   60 FORMAT (/,7X,'NNPRED   NBCORR  NBFULL  NBVOID  NMTRY  NMERG',
     &             '  NESC   NL')
      WRITE (6,65)  NNPRED, NBCORR, NBFULL, NBVOID, NMTRY, NMERG, NESC,
     &              NLIST(1)
   65 FORMAT (' #2',I10,I9,2I8,2I7,I6,I5)
*
*       Check optional output of single bodies & binaries.
      CALL BODIES
*
*       Update next output time and reset minimum neighbour sphere.
      TNEXT = TNEXT + DELTAT
      RSMIN = 100.0
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
      A(6) = RC
      A(7) = RHOM
      A(8) = CMAX
      A(9) = BE(3)
      A(10) = ERRTOT
      A(11) = Q
      A(12) = 0.0
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
      NK = 12
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
      IF (KZ(18).GT.0) THEN
          CALL CMCORR
      END IF
*
*       Restore current velocities to original values.
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
*       See whether COMMON should be saved on unit 2 after energy check.
      IF (KZ(2).GE.1) CALL MYDUMP(1,2)
*
*       Check termination time.
      IF (TIME.LT.TCRIT) RETURN
*       Save COMMON each time by CALL MYDUMP(1,1) before TIME < TCRIT test.
*
*       Terminate after optional COMMON save.
      WRITE (6,190)  TIME, CPUTOT/60.0, ERRTOT
  190 FORMAT (/,9X,'END RUN',3X,'TIME =',F8.2,'  CPUTOT =',F6.1,
     &                                                '  ERRTOT =',F9.5)
      IF (KZ(1).GT.0) CALL MYDUMP(1,1)
*
      STOP
*
      END
