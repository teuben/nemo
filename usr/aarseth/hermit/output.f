      SUBROUTINE OUTPUT
*
*
*       Energy check and output.
*       ------------------------
*
      INCLUDE 'commonp.h'
      REAL*8  XS(3),XSD(3),SEMI(NMAX),ECC(NMAX)
*
*
*       Predict X & XDOT for all particles.
      CALL XVPRED(1,N)
*
*       Obtain the total energy at current time.
      CALL ENERGY
*
*       Initialize c.m. terms.
      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
*
*       Obtain c.m. & angular momentum integrals and Z-moment of inertia.
      AZ = 0.0D0
      ZM = 0.0D0
      ZMASS = 0.0D0
      DO 20 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 15 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   15     CONTINUE
          AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
   20 CONTINUE
*
*       Form c.m. coordinates & velocities (vectors & scalars).
      DO 25 K = 1,3
          CMR(K) = CMR(K)/ZMASS
          CMRDOT(K) = CMRDOT(K)/ZMASS
   25 CONTINUE
*
      CMR(4) = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)
      CMRDOT(4) = SQRT(CMRDOT(1)**2 + CMRDOT(2)**2 + CMRDOT(3)**2)
*
*       Define crossing time and save single particle energy.
      ETOT = ZKIN - POT
      TCR = ZMASS**2.5/(2.0*ABS(ETOT))**1.5
*
*       Update energies and form the relative error (divide by ZKIN or ETOT).
      IF (TIME.LE.0.0D0) THEN
          DE = 0.0D0
          BE(1) = ETOT
          BE(3) = ETOT
      ELSE
          BE(2) = BE(3)
          BE(3) = ETOT
          DE = BE(3) - BE(2)
          DETOT = DETOT + DE
          DE = DE/MAX(ZKIN,ABS(ETOT))
          ERRTOT = ERRTOT + DE
      END IF
*
*       Print diagnostic information.
      WRITE (6,40)  TIME, NSTEPS, BE(3), DE, AZ
   40 FORMAT (/,' YRS =',1P,E9.1,'  # =',0P,I10,'  E =',F10.6,
     &          '  DE =',1P,E10.2,'  AZ =',F12.8)
*
      TPRINT = TPRINT + DELTAT
      CALL CPUTIM(TCOMP)
      CPUTOT = CPUTOT + TCOMP - CPU0
      CPU0 = TCOMP
*
*       Save COMMON after energy check.
      TDUMP = TIME
      IF (KZ(2).GE.1) CALL MYDUMP(1,2)
*
*       Check termination criterion.
      IF (TIME.GE.TCRIT) THEN
*       Terminate after optional COMMON save.
          WRITE (6,60)  TIME, CPUTOT/60.0, ERRTOT, DETOT
   60     FORMAT (//,9X,'END RUN',3X,'TIME =',F7.1,'  CPUTOT =',F7.1,
     &                  '  ERRTOT =',1P,E10.2,'  DETOT =',E10.2)
          IF (KZ(1).GT.0) CALL MYDUMP(1,1)
          STOP
      END IF
*
      RETURN
*
      END
