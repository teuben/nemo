      SUBROUTINE OUTPUT
*
*
*       Output and energy check.
*       ------------------------
*
      INCLUDE 'common1.h'
      COMMON/VELOC/  X2DOT(3,NMAX)
      SAVE AZ0
*
*
*       Obtain the total kinetic & potential energy at current time.
      CALL ENERGY
*
*       Form virial theorem ratio (NB! VIR differs from POT if EPS2 > 0).
      Q = ZKIN/POT
*
*       Define the standard crossing time.
      ETOT = ZKIN - POT
      TCR = ZMASS**2.5/(2.0*ABS(ETOT))**1.5
*
*       Update energies and form the relative error using MAX(ZKIN,POT).
      IF (TIME.EQ.0.0D0) THEN
          DE = 0.0
          BE(1) = ETOT
          BE(3) = ETOT
      ELSE
          BE(2) = BE(3)
          BE(3) = ETOT
          DE = (BE(3) - BE(2))/BE(3)
      END IF
*
*       Initialize c.m. integrals.
      DO 10 K = 1,3
          CMR(K) = 0.0
          CMRDOT(K) = 0.0
   10 CONTINUE
*
*       Obtain Z angular momentum & c.m. integrals.
      AZ = 0.0
      DO 30 I = 1,N
          AZ = AZ + BODY(I)*(X(1,I)*X2DOT(2,I) - X(2,I)*X2DOT(1,I))
          DO 25 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*X2DOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
      RCM = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)/ZMASS
      VCM = SQRT(CMRDOT(1)**2 + CMRDOT(2)**2 + CMRDOT(3)**2)/ZMASS
*
      TC = TIME/TCR
      WRITE (6,40)  TIME, Q, NSTEPI, DE, BE(3), TC
   40 FORMAT (/,' T =',F6.1,'  Q =',F5.2,'  STEPS =',I7,'  DE =',F10.6,
     &                                      '  E =',F10.6,'  TC =',F6.1)
      IF (TIME.EQ.0.0D0) AZ0 = AZ
      ERROR = (BE(3) - BE(1))/BE(1)
      WRITE (6,42)  RCM, VCM, ERROR, (AZ-AZ0)
   42 FORMAT (' ERRORS    RCM VCM DE/E DZ   ',1P,4E10.2)
*
*       Include diagnostics for the closest binary.
      RIJ2 = 0.0
      VIJ2 = 0.0
      RDOT = 0.0
      DO 45 K = 1,3
          RIJ2 = RIJ2 + (X(K,ICL) - X(K,JCL))**2
          VIJ2 = VIJ2 + (X2DOT(K,ICL) - X2DOT(K,JCL))**2
          RDOT = RDOT + (X(K,ICL)-X(K,JCL))*(X2DOT(K,ICL)-X2DOT(K,JCL))
   45 CONTINUE
      RIJ = SQRT(RIJ2)
      SEMI = 2.0/RIJ - VIJ2/(BODY(ICL) + BODY(JCL))
      SEMI = 1.0/SEMI
      ECC2 = (1.0 - RIJ/SEMI)**2 + RDOT**2/(SEMI*(BODY(ICL)+BODY(JCL)))
      IF ((SEMI.GT.0.0.AND.SEMI.LT.0.1).OR.N.EQ.2) THEN
          WRITE (6,50)  ICL, JCL, SQRT(ECC2), SEMI
   50     FORMAT (' BINARY    I J E A  ',2I5,F8.4,F9.4)
      END IF
      CALL FLUSH(6)
*
*       Update next output time.
      TNEXT = TNEXT + DELTAT
*
      IF (TIME.GT.TCRIT) STOP
*
      RETURN
*
      END
