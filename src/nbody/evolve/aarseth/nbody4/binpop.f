      SUBROUTINE BINPOP
*
*
*       Initial binary population.
*       --------------------------
*
      INCLUDE 'common4.h'
      REAL*8  XORB(2),VORB(2),XREL(3),VREL(3),PX(3),QX(3),BS(NMAX)
      REAL*4  RAN2
      DATA  ETA1,ETA2 /2.5,45.0/
*
*
      READ (5,*)  SEMI0, ECC0, RATIO, RANGE, NSKIP, IDORM, ICIRC
      NBIN = NBIN0
      NBIN1 = NBIN + 1
      WRITE (6,1)  NBIN, SEMI0, ECC0, RATIO, RANGE, NSKIP, IDORM, ICIRC
    1 FORMAT (/,12X,'BINARIES:   NBIN =',I4,'  A =',F9.6,'  E =',F6.2,
     &              '  RATIO =',F4.1,'  RANGE =',F6.1,'  NSKIP =',I3,
     &              '  IDORM =',I2,'  ICIRC =',I2,/)
*
*       Check type of binary mass distribution (NSKIP, IMF2 or split c.m.).
      IF (NSKIP.EQ.0.OR.KZ(20).GE.2) GO TO 10
      IF (RATIO.EQ.1.0) GO TO 20
*
*       Select binaries from the most massive bodies (frequency NSKIP).
      ILAST = (1 + NSKIP)*NBIN
      JSKIP = 0
      JS = 0
      JB = 1
*
*       Transfer binary masses to first NBIN locations.
      DO 6 I = 2,ILAST
          JSKIP = JSKIP + 1
*       Copy binary mass of body #I to new global location.
          IF (JSKIP.GT.NSKIP) THEN
              JSKIP = 0
              JB = JB + 1
              BODY(JB) = BODY(I)
          ELSE
*       Save next NSKIP masses of single bodies.
              JS = JS + 1
              BS(JS) = BODY(I)
          END IF
    6 CONTINUE
*
*       Restore the single bodies in subsequent locations.
      JS = 0
      DO 8 I = NBIN1,ILAST
          JS = JS + 1
          BODY(I) = BS(JS)
    8 CONTINUE
*
*       Move main variables of all single bodies.
   10 DO 15 I = N,NBIN1,-1
          J = I + NBIN
          BODY(J) = BODY(I)
          DO 12 K = 1,3
              X(K,J) = X(K,I)
              XDOT(K,J) = XDOT(K,I)
   12     CONTINUE
   15 CONTINUE
*
*       Create space for each binary component next to primary.
   20 DO 30 I = NBIN,2,-1
          J = 2*I - 1
          BODY(J) = BODY(I)
          DO 25 K = 1,3
              X(K,J) = X(K,I)
              XDOT(K,J) = XDOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
*       Introduce binary components from relative motion.
      DO 60 I = 1,NBIN
*
*       Randomize perihelion, node & inclination (ZI = 0.25 before 3/99).
          PI = TWOPI*RAN2(IDUM1)
          OMEGA = TWOPI*RAN2(IDUM1)
          ZI = 0.5*TWOPI*RAN2(IDUM1)
*
*       Set transformation elements (Brouwer & Clemence p. 35).
          PX(1) = COS(PI)*COS(OMEGA) - SIN(PI)*SIN(OMEGA)*COS(ZI)
          QX(1) =-SIN(PI)*COS(OMEGA) - COS(PI)*SIN(OMEGA)*COS(ZI)
          PX(2) = COS(PI)*SIN(OMEGA) + SIN(PI)*COS(OMEGA)*COS(ZI)
          QX(2) =-SIN(PI)*SIN(OMEGA) + COS(PI)*COS(OMEGA)*COS(ZI)
          PX(3) = SIN(PI)*SIN(ZI)
          QX(3) = COS(PI)*SIN(ZI) 
*
*       Specify component masses (copy BODY0 from IMF2 or use RATIO).
          I1 = 2*I - 1
          I2 = 2*I
          IF (KZ(20).GE.2) THEN
              BODY(I1) = BODY0(I1)
              BODY(I2) = BODY0(I2)
          ELSE IF (RATIO.EQ.1.0) THEN
              BODY(I2) = BODY(I1) 
          ELSE
              BODY(I1) = RATIO*BODY(I1)
              BODY(I2) = BODY(I1)*(1.0 - RATIO)/RATIO
          END IF
*
*       Choose random (thermalized) or fixed eccentricity.
          IF (ECC0.LT.0.0) THEN
              ECC2 = RAN2(IDUM1)
              ECC = SQRT(ECC2)
          ELSE
              ECC = ECC0
          END IF
*
*       Select semi-major axis from uniform distribution in log(A) or SEMI0.
          IF (RANGE.GT.0.0) THEN
              EXP1 = RAN2(IDUM1)*LOG10(RANGE)
              SEMI = SEMI0/10.0**EXP1
          ELSE
              SEMI = SEMI0
          END IF
*
*       Check for eigen-evolution (Pavel Kroupa & Rosemary Mardling).
          IF (ICIRC.NE.0) THEN
              ZMB = (BODY(I1) + BODY(I2))*ZMBAR
*       Include minimum period (copy RANGE; at least 1 day).
              PMIN = MAX(RANGE,1.0D0)
              IT = 0
   35         XR = RAN2(IDUM1)
*       Generate period distribution (Pavel Kroupa: MN 277, 1491, eq.11b).
              P0 = LOG10(PMIN) + SQRT(ETA2*(EXP(2.0*XR/ETA1) - 1.0))
              TK = 10.0**P0
*       Invert eccentricity from thermal distribution (XR = E**2).
              XR = RAN2(IDUM1)
              ES0 = SQRT(XR)
*       Set pericentre distance in AU with period in days & mass in SU.
              RP0 = (1.0 - ES0)*((TK/365.0)**2*ZMB)**0.3333
              RP0 = RP0/RAU
              A0 = RP0/(1.0 - ES0)
              E0 = ES0
*       Define K* = 0/1 and enhanced radii for pre-main sequence.
              KSTAR(I1) = 1
              KSTAR(I2) = 1
              IF (BODY(I1)*ZMBAR.LT.0.7) KSTAR(I1) = 0
              IF (BODY(I2)*ZMBAR.LT.0.7) KSTAR(I2) = 0
              RADIUS(I1) = 5.0*SQRT(BODY(I1)*ZMBAR)/SU
              RADIUS(I2) = 5.0*SQRT(BODY(I2)*ZMBAR)/SU
              IF (RP0.LT.MAX(RADIUS(I1),RADIUS(I2))) THEN
                  WRITE (6,38)  I1, ZMB, ES0, A0, RP0
   38             FORMAT (12X,'COLLISION:    I1 MB E A RP ',
     &                                       I6,F6.1,F7.3,1P,2E10.2)
                  ICOLL = ICOLL + 1
                  GO TO 35
              END IF
*       Perform eigen-evolution of pericentre & eccentricity for 10^6 yrs.
              TC = -1.0/TSCALE
              CALL TCIRC(RP0,ES0,I1,I2,ICIRC,TC)
*       Copy modified eccentricity and re-evaluate the semi-major axis.
              ECC = ES0
              SEMI = RP0/(1.0 - ECC)
              IT = IT + 1
              IF (SEMI.GT.SEMI0.AND.IT.LT.25) GO TO 35
              TK = 365.0*SQRT((SEMI*RAU)**3/ZMB)
              IF (ECC.LE.0.002) IC0 = IC0 + 1
              IF (TK.LT.PMIN) IC1 = IC1 + 1
              IF (TK.LT.2.0*PMIN) IC2 = IC2 + 1
              IF (TK.LT.5.0*PMIN) IC3 = IC3 + 1
              WRITE (23,40)  IT, I1, ZMB, E0, ECC, A0, SEMI, TK
   40         FORMAT (12X,'BINARY:   IT I1 MB E0 E A0 A P ',
     &                               I2,I5,F5.1,2F7.3,1P,3E10.2)
              CALL FLUSH(23)
          ELSE IF (KZ(27).EQ.1) THEN
*       Obtain tidal encounter distance (4*RADIUS) from square root relation.
              RSUN = 1.0/SU
              ZM = MAX(BODY(I1),BODY(I2))*ZMBAR
              RT = 4.0*RSUN*SQRT(ZM)
*       Modify orbital elements to avoid early tidal interaction.
   42         IF (SEMI*(1.0 - ECC).LT.2.0*RT) THEN
*       Increase semi-major axis or reduce eccentricity until peri > 2*RT.
   44             IF (SEMI.LT.2.0*RT) THEN
                      SEMI = 2.0*SEMI
                      GO TO 44
                  ELSE
                      ECC = 0.9*ECC
                  END IF
                  WRITE (17,46)  I1, I2, ECC, SEMI, SEMI*(1.0-ECC), RT
   46             FORMAT (12X,'REDUCE ECC:    I1 I2 E A PM RT ',
     &                                        2I5,F7.3,1P,3E10.2)
                  CALL FLUSH(17)
                  GO TO 42
              END IF
          END IF
*
*       Specify relative motion at apocentre and sum binding energy.
          XORB(1) = SEMI*(1.0 + ECC)
          XORB(2) = 0.0
          VORB(1) = 0.0
          ZMBIN = BODY(I1) + BODY(I2)
          VORB(2) = SQRT(ZMBIN*(1.0D0 - ECC)/(SEMI*(1.0D0 + ECC)))
          EBIN0 = EBIN0 - 0.5*BODY(I1)*BODY(I2)/SEMI
*
*       Transform to relative variables.
          DO 50 K = 1,3
              XREL(K) = PX(K)*XORB(1) + QX(K)*XORB(2)
              VREL(K) = PX(K)*VORB(1) + QX(K)*VORB(2)
   50     CONTINUE
*
*       Set global variables for each component.
          DO 55 K = 1,3
              X(K,I1) = X(K,I1) + BODY(I2)*XREL(K)/ZMBIN
              X(K,I2) = X(K,I1) - XREL(K)
              XDOT(K,I1) = XDOT(K,I1) + BODY(I2)*VREL(K)/ZMBIN
              XDOT(K,I2) = XDOT(K,I1) - VREL(K)
   55     CONTINUE
   60 CONTINUE
*
*       Update the total particle number after primary splitting or IMF2.
      IF (RATIO.LT.1.0.OR.KZ(20).GE.2) THEN
          N = N + NBIN
          NZERO = N
          NTOT = N
          BODYM = ZMASS/FLOAT(N)
          IF (NSKIP.GT.0) THEN
              WRITE (6,62)  (BODY(J),J=1,10)
              WRITE (6,64)  (BODY(J),J=2*NBIN+1,2*NBIN+10)
   62         FORMAT (/,12X,'BINARY MASSES (1-10):  ',10F9.5)
   64         FORMAT (/,12X,'SINGLE MASSES (1-10):  ',10F9.5,/) 
          END IF
      END IF
*
*       Include procedure for introducing dormant binaries.
      IF (IDORM.GT.0) THEN
          DO 66 I = 1,NBIN
              I1 = 2*I - 1
              I2 = I1 + 1
              ZMBIN = BODY(I1) + BODY(I2)
              DO 65 K = 1,3
                  X(K,I) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZMBIN
                  XDOT(K,I) = (BODY(I1)*XDOT(K,I1) +
     &                         BODY(I2)*XDOT(K,I2))/ZMBIN
   65         CONTINUE
              BODY(I) = ZMBIN
   66     CONTINUE
*
*       Move the original single particles up to form compact array.
          I1 = 2*NBIN + 1
          I2 = NBIN
          DO 68 I = I1,N
              I2 = I2 + 1
              BODY(I2) = BODY(I)
              DO 67 K = 1,3
                  X(K,I2) = X(K,I)
                  XDOT(K,I2) = XDOT(K,I)
   67         CONTINUE
   68     CONTINUE
*
*       Reset particle membership and turn off binary output option (if = 1).
          N = N - NBIN
          NZERO = N
          NTOT = N
          NBIN0 = 0
          EBIN0 = 0.0
          IF (KZ(8).EQ.1) KZ(8) = 0
      END IF
*
*       Set coordinates & velocities in c.m. rest frame.
      DO 70 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   70 CONTINUE
*
      DO 80 I = 1,N
          DO 75 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   75     CONTINUE
   80 CONTINUE
*
      DO 90 I = 1,N
          DO 85 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
   85     CONTINUE
   90 CONTINUE
*
      RETURN
*
      END
