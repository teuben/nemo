      SUBROUTINE SCALE
*
*
*       Scaling to new units.
*       ---------------------
*
      INCLUDE 'common4.h'
*
*
*       Read virial ratio, rotation scaling factors & boundary radius.
      READ (5,*)  Q, VXROT, VZROT, RTIDE
*
      ZMASS = 0.0D0
      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
*
*       Form total mass and centre of mass displacements.
      RIJ2 = 0.0
      DO 30 I = 1,N
          ZMASS = ZMASS + BODY(I)
          RIJ = 0.0
          DO 25 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
              RIJ = RIJ + X(K,I)*X(K,I)
   25     CONTINUE
          RIJ2 = MAX(RIJ2,RIJ)
   30 CONTINUE
      RIJ2 = SQRT(RIJ2)
      IF(KZ(22).EQ.4) RTIDE = RIJ2
*
*       Adjust coordinates and velocities to c.m. rest frame.
      IF(KZ(22).GT.0)THEN
          DO 40 I = 1,N
              DO 35 K = 1,3
                  X(K,I) = X(K,I) - CMR(K)/ZMASS
                  XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
                  F(K,I) = 0.D0
                  FDOT(K,I) = 0.D0
                  D2(K,I) = 0.D0
   35         CONTINUE
   40     CONTINUE
      ENDIF
      IF (KZ(22).GT.2.OR.KZ(5).EQ.3) GO TO 52
*
*       Scale masses to standard units of <M> = 1/N and set total mass.
      DO 50 I = 1,N
          BODY(I) = BODY(I)/ZMASS
   50 CONTINUE
      ZMASS = 1.0
*
*       Obtain the total kinetic & potential energy.
   52 CALL ENERGY
*
*       Use generalized virial theorem for external tidal field.
      IF (KZ(14).GT.0) THEN
          AZ = 0.0D0
          DO 55 I = 1,N
              AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
   55     CONTINUE
          IF (KZ(14).EQ.1) THEN
*       Use Chandrasekhar eq. (5.535) for virial ratio (rotating frame only).
              VIR = POT - 2.0*(ETIDE + 0.5*TIDAL(4)*AZ)
          ELSE
              VIR = POT - 2.0*ETIDE
          END IF
      ELSE
          VIR = POT
      END IF
*
*       Allow two optional ways of skipping standard velocity scaling.
      IF (KZ(22).GT.2.OR.KZ(5).EQ.2.OR.KZ(5).EQ.3) THEN
          QV = SQRT(Q*VIR/ZKIN)
          E0 = ZKIN*QV**2 - POT + ETIDE
          SX = 1.0
*       Rescale velocities to new masses for two Plummer spheres.
          IF (KZ(5).EQ.2) THEN
              ZKIN = 0.0
              DO 57 I = 1,N
                  DO 56 K = 1,3
                      XDOT(K,I) = XDOT(K,I)*QV
                      ZKIN = ZKIN + 0.5*BODY(I)*XDOT(K,I)**2
   56             CONTINUE
   57         CONTINUE
              E0 = ZKIN - POT + ETIDE
              Q = ZKIN/POT
              WRITE (6,59)  E0, ZKIN/POT
   59         FORMAT (/,12X,'UNSCALED ENERGY    E =',F10.6,
     &                                       '  Q =',F6.2)
          ELSE
              IF (KZ(5).EQ.3) E0 = ZKIN - POT
              WRITE (6,54)  E0
   54         FORMAT (/,12X,'UNSCALED ENERGY    E =',F10.6)
          END IF
      ELSE
*       Scale non-zero velocities by virial theorem ratio.
          IF (ZKIN.GT.0.0D0) THEN
              QV = SQRT(Q*VIR/ZKIN)
              DO 60 I = 1,N
                  DO 58 K = 1,3
                      XDOT(K,I) = XDOT(K,I)*QV
   58             CONTINUE
   60         CONTINUE
          ELSE
              QV = 1.0
          END IF
*
*       Scale total energy to standard units (E = -0.25 for Q < 1).
          E0 = -0.25
*       Include case of hot system inside reflecting boundary.
          IF (KZ(29).GT.0.AND.Q.GT.1.0) E0 = 0.25
*         ETOT = (Q - 1.0)*POT
          ETOT = ZKIN*QV**2 - POT + ETIDE
*       Note that final ETOT will differ from -0.25 since ETIDE = 0.
          IF (Q.LT.1.0) THEN
              SX = E0/ETOT
          ELSE
              SX = 1.0
          END IF
*
          WRITE (6,65)  SX, ETOT, BODY(1), BODY(N), ZMASS/FLOAT(N)
   65     FORMAT (//,12X,'SCALING:    SX =',F6.2,'  E =',1PE10.2,
     &                   '  M(1) =',E9.2,'  M(N) =',E9.2,'  <M> =',E9.2)
*
*       Scale coordinates & velocities to the new units.
          RIJ2 = 0.0
          DO 70 I = 1,N
              RIJ = 0.0
              DO 68 K = 1,3
                  X(K,I) = X(K,I)/SX
                  XDOT(K,I) = XDOT(K,I)*SQRT(SX)
                  RIJ = RIJ + X(K,I)*X(K,I)
   68         CONTINUE
              RIJ2 = MAX(RIJ2,RIJ)
   70     CONTINUE
          RIJ2 = SQRT(RIJ2)
      ENDIF
*
*       Scale tidal radius given by King model.
      IF (KZ(22).GE.4) THEN
*       In this case Rtide is actually the maximum radius for the model. 
          RBAR = RTIDE/RIJ2
          WRITE (6,72) RBAR
 72       FORMAT (/,12X,'KING MODEL: NEW RBAR =',F6.2)
          IF(KZ(14).EQ.3) RTIDE = 1.2*RIJ2
      ELSEIF (KZ(14).GT.3.AND.KZ(22).NE.2) THEN
          RTIDE = 1.2*RIJ2
      ELSEIF (KZ(14).GT.3) THEN
          RTIDE = RTIDE/SX
*         RBAR = SX
      END IF
*       Radius of reflecting sphere for hot system ( KZ(29) )
      RSPH2 = RTIDE
*
*       Check whether to include rotation (VXROT = 0 in standard case). 
      IF (VXROT.GT.0.0D0) THEN
*       Set angular velocity for retrograde motion (i.e. star clusters).
          OMEGA = -SX*SQRT(ZMASS*SX)
          WRITE (6,75)  VXROT, VZROT, OMEGA
   75     FORMAT (/,12X,'ROTATION:    VXROT =',F6.2,'  VZROT =',F6.2,
     &                                                 '  OMEGA =',F7.2)
*
*       Add solid-body rotation about Z-axis (reduce random velocities).
          DO 80 I = 1,N
              XDOT(1,I) = XDOT(1,I)*VXROT - X(2,I)*OMEGA
              XDOT(2,I) = XDOT(2,I)*VXROT + X(1,I)*OMEGA
              XDOT(3,I) = XDOT(3,I)*VZROT
   80     CONTINUE
      END IF
*
*       Check option for writing the initial conditions on unit 10.
      IF (KZ(22).EQ.1.AND.KZ(5).LE.1) THEN
          DO 85 I = 1,N
              WRITE (10,84)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
   84         FORMAT (1P,7E14.6)
   85     CONTINUE
      END IF
*
*       Check option for reading initial subsystems.
      IF (KZ(24).GT.0) THEN
          K = KZ(24)
          DO 88 I = 1,K
              READ (5,*)  (X(J,I),J=1,3), (XDOT(J,I),J=1,3)
   88     CONTINUE
      END IF
*
*       Set initial crossing time in scaled units.
   90 TCR = ZMASS**2.5/(2.0D0*ABS(E0))**1.5
      TCR0 = TCR
*
*       Form approximate half-mass radius and define hard binary energy.
      RSCALE = 0.4*ZMASS**2/(SX*POT)
      EBH = -0.5*ZMASS*ECLOSE/FLOAT(N)
*       Set square radius of reflecting sphere (scaled by new RSCALE).
      RSPH2 = (RSPH2*RSCALE)**2
*       Form equilibrium rms velocity (temporarily defined as VC).
      VC = SQRT(2.0D0*ABS(E0)/ZMASS)
*
*       Check for general binary search of initial condition.
      IF (KZ(4).GT.0) THEN
          CALL EVOLVE(0,0)
      END IF
*
*       Obtain actual half-mass radius (suppress, done in adjust.f).
*     IF (KZ(7).GT.0) THEN
*         CALL LAGR(RDENS)
*     END IF
*
*       Print half-mass relaxation time & equilibrium crossing time.
      A1 = FLOAT(N)
      TRH = 4.0*TWOPI/3.0*(VC*RSCALE)**3/(15.4*ZMASS**2*LOG(0.4*A1)/A1)
      WRITE (6,95)  TRH, TCR, 2.0*RSCALE/VC
   95 FORMAT (/,12X,'TIME SCALES:    TRH =',1PE8.1,'  TCR =',E8.1,
     &                                            '  2<R>/<V> =',E8.1,/)
*
      RETURN
*
      END
