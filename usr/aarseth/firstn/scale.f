      SUBROUTINE SCALE
*
*
*       Scaling to new units.
*       ---------------------
*
      INCLUDE 'common1.h'
      REAL*4  A(3)
*
*
*       Read virial ratio and scaling option.
      READ (5,*)  Q, ISCALE
*
      ZMASS = 0.0
      DO 10 K = 1,3
          CMR(K) = 0.0
          CMRDOT(K) = 0.0
   10 CONTINUE
*
*       Form total mass and centre of mass displacements.
      DO 30 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 25 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   25     CONTINUE
   30 CONTINUE
*
*       Adjust coordinates and velocities to c.m. rest frame.
      DO 40 I = 1,N
          DO 35 K = 1,3
              X(K,I) = X(K,I) - CMR(K)/ZMASS
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)/ZMASS
   35     CONTINUE
   40 CONTINUE
*
*       Scale masses to standard units of <M> = 1/N.
      DO 50 I = 1,N
          BODY(I) = BODY(I)/ZMASS
   50 CONTINUE
      ZMASS = 1.0
*
*       Obtain the total kinetic & potential energy.
      CALL ENERGY
*
*       Scale non-zero velocities by virial theorem ratio (assume VIR = POT).
      IF (ZKIN.GT.0.0) THEN
          QV = SQRT(Q*POT/ZKIN)
          DO 55 I = 1,N
              DO 52 K = 1,3
                  XDOT(K,I) = XDOT(K,I)*QV
   52         CONTINUE
   55     CONTINUE
      END IF
*
*       Scale total energy to standard units E = -0.25 (skip if ISCALE = 0).
      ETOT = (Q - 1.0)*POT
*
*       Check scaling option.
      IF (ISCALE.EQ.0) THEN
          E0 = ETOT
          GO TO 80
      END IF
*
*       Define scaling factor (set E0 = ETOT if energy scaling not desired).
      E0 = -0.25
      SX = E0/ETOT
*
      WRITE (6,65)  SX, ETOT, BODY(1), BODY(N), ZMASS/FLOAT(N)
   65 FORMAT (//,12X,'SCALING:   SX  =',F6.2,'  E =',1PE10.2,
     &                 '  M(1) =',E9.2,'  M(N) =',E9.2,'  <M> =',E9.2,/)
*
*       Scale coordinates & velocities to the new units.
      DO 70 I = 1,N
          DO 68 K = 1,3
              X(K,I) = X(K,I)/SX
              XDOT(K,I) = XDOT(K,I)*SQRT(SX)
   68     CONTINUE
   70 CONTINUE
*
*       Set initial crossing time in scaled units.
   80 TCR = ZMASS**2.5/(2.0*ABS(E0))**1.5
*
*       Obtain forces and assign a small step to reduce error.
      CALL FORCE
      STEP = 0.2*STEP
*
*       Initialize force differences and copy b1.
      DO 90 I = 1,N
          DO 85 K = 1,3
              a1(K,I) = 0.0
              d2(K,I) = 0.0
              b1(K,I) = F(K,I)
   85     CONTINUE
   90 CONTINUE
*
*       Specify the velocities at T = -STEP/2 (half the initial step).
      h1 = 0.5*STEP
      DO 100 I = 1,N
          DO 95 K = 1,3
              XDOT(K,I) = XDOT(K,I) - h1*F(K,I)
   95     CONTINUE
  100 CONTINUE
*
*       Scale output time interval & termination time by initial TCR.
      DELTAT = DELTAT*TCR
      TCRIT = TCRIT*TCR
*
      RETURN
*
      END
