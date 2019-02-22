      SUBROUTINE SCALE
*
*
*       Scaling to new units.
*       ---------------------
*
      INCLUDE 'common2.h'
      REAL*4  A(3)
*
*
*       Read virial ratio and scaling factors.
      READ (5,*)  Q, VXROT, VZROT, RBAR, ZMBAR
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
*       Check option for no scaling of initial conditions.
      IF (KZ(16).GT.0) THEN
          CALL ENERGY
          E0 = ZKIN - POT - ETIDE
          GO TO 90
      END IF
*
*       Scale masses to standard units of <M> = 1/N and redefine total mass.
      DO 50 I = 1,N
          BODY(I) = BODY(I)/ZMASS
   50 CONTINUE
      ZMASS = 1.0
*
*       Obtain the total kinetic & potential energy (including virial).
      CALL ENERGY
*
*       Scale non-zero velocities by virial theorem ratio.
      IF (ZKIN.GT.0.0) THEN
          QV = SQRT(Q*VIR/ZKIN)
          ZKIN = ZKIN*QV**2
          DO 55 I = 1,N
              DO 54 K = 1,3
                  XDOT(K,I) = XDOT(K,I)*QV
   54         CONTINUE
   55     CONTINUE
      END IF
*
*       Specify the total energy in standard units (E = -0.25).
      E0 = -0.25
      ETOT = ZKIN - POT - ETIDE
*       Retain old energy with external potential (equation not homologous).
      IF (KZ(15).GT.0) THEN
          E0 = ETOT
      END IF
*       Define scaling factor (set E0 = ETOT if energy scaling not desired).
      SX = E0/ETOT
*       Scale EPS2 by SX**2 if exact energy of -0.25 is desired.
*     EPS2 = EPS2/SX**2
*
      WRITE (6,65)  SX, ETOT, BODY(1), BODY(N), ZMASS/FLOAT(N)
   65 FORMAT (/,12X,'SCALING:   SX =',F6.2,'  E =',1P,E10.2,
     &                   '  M(1) =',E9.2,'  M(N) =',E9.2,'  <M> =',E9.2)
*
*       Scale coordinates & velocities to the new units.
      DO 70 I = 1,N
          DO 68 K = 1,3
              X(K,I) = X(K,I)/SX
              XDOT(K,I) = XDOT(K,I)*SQRT(SX)
   68     CONTINUE
   70 CONTINUE
*
*       Check whether to include rotation (VXROT = 0 in standard case). 
      IF (VXROT.GT.0.0) THEN
*       Set angular velocity for retrograde motion (i.e. star clusters).
          OMEGA = -SX*SQRT(ZMASS*SX)
          WRITE (6,75)  VXROT, VZROT, OMEGA
   75     FORMAT (/,12X,'VXROT =',F6.2,'  VZROT =',F6.2,
     &                                                 '  OMEGA =',F7.2)
*
*       Add solid-body rotation about Z-axis (reduce random velocities).
          DO 80 I = 1,N
              XDOT(1,I) = XDOT(1,I)*VXROT - X(2,I)*OMEGA
              XDOT(2,I) = XDOT(2,I)*VXROT + X(1,I)*OMEGA
              XDOT(3,I) = XDOT(3,I)*VZROT
   80    CONTINUE
      END IF
*
*       Check option for writing the initial conditions on unit 4 (REAL*4).
      IF (KZ(4).EQ.1) THEN
          DO 85 I = 1,N
              DO 82 K = 1,3
                  A(K) = X(K,I)
   82         CONTINUE
              WRITE (4)  BODY(I), (A(K),K=1,3), (XDOT(K,I),K=1,3)
   85     CONTINUE
      END IF
*
*       Set initial crossing time in scaled units.
   90 TCR = ZMASS**2.5/(2.0*ABS(E0))**1.5
      IF(KZ(19).EQ.1) TCR=1
      TCR0 = TCR
*
*       Scale output time interval & termination time by initial TCR.
      DELTAT = DELTAT*TCR
      TCRIT = TCRIT*TCR
*
*       Define GM in cgs units and length scale in pc.
      GM = 6.67E-08*1.989E+33
      PC = 3.0857E+18
*
*       Form time scale in seconds and velocity scale in km/sec.
      TSTAR = SQRT(PC/GM)*PC
      VSTAR = 1.0E-05*SQRT(GM/PC)
*
*       Convert time scale from units of seconds to million years.
      TSTAR = TSTAR/(3.15E+07*1.0E+06)
*
*       Ensure ZMBAR & RBAR > 0.
      IF (ZMBAR.LE.0.0) ZMBAR = 1.0
      IF (RBAR.LE.0.0) RBAR = 1.0
*
*       Rescale solar mass conversion factor by current mean mass.
      ZMBAR = ZMBAR*FLOAT(N)/ZMASS
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      TSTAR = TSTAR*SQRT(RBAR**3/(ZMASS*ZMBAR))
      VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
*       Physical scaling: X, M, V, T from RBAR*X, ZMBAR*M, VSTAR*V, TSTAR*T.
*
      WRITE (6,100)  RBAR, ZMBAR, VSTAR, TSTAR
  100 FORMAT (/,12X,'SCALING PARAMETERS:   R* =',1P,E9.2,
     &                      '  M* =',E9.2,'  V* =',E9.2,'  T* =',E9.2,/)
*
      RETURN
*
      END
