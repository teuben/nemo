*       Program NBODY0. Hermite block-steps. Coded by Sverre Aarseth 10/06.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON  X(3,50),XDOT(3,50),BODY(50),D2(3,50),D3(3,50),EPS2,N
      REAL*8  X0(3,50),X0DOT(3,50),T0(50),STEP(50),F(3,50),FDOT(3,50),
     &        FIRR(3),FD(3)
      INTEGER NEXT(50)
      DATA TIME,TNEXT,SMAX,TMIN,NSTEPS /0.0D0,0.0D0,64.0D0,1.0D+10,0/
      READ (5,*)  N, ETA, DELTAT, TCRIT, EPS2
      READ (5,*)  (BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3),I=1,N)
*       Initialize forces and time-steps.
      DO 5 J = 1,N
          CALL FFDOT(J,FIRR,FD)
          DO 4 K = 1,3
              F(K,J) = 0.5*FIRR(K)
              FDOT(K,J) = FD(K)/6.0
              X0(K,J) = X(K,J)
              X0DOT(K,J) = XDOT(K,J)
    4     CONTINUE
          STEP(J) = SQRT(ETA/SQRT(FIRR(1)**2 + FIRR(2)**2 + FIRR(3)**2))
*       Create a truncated block time-step according to STEP = 1/2**(K-1).
          K = 2 + LOG(SMAX/STEP(J))/LOG(2.0D0)
          STEP(J) = SMAX/2.0D0**(K-1)
          TMIN = MIN(STEP(J),TMIN)
          T0(J) = 0.0D0
    5 CONTINUE
*       Obtain the total energy and generate some output.
   10 EK = 0.0D0
      POT = 0.0D0
      DO 15 I = 1,N
          WRITE (6,12)  BODY(I),(X(K,I),K=1,3),(XDOT(K,I),K=1,3),STEP(I)
   12     FORMAT (' M R V STEP ',F8.3,2X,3F8.3,2X,3F7.2,F10.6)
          EK = EK + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
          DO 14 J = I+1,N
              POT = POT + BODY(I)*BODY(J)/SQRT((X(1,I) - X(1,J))**2 +
     &              (X(2,I) - X(2,J))**2 + (X(3,I) - X(3,J))**2 + EPS2)
   14     CONTINUE
   15 CONTINUE
      WRITE (6,18)  TIME, NSTEPS, 0.5*EK - POT
   18 FORMAT (' TIME =',F6.2,'  NSTEPS =',I6,'  ENERGY =',F12.8)
      IF (TIME.GT.TCRIT) STOP
      TNEXT = TNEXT + DELTAT
*       Determine block of particles to be advanced at TMIN and set new time.
   20 LENGTH = 0
      DO 25 J = 1,N
          IF (T0(J) + STEP(J).EQ.TMIN) THEN
              LENGTH = LENGTH + 1
              NEXT(LENGTH) = J
          END IF
   25 CONTINUE
      TIME = TMIN
*       Predict all coordinates and velocities to order FDOT.
      DO 30 J = 1,N
          S = TIME - T0(J)
          DO 28 K = 1,3
          X(K,J) = ((FDOT(K,J)*S + F(K,J))*S + X0DOT(K,J))*S + X0(K,J)
          XDOT(K,J) = (FDOT(K,J)*1.5*S + F(K,J))*2.0*S + X0DOT(K,J)
   28     CONTINUE
   30 CONTINUE
      TMIN = 1.0D+10
      DO 50 L = 1,LENGTH
          I = NEXT(L)
*       Evaluate new force and first derivative for body #I.
          CALL FFDOT(I,FIRR,FD)
*       Include Hermite corrector and set new F/2, FDOT/6, D2, & D3.
          DT = TIME - T0(I)
          T0(I) = TIME
          DO 40 K = 1,3
	      DF = 2.0*F(K,I) - FIRR(K)
	      FID = 6.0*FDOT(K,I)
	      SUM = FID + FD(K)
	      AT3 = 2.0*DF + DT*SUM
	      BT2 = -3.0*DF - DT*(SUM + FID)
              X0(K,I) = (0.6*AT3 + BT2)*DT**2/12.0 + X(K,I)
              X0DOT(K,I) = (0.75*AT3 + BT2)*DT/3.0 + XDOT(K,I)
	      F(K,I) = 0.5*FIRR(K)
	      FDOT(K,I) = FD(K)/6.0
              D2(K,I) = (3.0*AT3 + BT2)*2.0/DT**2
              D3(K,I) = AT3*6.0/DT**3
   40     CONTINUE
*       Form new time-step from force derivatives and check change.
          FI2 = FIRR(1)**2 + FIRR(2)**2 + FIRR(3)**2
          FD2 = FD(1)**2 + FD(2)**2 + FD(3)**2
          FD22 = D2(1,I)**2 + D2(2,I)**2 + D2(3,I)**2
          FD32 = D3(1,I)**2 + D3(2,I)**2 + D3(3,I)**2
          SI = SQRT(ETA*(SQRT(FI2*FD22) + FD2)/(SQRT(FD2*FD32) + FD22))
          IF (SI.GT.2.0*STEP(I).AND.DMOD(TIME,2.0*STEP(I)).EQ.0.D0) THEN
              STEP(I) = MIN(2.0*STEP(I),SMAX)
          ELSE IF (SI.LT.STEP(I)) THEN
              STEP(I) = 0.5*STEP(I)
          END IF
          TMIN = MIN(T0(I) + STEP(I),TMIN)
          NSTEPS = NSTEPS + 1
   50 CONTINUE
      IF (TIME.LT.TNEXT) GO TO 20
      GO TO 10
      END
      SUBROUTINE FFDOT(I,FIRR,FD)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON  X(3,50),XDOT(3,50),BODY(50),D2(3,50),D3(3,50),EPS2,N
      REAL*8  FIRR(3),FD(3),A(3),DV(3)
      DO 1 K = 1,3
          FIRR(K) = 0.0D0
          FD(K) = 0.0D0
    1 CONTINUE
      DO 20 J = 1,N
          IF (J.EQ.I) GO TO 20
          DO 5 K = 1,3
              A(K) = X(K,J) - X(K,I)
              DV(K) = XDOT(K,J) - XDOT(K,I)
    5     CONTINUE
          RIJ2 = A(1)**2 + A(2)**2 + A(3)**2 + EPS2
          DR3I = BODY(J)/(RIJ2*SQRT(RIJ2))
          DRDV = 3.0*(A(1)*DV(1) + A(2)*DV(2) + A(3)*DV(3))/RIJ2
          DO 10 K = 1,3
              FIRR(K) = FIRR(K) + A(K)*DR3I
              FD(K) = FD(K) + (DV(K) - A(K)*DRDV)*DR3I
   10     CONTINUE
   20 CONTINUE
      RETURN
      END
