      SUBROUTINE ESCAPE
*
*
*       Escaper detection. 
*       ------------------
*
      INCLUDE 'common1.h'
      INTEGER  JLIST(20)
      REAL*4  RDENS(3)
*
*
*       Adopt twice the tidal radius as escape condition.
      RESC2 = 4.0*RTIDE**2
      RTIDE2 = RTIDE**2
      NCORR = 0
      NCRIT1 = 0
      DO 1 K = 1,3
          RDENS(K) = 0.0
          CMR(K) = 0.0
          CMRDOT(K) = 0.0
    1 CONTINUE
*
      I = 1
*       Set the distance (squared) with respect to the density centre.
    5 RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 + 
     &                               (X(3,I) - RDENS(3))**2
      IF (RI2.LT.RTIDE2) NCRIT1 = NCRIT1 + 1
*
*       See whether escape is indicated.
      IF (RI2.GT.RESC2) THEN
*       Find distance to the nearest neighbour and calculate potential.
          RJMIN2 = 1.0E+06
          POTI = 0.0
          DO 8 J = 1,N
              IF (J.EQ.I) GO TO 8
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              POTI = POTI + BODY(J)/SQRT(RIJ2 + EPS2)
              IF (RIJ2.LT.RJMIN2) THEN
                  RJMIN2 = RIJ2
                  JMIN = J
              END IF
    8     CONTINUE
          VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
*       Check escape criterion for isolated system.
          EI = 0.5*VI2 - POTI
          IF (EI.GT.0.0) GO TO 30
      END IF
*
   10 I = I + 1
   11 IF (I.LE.N) GO TO 5
*
      IF (NCORR.EQ.0) GO TO 25
*
*       Form centre of mass terms.
      DO 14 I = 1,N
          DO 13 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)/ZMASS
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)/ZMASS
   13     CONTINUE
   14 CONTINUE
      CMR(4) = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)
      A1 = ZMASS/FLOAT(N)
*
      WRITE (6,15)  N, ZMASS, BE(3), CMR(4), RESC, STEPI, A1,
     &              NCRIT1, (JLIST(J),J=1,NCORR)
   15 FORMAT (/,3X,'ESCAPE   N =',I4,F7.3,F10.6,2F7.2,F7.3,F7.2,I5,10I4)
*
*       Enforce new NLIST in routine INTGRT.
      TLIST = TIME
*
   25 RETURN
*
   30 A2 = (X(1,JMIN) - RDENS(1))**2 + (X(2,JMIN) - RDENS(2))**2 +
     &                                 (X(3,JMIN) - RDENS(3))**2
*       See whether nearest body also satisfies the escape condition.
      IF (A2.LT.RESC2) THEN
          A3 = XDOT(1,JMIN)**2 + XDOT(2,JMIN)**2 + XDOT(3,JMIN)**2
          A4 = (XDOT(1,I) - XDOT(1,JMIN))**2 +
     &         (XDOT(2,I) - XDOT(2,JMIN))**2 +
     &         (XDOT(3,I) - XDOT(3,JMIN))**2
          A5 = (BODY(I) + BODY(JMIN))/SQRT(RJMIN2 + EPS2)
*       Check velocity of binary component in case of bound pair.
          A6 = 2.0*A5/A4
          IF (A6.LT.1.0.OR.A3.LT.2.0*ZMASS/SQRT(A2)) THEN
*       Retain escaper if dynamical effect on neighbour is significant.
              IF (A6.GT.0.01) GO TO 10
          END IF
      END IF
*
*       Form output diagnostics.
      RESC = SQRT(RI2)
      STEPI = STEP(I)
      NCORR = NCORR + 1
      JLIST(NCORR) = NAME(I)
*
*       Correct total energy & mass.
      BE(3) = BE(3) - BODY(I)*EI
      ZMASS = ZMASS - BODY(I)
*
*       Reduce particle number and update COMMON arrays.
      N = N - 1
      CALL REMOVE(I)
*
*       Check next particle (NB! Same location).
      GO TO 11
*
      END
