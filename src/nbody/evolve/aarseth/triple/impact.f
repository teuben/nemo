      SUBROUTINE IMPACT(PC)
*
*
*       Impact parameter.
*       -----------------
*
      IMPLICIT  REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  Q(8),P(8),R,R1,R2,ENERGY,M(3),X(3,3),XDOT(3,3),
     &               RCOLL,ERROR,C11,C12,C19,C20,C24,C25,NSTEPS,NAME(3)
      COMMON/CLOSE/  RIJ(3,3),ICALL
      COMMON/ANGLES/  TWOPI,ALPHA,DEGREE
      REAL*8  VI2(3)
*
*
*       Transform to physical variables.
      CALL TRANSF(3)
*       Identify index of second binary component & escaper.
      IMIN = 1
      IF (R2.LT.R1) IMIN = 2
      I = 3 - IMIN
*
      RIDOT = X(1,I)*XDOT(1,I) + X(2,I)*XDOT(2,I) + X(3,I)*XDOT(3,I)
*
*       Set distance & radial velocity of body #I with respect to binary.
      MB = M(3) + M(IMIN)
      FAC = (MB + M(I))/MB
      RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
      RIDOT = FAC*RIDOT/RI
      RI = FAC*RI
*
      DO 10 I = 1,3
          VI2(I) = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
   10 CONTINUE
*
*       Define probable escaper from maximum velocity.
      VMAX = 0.0
      DO 15 I = 1,3
          IF (VI2(I).GT.VMAX) THEN
              VMAX = VI2(I)
              IESC = I
          END IF
   15 CONTINUE
*
*       Determine particle IMIN closest to IESC.
      RMIN = 1000.0
      DO 20 I = 1,3
          IF (I.NE.IESC) THEN
              RIJ2 = (X(1,I) - X(1,IESC))**2 + (X(2,I) - X(2,IESC))**2
              IF (RIJ2.LT.RMIN) THEN
                  RMIN = RIJ2
                  IMIN = I
              END IF
          END IF
   20 CONTINUE
*
*       Determine the third particle IMAX.
      DO 25 I = 1,3
          IF (I.NE.IMIN.AND.I.NE.IESC) THEN
              IMAX = I
          END IF
   25 CONTINUE
*
*       Find sign of radial velocity of IESC & IMAX.
      RDOT = 0.0
      DO 30 K = 1,3
          RDOT = RDOT + (X(K,IESC) - X(K,IMAX))*
     &                                     (XDOT(K,IESC) - XDOT(K,IMAX))
   30 CONTINUE
*
*       Determine velocity angles of IESC & IMAX.
      AESC = ATAN2(XDOT(2,IESC),XDOT(1,IESC))
      AMAX = ATAN2(XDOT(2,IMAX),XDOT(1,IMAX))
      IF (AESC.LT.0.0) AESC = AESC + TWOPI
      IF (AMAX.LT.0.0) AMAX = AMAX + TWOPI
      AESC = DEGREE*AESC
      AMAX = DEGREE*AMAX
*
*       Switch IESC & IMIN if IESC & IMAX are moving in same direction.
      IF (RDOT.GT.0.0) THEN
          K = IESC
          IESC = IMIN
          IMIN = K
      END IF
*
*     WRITE (6,50) NAME(IESC), NAME(IMIN), NAME(IMAX), RDOT, AESC, AMAX
*  50 FORMAT ('  NAME(IESC;IMIN;IMAX) RDOT AESC AMAX ',3I4,F8.2,2F7.1)
      CALL FLUSH(6)
*
      RETURN
*
      END
