      SUBROUTINE NBLIST(I,RS0)
*
*
*       Neighbour list & radius.
*       ------------------------
*
      INCLUDE 'common2.h'
*
*
*       Form square neighbour radius.
      RS2 = RS0**2
*
      IF (KZ(5).GE.1) THEN
*       Modify initial guess for Plummer model.
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          RS2 = RS2*(1.0 + RI2)
      END IF
*
      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
      IF (VI2.GT.0.0) THEN
          DTR = 0.1*RS0/SQRT(VI2)
          VRFAC = -0.1*RS2/DTR
*       Estimated radial velocity factor for outer sphere.
      ELSE
          VRFAC = 0.0
      END IF
*
*       Search all particles.
    1 RCRIT2 = 1.59*RS2
      NNB = 1
      DO 10 J = 1,N
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          IF (RIJ2.GT.RCRIT2) GO TO 10
          IF (RIJ2.GT.RS2) THEN
              RIJDOT = 0.0
              DO 4 K = 1,3
                  RIJDOT = RIJDOT +
     &                         (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
    4         CONTINUE
*       Accept member if maximum penetration factor is significant.
              IF (RIJDOT.GE.VRFAC) GO TO 10
          ELSE
              IF (J.EQ.I) GO TO 10
          END IF
          NNB = NNB + 1
          ILIST(NNB) = J
   10 CONTINUE
*
      IF (NNB.EQ.1) THEN
*       Double the neighbour sphere and try again.
          RS2 = 1.59*RS2
          NBVOID = NBVOID + 1
          GO TO 1
      ELSE IF (NNB.GT.NNBMAX) THEN
*       Reduce neighbour sphere but avoid possible looping.
          A1 = 0.75*FLOAT(NNBMAX)/FLOAT(NNB)
          IF (ABS(A1 - 0.5).LT.0.05) A1 = 1.2*A1
          RS2 = RS2*A1**0.66667
          NBFULL = NBFULL + 1
          GO TO 1
      END IF
*
*       Set neighbour radius and copy list membership.
      RS(I) = SQRT(RS2)
      LIST(1,I) = NNB - 1
      DO 20 L = 2,NNB
          LIST(L,I) = ILIST(L)
   20 CONTINUE
*
      RETURN
*
      END
