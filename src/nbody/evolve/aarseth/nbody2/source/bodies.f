      SUBROUTINE BODIES
*
*
*       Output of single bodies or binaries.
*       ------------------------------------
*
      INCLUDE 'common2.h'
      REAL*4  A(3)
*
*
*       Check option for printing single bodies.
      IF (KZ(9).EQ.0) GO TO 20
      K = KZ(9)
      IBODY = MIN(5**K,N)
*
      DO 10 I = 1,IBODY
          FIRR = SQRT(FI(1,I)**2 + FI(2,I)**2 + FI(3,I)**2)
          FREG = SQRT(FR(1,I)**2 + FR(2,I)**2 + FR(3,I)**2)
          EI = 0.5*(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
          DO 4 J = 1,N
              IF (J.EQ.I) GO TO 4
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              EI = EI - BODY(J)/SQRT(RIJ2 + EPS2)
    4     CONTINUE
          DO 5 K = 1,3
             A(K) = X(K,I) - RDENS(K)
    5     CONTINUE
          RI = SQRT(A(1)**2 + A(2)**2 + A(3)**2)
          WRITE (6,6)  I, NAME(I), BODY(I), STEP(I), STEPR(I), EI, RI,
     &                 LIST(1,I), (A(K),K=1,3), (XDOT(K,I),K=1,3),
     &                 FIRR, FREG, RS(I)
    6     FORMAT (I6,I5,2F8.4,F7.3,F7.1,F7.2,I6,3X,3F7.2,3X,3F6.2,3X,
     &                                                   F7.1,F6.1,F7.2)
   10 CONTINUE
*
*       Optional search for binaries (frequency NFIX with KZ(6) = 2).
   20 IF (KZ(6).EQ.0) GO TO 50
      IF (KZ(6).EQ.2.AND.NPRINT.NE.1) GO TO 50
      SMAX = 0.01*TCR
*
      DO 40 I = 1,N
          IF (STEP(I).GT.SMAX) GO TO 40
          JMIN = 0
          RJMIN2 = RSCALE**2
          NNB = LIST(1,I)
          DO 30 L = 1,NNB
              J = LIST(L+1,I)
              IF (STEP(J).GT.SMAX) GO TO 30
              A1 = X(1,I) - X(1,J)
              A2 = X(2,I) - X(2,J)
              A3 = X(3,I) - X(3,J)
              RIJ2 = A1**2 + A2**2 + A3**2 + EPS2
              IF (RIJ2.LT.RJMIN2) THEN
                  RJMIN2 = RIJ2
                  JMIN = J
              END IF
   30     CONTINUE
          IF (JMIN.LT.I) GO TO 40
          RIJMIN = SQRT(RJMIN2)
          VR2 = (XDOT(1,I) - XDOT(1,JMIN))**2 +
     &          (XDOT(2,I) - XDOT(2,JMIN))**2 +
     &          (XDOT(3,I) - XDOT(3,JMIN))**2
          EREL = 0.5*VR2 - (BODY(I) + BODY(JMIN))/RIJMIN
          IF (EREL.GT.0.0) GO TO 40
          SEMI = -0.5*(BODY(I) + BODY(JMIN))/EREL
*       Only print significant binaries.
          IF (SEMI.GT.0.1*RSCALE) GO TO 40
          ZN = SQRT((BODY(I) + BODY(JMIN))/SEMI**3)
          RDOT = (X(1,I) - X(1,JMIN))*(XDOT(1,I) - XDOT(1,JMIN)) +
     &           (X(2,I) - X(2,JMIN))*(XDOT(2,I) - XDOT(2,JMIN)) +
     &           (X(3,I) - X(3,JMIN))*(XDOT(3,I) - XDOT(3,JMIN))
          ECC2 = (1.0 - RIJMIN/SEMI)**2 +
     &                             RDOT**2/(SEMI*(BODY(I) + BODY(JMIN)))
          ECC = SQRT(ECC2)
          RI = SQRT((X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 + 
     &                                       (X(3,I) - RDENS(3))**2)
          WRITE (6,35)  NAME(I), NAME(JMIN), BODY(I), BODY(JMIN), EREL,
     &                  SEMI, ZN, RIJMIN, RI, ECC, LIST(1,I)
   35     FORMAT (3X,'BINARY ',2I5,2F8.4,F7.1,F9.5,F8.1,F9.5,F7.2,F8.3,
     &                                                               I5)
   40 CONTINUE
*
   50 RETURN
*
      END
