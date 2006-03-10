      SUBROUTINE ORBIT(I,J,SEMI,ECC,GI)
*
*
*       Two-body elements.
*       -----------------
*
      INCLUDE 'common4.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Find the dominant neighbour from LSHORT if body #J not specified.
      GI = 0.0
      NNB = LSHORT(1)
      IF (J.LE.0.AND.NNB.GT.0) THEN
          FMAX = 1.0D-10
          JM = LSHORT(2)
          IF (JM.EQ.I) THEN
              JM = LSHORT(3)
              IF (JM.EQ.I) THEN
                  J = N + 1
                  GO TO 100
              END IF
          END IF
          DO 2 L = 1,NNB
              JJ = LSHORT(L+1)
              JLIST(L) = JJ
              IF (JJ.EQ.I) GO TO 2
              RIJ2 = 0.0
              DO 1 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,JJ))**2
    1         CONTINUE
              FIJ = (BODY(I) + BODY(JJ))/RIJ2
*       Exclude any c.m. bodies from dominant motion.
              IF (FIJ.GT.FMAX.AND.JJ.LE.N) THEN
                  FMAX = FIJ
                  JM = JJ
              END IF
    2     CONTINUE
*
*       Obtain the relative perturbation for decision-making.
          CALL FPERT(I,JM,NNB,PERT)
          GI = PERT/FMAX
          J = JM
      END IF
      IF(J.EQ.0.OR.J.GT.N) GO TO 100
*
*       Skip any subsystem c.m. body.
      DO 3 ISUB = 1,NSUB
          NAMEI = NAMES(1,ISUB)
          IF (NAMEI.EQ.NAME(J)) THEN
              J = 0
              GO TO 100
          END IF
    3 CONTINUE
*
*       Determine the semi-major axis and eccentricity.
      RIJ2 = 0.0
      VIJ2 = 0.0
      RDOT = 0.0
      DO 5 K = 1,3
          RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
          VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
          RDOT = RDOT + (X(K,I) - X(K,J))*(XDOT(K,I) - XDOT(K,J))
    5 CONTINUE
      RIJ = SQRT(RIJ2)
      A1 = 2.0/RIJ - VIJ2/(BODY(I) + BODY(J))
      SEMI = 1.0/A1
      ECC2 = (1.0 - RIJ/SEMI)**2 + RDOT**2/(SEMI*(BODY(I) + BODY(J)))
      ECC = SQRT(ECC2)
*
*       Skip full neighbour list test on large apocentre or perturbation.
      IF (ECC.GT.1.0.OR.SEMI*(1.0 + ECC).GT.5.0*RMIN.OR.
     &    GI.GT.0.25) THEN
          GI = 1.0
          GO TO 100
      END IF
*
*       Obtain neighbour list (not all perturbers are in LSHORT).
      H2 = 25.0*CMSEP2*RMIN2
   10 CALL NBLIST(I,H2,NNB)
      IF ((NNB.LT.10.AND.N.GT.100).OR.
     &    (NNB.LT.1 + SQRT(FLOAT(N)/100))) THEN
          H2 = 2.0*H2
          GO TO 10
      END IF
*
*       Copy neighbour list and evaluate actual perturbing force.
      FMAX = 0.0
      DO 20 L = 1,NNB
          J = ILIST(L+1)
          JLIST(L) = J
          RIJ2 = 0.0
          DO 15 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   15     CONTINUE
          FIJ = (BODY(I) + BODY(J))/RIJ2
*       Exclude any c.m. bodies from dominant motion.
          IF (FIJ.GT.FMAX.AND.J.LE.N) THEN
              FMAX = FIJ
              JM = J
          END IF
   20 CONTINUE
      CALL FPERT(I,JM,NNB,PERT)
      GI = PERT/FMAX
*
  100 RETURN
*
      END
