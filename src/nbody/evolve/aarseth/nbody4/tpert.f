      SUBROUTINE TPERT(IPAIR,GA,DT)
*
*
*       Perturbation time scale.
*       ------------------------
*
      INCLUDE 'common4.h'
*
*
*       Set c.m. index and initialize scalars.
      I = N + IPAIR
      FMAX = 0.0
      DTIN = 1.0E+20
      JCLOSE = 0
      NTPERT = NTPERT + 1
*
*       Define neighbour distance from density fitting expression.
      RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                               (X(3,I) - RDENS(3))**2
      RS2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
*
*       Ensure minimum search distance of four times the unperturbed size.
      IF (LIST(1,2*IPAIR-1).EQ.0) THEN
          FACM = (2.0*BODYM/BODY(I))**0.66667
          RS2 = MAX(RS2,16.0*FACM*CMSEP2*R(IPAIR)**2)
      END IF
*
*       Include modification by large c.m. step (typical rms velocity 1.0).
      RS2 = MAX(RS2,4.0*STEP(I)**2)
*
*       Adopt cutoff from n(r) = N/2*(r/RSCALE)**2 to limit GRAPE list < 10.
      IF (RS2.GT.0.01*RSCALE**2) THEN
          RS2 = MIN(RS2,20.0D0*RSCALE**2/FLOAT(N))
      END IF
*
*       Use option #39 = 2 for fast neighbour list (one member only).
      KZ39 = KZ(39)
      IF (LIST(1,2*IPAIR-1).EQ.0.AND.KZ(39).EQ.0) THEN
          KZ(39) = 2
      END IF
*
*       Obtain neighbour list inside modified interparticle distance.
      CALL NBLIST(I,RS2,NNB)
      KZ(39) = KZ39
*
      IF (NNB.EQ.0) THEN
          DT = MIN(2.0*STEP(I),1.0D0)
          GO TO 20
      END IF
*
*       Find the most likely perturbers (first approach & maximum force).
      DO 10 L = 2,NNB+1
          J = ILIST(L)
          RIJ2 = 0.0
          RDOT = 0.0
*
          DO 6 K = 1,3
              XREL = X(K,J) - X(K,I)
              VREL = XDOT(K,J) - XDOT(K,I)
              RIJ2 = RIJ2 + XREL**2
              RDOT = RDOT + XREL*VREL
    6     CONTINUE
*
          VR = RDOT/RIJ2
          IF (VR.LT.DTIN) THEN
              DTIN = VR
*       Note DTIN is inverse travel time to include case of no RDOT < 0.
              RCRIT2 = RIJ2
              JCRIT = J
          END IF
          FIJ = (BODY(I) + BODY(J))/RIJ2
          IF (FIJ.GT.FMAX) THEN
              FMAX = FIJ
              RJMIN2 = RIJ2
              JCLOSE = J
          END IF
   10 CONTINUE
*
*       Form radial velocity of body with shortest approach time (if any).
      RCRIT = SQRT(RCRIT2)
      RDOT = RCRIT*ABS(DTIN)
      A1 = 2.0/(BODY(I)*GA)
      SEMI = -0.5*BODY(I)/H(IPAIR)
*
*       Distinguish between actual and estimated apocentre.
      IF (LIST(1,2*IPAIR-1).EQ.0.AND.R(IPAIR).GT.SEMI) THEN
          RI = R(IPAIR)
      ELSE
          RI = 2.0*SEMI
      END IF
*
*       Estimate time interval to reach tidal perturbation of GA.
      IF (ABS(RDOT).LT.1.0E-10) RDOT = 1.0E-10
      DT = (RCRIT - RI*(BODY(JCRIT)*A1)**0.3333)/RDOT
*
*       Include a second test using typical heavy body and high velocity.
      RDOT2 = MAX(RDOT,4.0D0)
      ZM = MIN(5.0*BODYM,BODY1)
      DT2 = (RCRIT - RI*(ZM*A1)**0.3333)/RDOT2
*
*       Compare the travel time based on acceleration only.
      DTMAX = SQRT(2.0D0*ABS(DT)*RDOT*RCRIT2/(BODY(I) + BODY(JCRIT)))
      DT = MIN(DT,DTMAX,DT2)
*
*       Skip dominant force test if there is only one critical body.
      IF (JCRIT.NE.JCLOSE) THEN
*       Form the return time of the dominant body and choose the minimum.
          DR = SQRT(RJMIN2) - RI*(BODY(JCLOSE)*A1)**0.3333
          DTMAX = SQRT(2.0*ABS(DR)/FMAX)
          DT = MIN(DT,DTMAX)
      END IF
*
*       Apply safety test in case background force dominates c.m. motion.
      DT = MIN(DT,4.0D0*STEP(I))
      DT = MAX(DT,0.0D0)
*
   20 RETURN
*
      END
