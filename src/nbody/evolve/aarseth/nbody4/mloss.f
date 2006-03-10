      SUBROUTINE MLOSS
*
*
*       Mass loss from evolving stars.
*       ------------------------------
*
*       Original scheme for Elena Terlevich thesis 1983.
*       ------------------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8  A(9)
*
*
*       Find the heaviest body (exclude any merged binary components).
   99 IMAX = 0
    1 BMAX = 0.0D0
      KS = 0
      DO 2 I = 1,N
          IF (BODY(I).LT.BMAX.OR.I.LE.IMAX) GO TO 2
          BMAX = BODY(I)
          IMAX = I
    2 CONTINUE
*
      ZMSTAR = ZMBAR*BMAX
*       Obtain total evolution time for ZMSTAR in solar masses.
      TMS = (2.55D+03 + 6.69D2*ZMSTAR**2.5D0 + ZMSTAR**4.5D0)/
     &      (3.27D-02*ZMSTAR**1.5D0 + 3.46D-01*ZMSTAR**4.5D0)
      TG = 0.15*TMS
      THE = TMS*1.37*ZMSTAR**(-0.881D0)
      TEV1 = TMS + TG + THE
*       Scale the evolution time to model units.
      TMDOT = TEV1/TSCALE
*
      IF (TMDOT.GT.TIME) THEN
*       Update the maximum mass.
          BODY1 = BMAX
*       Set phase indicator = -1 for new NLIST in routine INTGRT.
          IPHASE = -1
          ISEND = -1
          RETURN
      END IF
*
      IBODY = IMAX
*       Include special treatment for regularized components.
      IF (IMAX.GT.2*NPAIRS) GO TO 4
      IPAIR = KVEC(IMAX)
      IBODY = N + IPAIR
*       Obtain current coordinates and velocities.
      CALL RESOLV(IPAIR,2)
*
*       Check merger termination.
      IF (NAME(IBODY).LT.0) THEN
          KSPAIR = IPAIR
          IPHASE = 7
          CALL RESET
          GO TO 1
      END IF
*
*       Set variable white dwarf mass (Iben & Renzini, Ann. Rev. 21, 298).
    4 ZMSTAR = 0.38 + 0.15*ZMSTAR
*       Assume neutron star instead if mass > 6 solar masses.
      IF (ZMBAR*BMAX.GT.6.0) ZMSTAR = 1.5
      DM = BODY(IMAX) - ZMSTAR/ZMBAR
      ZMASS = ZMASS - DM
      ZMDOT = ZMDOT + DM*ZMBAR
      NMDOT = NMDOT + 1
      IF (TIDAL(1).NE.0.0D0.AND.KZ(14).LT.3)THEN
          RTIDE = (ZMASS/TIDAL(1))**0.3333
      END IF
*
*       Save old velocity for FDOT correction.
      VI2 = 0.0
      DO 5 K = 1,3
          A(K+6) = XDOT(K,IMAX)
          VI2 = VI2 + XDOT(K,IMAX)**2
    5 CONTINUE
*
*       Set first part of energy correction.
      DE = 0.5*DM*(VI2 - TIDAL(1)*X(1,IMAX)**2 - TIDAL(3)*X(3,IMAX)**2)
      A0 = 1.0
*
*       Terminate any KS pair (bug fix May 2005).
      IF (IMAX.LT.IFIRST) THEN
          KSPAIR = IPAIR
          IMAX = IMAX + 2*(NPAIRS - KSPAIR)
          CALL RESOLV(KSPAIR,3)
          RR = R(IPAIR)
          IPHASE = 2
          JCOMP = 0
          CALL KSTERM
          KS = 1
      END IF
*       Update the mass after possible KS transformations.
      BODY(IMAX) = ZMSTAR/ZMBAR
      IF (ZMBAR*BMAX.LE.6.0) GO TO 8
*
*       Assign a high velocity to neutron star (at least 4*rms velocity).
      A0 = 4.0*SQRT(0.5D0*ZMASS/(RSCALE*VI2))
      IF (KS.EQ.0) GO TO 6
      A2 = 2.0*(BODY(IFIRST) + BODY(IFIRST+1))/RR
*       Add escape velocity from the regularized pair.
      A0 = SQRT(A0**2 + A2/VI2)
    6 DO 7 K = 1,3
          XDOT(K,IMAX) = A0*XDOT(K,IMAX)
          X0DOT(K,IMAX) = XDOT(K,IMAX)
          X0(K,IMAX) = X(K,IMAX)
    7 CONTINUE
*
*       Correct total energy, forces & first derivatives.
    8 POTI = 0.0D0
      DO 20 J = IFIRST,NTOT
          IF (J.EQ.IMAX) GO TO 20
          RIJDOT = 0.0D0
          RDVDOT = 0.0D0
*
          DO 10 K = 1,3
              A(K) = X(K,IMAX) - X(K,J)
              A(K+3) = A(K+6) - XDOT(K,J)
              RIJDOT = RIJDOT + A(K)*A(K+3)
              RDVDOT = RDVDOT + A(K)*(XDOT(K,IMAX) - A(K+6))
   10     CONTINUE
*
          RIJ2 = A(1)**2 + A(2)**2 + A(3)**2
          POTI = POTI + BODY(J)/SQRT(RIJ2)
          IF (J.EQ.IBODY) GO TO 20
          A3 = 1.0/(RIJ2*SQRT(RIJ2))
          A4 = BODY(IMAX)*A3
          A5 = DM*A3
          A6 = 3.0*RIJDOT/RIJ2
          A7 = 3.0*RDVDOT/RIJ2
*
          DO 11 K = 1,3
              A(K+3) = (A(K+3) - A(K)*A6)*A5
              IF (A0.GT.1.0) THEN
*       Include FDOT corrections due to increased velocity.
                  A(K+3) = A(K+3) + (XDOT(K,IMAX) - A(K+6))*A4
                  A(K+3) = A(K+3) - A7*A(K)*A4
              END IF
   11     CONTINUE
*
          DO 18 K = 1,3
              F(K,J) = F(K,J) - 0.5*A(K)*A5
              FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
   18     CONTINUE
   20 CONTINUE
*
*       Initialize the mass-losing body (except any c.m.).
      IF (IBODY.LE.N.AND.KS.EQ.0) THEN
          DO 25 K = 1,3
              X0DOT(K,IBODY) = XDOT(K,IBODY)
   25     CONTINUE
          CALL GPSEND
          T0(IBODY) = TIME
          IPHASE = -3
          CALL FPOLYI(IBODY)
      END IF
*
*       Correct total energy for mass loss effect.
      DE = DE - DM*POTI
      EMDOT = EMDOT + DE
      RI = SQRT((X(1,IMAX) - RDENS(1))**2 + (X(2,IMAX) - RDENS(2))**2 +
     &                                      (X(3,IMAX) - RDENS(3))**2)
      WRITE (6,30)  NAME(IMAX), BMAX*ZMBAR, DM*ZMBAR, ZMDOT, DE, BE(3),
     &              TMDOT, RI
   30 FORMAT (/,'   MASS LOSS   NAME =',I5,'  M* =',F5.1,' DMS =',F7.3,
     &                       '  ZM* =',F6.1,'  DE =',F9.5,'  E =',F10.6,
     &                       '  TMDOT =',F6.1,'  RI =',F5.2)
*
      IF (ZMBAR*BMAX.LE.6.0) GO TO 50
*       Include correction to total energy from the increased velocity.
      DE = 0.5*BODY(IMAX)*VI2*(A0**2 - 1.0)
      ECDOT = ECDOT - DE
      WRITE (6,45)  IMAX, KS, DE, TEV1, A0, SQRT(VI2), BE(3)
   45 FORMAT ('   RECOIL   I =',I5,'  KS =',I3,'  DE =',F9.5,
     &        '  T* =',F5.1,'  V/V0 =',F6.2,'  V0 =',F5.2,'  E =',F10.6)
*
   50 IF (KS.GT.0) THEN
          ICOMP = IFIRST
          JCOMP = IFIRST + 1
          IPHASE = 1
          CALL GPSEND
          CALL KSREG
      END IF
*
*       Set next mass loss time & current maximum mass before returning.
      GO TO 99
*
      END
