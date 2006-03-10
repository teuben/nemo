      SUBROUTINE HMDOT(IGHOST,IMERGE,M1,KW,MC,DMS,RNEW,ITERM)
*
*
*       Mass loss from inner hierarchical binary.
*       -----------------------------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8  M1,MC
      LOGICAL IQUIT
*
*
*       Define merger termination index and relevant KS/c.m. indices.
      ITERM = 0
      IPAIR = KSPAIR
      I1 = 2*IPAIR - 1
      ICM = N + IPAIR
*
*       Decide between #I1 & IGHOST (NB! choose I1 first if TEV equal).
      J = IGHOST
      IF (TEV(I1).LE.TEV(J)) J = I1
*
*       Include quitting conditions to avoid large changes at end-point.
      IQUIT = .FALSE.
*       Note that tidal dissipation needs full KS for rectification in MDOT.
      IF (KSTARM(IMERGE).LT.0.OR.BODY(ICM).EQ.0.0D0) IQUIT = .TRUE.
*
*       Quit for mis-identification, tidal evolution or double merger.
      IF (J.LE.0.OR.IQUIT) THEN
          ITERM = 1
          GO TO 50
      END IF
*
*       Update radius for #J and copy #I1 for MDOT & HCORR in case of ghost.
      RSTAR = RNEW
      IF (J.GT.I1) THEN
          RADIUS(J) = RNEW
          RSTAR = RADIUS(I1)
*       Return RNEW for copying back to RADIUS(I1) on ghost & DM/M < 0.01.
      END IF
*
*       Skip further modifications on small mass loss of same type.
      IF (ABS(DMS/M1).LT.0.01.AND.KW.EQ.KSTAR(J)) GO TO 50
*
*       Set two-body elements for outer orbit and inner semi-major axis.
      SEMI = -0.5*BODY(ICM)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(ICM)*SEMI)
      ECC1 = SQRT(ECC2)
      ZMB0 = CM(1,IMERGE) + CM(2,IMERGE)
      SEMI0 = -0.5*ZMB0/HM(IMERGE)
*
*       Form modified mass ratio and semi-major axes from M*A = const.
      DM = DMS/ZMBAR
      ZMB = ZMB0 - DM
      SEMI1 = SEMI*(ZMB0 + BODY(2*IPAIR))/(ZMB + BODY(2*IPAIR))
      SEMI2 = SEMI0*ZMB0/ZMB
*
*       Evaluate old separation, square regularized velocity, t'' & ECC.
      RI = 0.0D0
      V20 = 0.0
      TD2 = 0.0
      DO 10 K = 1,4
          RI = RI + UM(K,IMERGE)**2
          V20 = V20 + UMDOT(K,IMERGE)**2
          TD2 = TD2 + 2.0*UM(K,IMERGE)*UMDOT(K,IMERGE)
   10 CONTINUE
      ECC2 = (1.0 - RI/SEMI0)**2 + TD2**2/(SEMI0*ZMB0)
      ECC = SQRT(ECC2)
*
*       Determine inclination (use #I1 as first KS component).
      CALL HIMAX(I1,IMERGE,ECC,SEMI0,EMAX,EMIN,ZI,TG,EDAV)
*
*       Obtain stability parameters of the new configuration.
      PCRIT = stability(CM(1,IMERGE),CM(2,IMERGE),BODY(2*IPAIR),ECC,
     &                                            ECC1,ZI)*SEMI2
      PMIN = SEMI1*(1.0 - ECC1)
*
*       Update pericentre distance on successful stability test or exit.
      IF (PMIN.GT.PCRIT) THEN
          R0(IPAIR) = PCRIT
      ELSE
          ITERM = 1
          GO TO 50
      END IF
*
*       Check condition for sequential circularization.
      RP = SEMI0*(1.0 - ECC)
      IF (KZ(27).GT.0.AND.RP.LT.4.0*RADIUS(J).AND.
     &    KSTARM(IMERGE).EQ.0) THEN
          WRITE (6,15)  NAME(J), KSTAR(J), KSTARM(IMERGE), SQRT(ECC2),
     &                  RP, RADIUS(J)
   15     FORMAT (' HMDOT TERM    NAM K* E RP R*',I6,2I4,F8.4,1P,2E10.2)
          ITERM = 1
          GO TO 50
      END IF
*
*       Obtain energy change from M*A = const and H = -M/(2*A) (DCH 8/96).
      DH = DM/SEMI0*(1.0 - 0.5*DM/ZMB0)
      HM0 = HM(IMERGE)
      HM(IMERGE) = HM(IMERGE) + DH
*
*       Form KS coordinate & velocity scaling factors (general point is OK).
      SEMI2 = -0.5*ZMB/HM(IMERGE)
      C2 = SQRT(SEMI2/SEMI0)
      V2 = 0.5*(ZMB + HM(IMERGE)*RI*(SEMI2/SEMI0))
      C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy (H < 0: constant eccentricity).
      DO 20 K = 1,4
          UM(K,IMERGE) = C2*UM(K,IMERGE)
          UMDOT(K,IMERGE) = C1*UMDOT(K,IMERGE)
*       Note: no need to update XREL & VREL for RESET (c.m. error cancels).
   20 CONTINUE
*
*       Reduce mass of relevant component (ghost is original second member).
      ZMU0 = CM(1,IMERGE)*CM(2,IMERGE)/ZMB0
      KM = 1
      IF (J.GE.IFIRST) KM = 2
      CM(KM,IMERGE) = CM(KM,IMERGE) - DM
*
*       Include corrections to EMERGE & EMDOT (consistency; no net effect!).
      ZMU = CM(1,IMERGE)*CM(2,IMERGE)/ZMB
      DECORR = ZMU*HM(IMERGE) - ZMU0*HM0
      EMERGE = EMERGE + DECORR
      EMDOT = EMDOT - DECORR
*
*       Correct outer orbit for mass loss (use #I1 in case #J is a ghost).
      CALL HCORR(I1,DM,RSTAR)
*
*       Print some diagnostics on significant mass loss.
      IF (DMS.GT.0.005*M1) THEN
          WRITE (6,30)  NAME(J), KSTAR(J), IMERGE, KM, M1, DMS, PMIN,
     &                  PCRIT, SEMI2, SEMI1, DECORR
   30     FORMAT (' HMDOT    NAM K* IM KM M1 DM PM PC A A1 DE ',
     &                       I6,3I4,2F6.2,1P,4E10.2,0P,F10.6)
      END IF
*
   50 RETURN
*
      END
