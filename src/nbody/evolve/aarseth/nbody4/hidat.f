      SUBROUTINE HIDAT
*
*
*       Hierarchical data bank.
*       -----------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  XX(3,3),VV(3,3),M1,M2,M3
      REAL*8  EB(KMAX),ECM(KMAX)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*
*       Write formatted data bank on unit 87.
      IF (FIRST) THEN
          OPEN (UNIT=87,STATUS='NEW',FORM='FORMATTED',FILE='HIDAT')
          FIRST = .FALSE.
          WRITE (87,1)
    1     FORMAT (/,'  NAM1  NAM2  NAM3  K*       M1   M2   M3    RI',
     &            '    EMAX   E0    E1      P0    P1')
      END IF
*
      IF (NMERGE.GT.0) THEN
          MULT = 0
*       Count higher-order systems.
          DO 2 I = N+1,NTOT
              IF (NAME(I).LT.-2*NZERO) MULT = MULT + 1
    2     CONTINUE
          WRITE (87,3)  NPAIRS, NRUN, N, NC, NMERGE, MULT, NEWHI, TTOT
    3     FORMAT (/,I6,I4,I6,3I4,I6,F9.1)
      END IF
*
*       Form output parameters for each merged KS pair (TRIPLE & QUAD).
      IPAIR = 0
      ISTAB = 0
      DO 30 JPAIR = 1,NPAIRS
          ICM = N + JPAIR
          IF (NAME(ICM).GT.0.OR.NAME(ICM).LT.-2*NZERO) GO TO 30
          J2 = 2*JPAIR
          J1 = J2 - 1
*       Determine merger & ghost index (delay ghost c.m.).
          CALL FINDJ(J1,J,IM)
          SEMI1 = -0.5*BODY(ICM)/H(JPAIR)
          ECC1 = (1.0 - R(JPAIR)/SEMI1)**2 +
     &                                 TDOT2(JPAIR)**2/(BODY(ICM)*SEMI1)
          IF (BODY(ICM).GT.0.0) THEN
              IPAIR = IPAIR + 1
*       Employ actual masses and two-body distance for energy & eccentricity.
              BODYCM = CM(1,IM) + CM(2,IM)
              EB(IPAIR) = CM(1,IM)*CM(2,IM)*HM(IM)/BODYCM
              SEMI = -0.5*BODYCM/HM(IM)
              RJ = SQRT(XREL(1,IM)**2 + XREL(2,IM)**2 + XREL(3,IM)**2)
              TD2 = 0.0
              DO 5 K = 1,4
                  TD2 = TD2 + 2.0*UM(K,IM)*UMDOT(K,IM)
    5         CONTINUE
              ECC2 = (1.0 - RJ/SEMI)**2 + TD2**2/(BODYCM*SEMI)
*       Include separate diagnostics for the hierarchy (inner comps J1 & J).
              M1 = CM(1,IM)*SMU
              M2 = CM(2,IM)*SMU
              M3 = BODY(J2)*SMU
*       Retain the next part just in case.
              E0 = SQRT(ECC2)
              E1 = SQRT(ECC1)
*       Determine EMAX and other relevant parameters (if needed).
              CALL HIMAX(J1,IM,E0,SEMI,EMAX,EMIN,ZI,TG,EDAV)
              Q = BODY(J2)/BODYCM
              XFAC = (1.0 + Q)*(1.0 + E1)/SQRT(ABS(1.0 - E1))
              FE = 2.8
              PCR = FE*XFAC**0.4*SEMI
              PMIN = SEMI1*(1.0 - E1)
              IF (J.LT.0) J = J1
              RM = SEMI*(1.0 - E0)/MAX(RADIUS(J1),RADIUS(J),1.0D-20)
              RM = MIN(RM,99.9D0)
*       Obtain inclination between inner relative motion and outer orbit.
              DO 10 K = 1,3
                  XX(K,1) = XREL(K,IM)
                  XX(K,2) = 0.0
                  XX(K,3) = X(K,J2)
                  VV(K,1) = VREL(K,IM)
                  VV(K,2) = 0.0
                  VV(K,3) = XDOT(K,J2)
   10         CONTINUE
              CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ALPHA)
              NAM2 = NAME(J)
              KSTAR3 = KSTAR(J2)
*       Perform stability check for new inclination (skip quadruples).
              YF = 1.0 - 0.3*ALPHA/180.0
              IF (PMIN*(1.0 - GAMMA(JPAIR)).LT.YF*PCR.AND.
     &            E1.LT.0.96) THEN
                  ISTAB = JPAIR
                  WRITE (6,15)  NAME(J1), NAME(J), ALPHA, E1, PMIN,
     &                          YF*PCR, GAMMA(JPAIR)
   15             FORMAT (' MERGE UNSTAB    NAM INC E1 PM YF*PCR G ',
     &                                      2I6,F7.1,F7.3,1P,3E10.2)
                  TMDIS(IM) = TIME
              END IF
          ELSE IF (BODY(J1).EQ.0.0D0) THEN
              IPAIR = IPAIR + 1
              BODYJ1 = CM(3,IM)
              BODYJ2 = CM(4,IM)
              BODYCM = BODYJ1 + BODYJ2
              BODYCM = MAX(BODYCM,1.0D-10)
              M1 = BODYJ1*SMU
              M2 = BODYJ2*SMU
              M3 = (BODY(ICM) - BODYCM)*SMU
              EB(IPAIR) = BODYJ1*BODYJ2*H(JPAIR)/BODYCM
              SEMI = -0.5*BODYCM/H(JPAIR)
              ECC2 = (1.0 - SEMI/R(JPAIR))**2
              ALPHA = 0.0
              EMAX = 0.0
              NAM2 = NAME(J)
              KSTAR3 = KSTARM(IM)
          END IF
*
          E0 = SQRT(ECC2)
          E1 = SQRT(ECC1)
          P0 = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
          P1 = DAYS*SEMI1*SQRT(ABS(SEMI1)/BODY(ICM))
          P0 = MIN(P0,9999.0D0)
*       Obtain binding energy (per unit mass) of c.m. motion.
          VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2
          IF (BODY(ICM).EQ.0.0D0) VJ2 = 0.0
          ECM(IPAIR) = 0.5*VJ2 + PHI(ICM)
          RI = SQRT((X(1,ICM) - RDENS(1))**2 + (X(2,ICM) - RDENS(2))**2
     &                                       + (X(3,ICM) - RDENS(3))**2)
          WRITE (87,25)  NAME(J1), NAM2, NAME(J2), KSTAR(J1), KSTAR(J),
     &                   KSTAR3, M1, M2, M3, RI, EMAX, E0, E1,
     &                   P0, P1
   25     FORMAT (3I6,3I3,3F5.1,F7.2,F7.3,2F6.2,F8.1,1P,E9.1)
   30 CONTINUE
      CALL FLUSH(87)
*
*       Check merger termination of unstable system.
*     IF (ISTAB.GT.0) THEN
*         KSPAIR = ISTAB
*         IPHASE = 7
*         CALL RESET
*     END IF
*
      RETURN
*
      END
