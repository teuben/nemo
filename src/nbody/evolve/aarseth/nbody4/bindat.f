      SUBROUTINE BINDAT
*
*
*       Binary data bank.
*       -----------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/ CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &               HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &               NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8 XX(3,3),VV(3,3)
      REAL*8 EB(KMAX),ECC(KMAX),RCM(KMAX),ECM(KMAX),PB(KMAX),AS(30)
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
*
*
*       Form binding energy and central distance for each KS pair.
      ZMBIN = 0.d0
      DO 10 JPAIR = 1,NPAIRS
          J2 = 2*JPAIR
          J1 = J2 - 1
          ICM = N + JPAIR
          ZMBIN = ZMBIN + BODY(ICM)
          BODYCM = BODY(ICM)
*       Determine merger & ghost index for negative c.m. name (skip ghost).
          IF (NAME(ICM).LT.0.AND.BODY(ICM).GT.0.0) THEN
              CALL FINDJ(J1,J,IMERGE)
*       Employ actual masses and two-body distance for energy & eccentricity.
              BODYCM = CM(1,IMERGE) + CM(2,IMERGE)
              EB(JPAIR) = CM(1,IMERGE)*CM(2,IMERGE)*HM(IMERGE)/BODYCM
              SEMI = -0.5d0*BODYCM/HM(IMERGE)
              RJ = SQRT(XREL(1,IMERGE)**2 + XREL(2,IMERGE)**2 +
     &                                      XREL(3,IMERGE)**2)
              TD2 = 0.d0
              DO 1 K = 1,4
                  TD2 = TD2 + 2.d0*UM(K,IMERGE)*UMDOT(K,IMERGE)
    1         CONTINUE
              ECC2 = (1.d0 - RJ/SEMI)**2 + TD2**2/(BODYCM*SEMI)
*       Include separate diagnostics for the hierarchy (inner comps J1 & J).
              SEMI1 = -0.5d0*BODY(ICM)/H(JPAIR)
              ECC1 = (1.d0 - R(JPAIR)/SEMI1)**2 +
     &                          TDOT2(JPAIR)**2/(BODY(ICM)*SEMI1)
              E0 = SQRT(ECC2)
              E1 = SQRT(ECC1)
              Q = BODY(J2)/BODYCM
              XFAC = (1.d0 + Q)*(1.d0 + E1)/SQRT(ABS(1.d0 - E1))
              FE = 2.8d0
              PCR = FE*XFAC**0.4d0*SEMI
              PM = SEMI1*(1.d0 - E1)/PCR
              IF (J.LT.0) J = J1
              RM = SEMI*(1.d0 - E0)/MAX(RADIUS(J1),RADIUS(J))
              RM = MIN(RM,99.9d0)
              P0 = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
              P1 = DAYS*SEMI1*SQRT(ABS(SEMI1)/BODY(ICM))
              DO 2 K = 1,3
                  XX(K,1) = XREL(K,IMERGE)
                  XX(K,2) = 0.0
                  XX(K,3) = X(K,J2)
                  VV(K,1) = VREL(K,IMERGE)
                  VV(K,2) = 0.0
                  VV(K,3) = XDOT(K,J2)
    2         CONTINUE
              CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ALPHA)
              WRITE (84,3) TTOT, NAME(J1), NAME(J), KSTAR(J1), KSTAR(J),
     &                     KSTARM(IMERGE), E0, E1, PM, RM, ALPHA, P0,
     &                     P1, SEMI1
    3         FORMAT (' BINDAT:    T NM K* E0 E1 PM/PC PM0/R* IN P0 P1',
     &                             ' A1 ',
     &                             F7.1,2I5,2I3,I4,2F7.3,3F6.1,1P,3E9.1)
              CALL FLUSH(84)
          ELSE IF (BODY(J1).GT.0.0D0) THEN
*       Form binding energy and eccentricity for standard case.
              EB(JPAIR) = BODY(J1)*BODY(J2)*H(JPAIR)/
     &                                             (BODY(J1) + BODY(J2))
              SEMI = -0.5d0*BODY(ICM)/H(JPAIR)
              ECC2 = (1.d0 - R(JPAIR)/SEMI)**2 +
     &                           TDOT2(JPAIR)**2/(BODY(ICM)*SEMI)
          ELSE
              IM = 0
*       Search merger table to identify corresponding index of c.m. name.
              DO 5 K = 1,NMERGE
                  IF (NAMEG(K).EQ.NAME(ICM)) THEN
                      IM = K
                  END IF
    5         CONTINUE
              BODYJ1 = CM(3,IM)
              BODYJ2 = CM(4,IM)
              BODYCM = BODYJ1 + BODYJ2
              IF (BODYCM.EQ.0.0D0) THEN
                  WRITE (6,8)  JPAIR, IM, BODYCM
    8             FORMAT (' WARNING!    BINDAT    KS IM M ',2I4,1P,E9.1)
                  BODYCM = 1.0d-10
              END IF
              EB(JPAIR) = BODYJ1*BODYJ2*H(JPAIR)/BODYCM
              SEMI = -0.5d0*BODYCM/H(JPAIR)
              ECC2 = (1.d0 - SEMI/R(JPAIR))**2
          END IF
          ECC(JPAIR) = SQRT(ECC2)
          EB(JPAIR) = MAX(EB(JPAIR),-9.99999d0)
          PB(JPAIR) = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
          PB(JPAIR) = MIN(PB(JPAIR),99999.9d0)
*       Obtain binding energy (per unit mass) of c.m. motion.
          VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + 
     &                           XDOT(3,ICM)**2
          IF (BODY(ICM).EQ.0.0D0) VJ2 = 0.0
          ECM(JPAIR) = 0.5d0*VJ2 + PHI(ICM)
          RCM(JPAIR) = SQRT((X(1,ICM) - RDENS(1))**2 +
     &                      (X(2,ICM) - RDENS(2))**2 +
     &                      (X(3,ICM) - RDENS(3))**2)
          RCM(JPAIR) = MIN(RCM(JPAIR),99.9d0)
   10 CONTINUE
*
*       Copy relevant binary diagnostics to single precision.
      AS(1) = TTOT
      AS(2) = RSCALE
      AS(3) = RTIDE
      AS(4) = RC
      AS(5) = TPHYS
      AS(6) = -1.5d0*(TIDAL(1)*ZMASS**2)**0.3333
      AS(7) = 0.d0
      DO 20 K = 1,10
          AS(K+7) = E(K)
   20 CONTINUE
      AS(18) = SBCOLL
      AS(19) = BBCOLL
      AS(20) = ZKIN
      AS(21) = POT
      AS(22) = EBIN0
      AS(23) = EBIN
      AS(24) = ESUB
      AS(25) = EMERGE
      AS(26) = BE(3)
      AS(27) = ZMASS
      AS(28) = ZMBIN
      AS(29) = CHCOLL
      AS(30) = ECOLL
*
*       Write formatted data bank on unit 9.
      IF (FIRST) THEN
          OPEN (UNIT=9,STATUS='NEW',FORM='FORMATTED',FILE='OUT9')
          FIRST = .FALSE.
      END IF
*
      WRITE (9,30)  NPAIRS, MODEL, NRUN, N, NC, NMERGE, (AS(K),K=1,7)
   30 FORMAT (3I4,I6,2I4,2X,F7.1,2F7.2,F7.3,3F9.4)
      WRITE (9,35)  (AS(K),K=8,17)
   35 FORMAT (10F11.6)
      WRITE (9,40)  (AS(K),K=18,30)
   40 FORMAT (13F10.5)
*
      DO 50 JPAIR = 1,NPAIRS
          J1 = 2*JPAIR - 1
          J2 = 2*JPAIR
          KCM = KSTAR(N+JPAIR)
          IF (NAME(N+JPAIR).LT.0) THEN
              KCM = -10
              IF (NAME(N+JPAIR).LT.-2*NZERO) KCM = -20
          END IF
          WRITE (9,45)  EB(JPAIR), ECC(JPAIR), ECM(JPAIR), RCM(JPAIR),
     &                  BODY(J1)*ZMBAR, BODY(J2)*ZMBAR, PB(JPAIR),
     &                  NAME(J1), NAME(J2), KSTAR(J1), KSTAR(J2), KCM
   45     FORMAT (F8.5,F7.3,F7.2,F6.2,2F5.1,F8.1,2I6,3I4)
   50 CONTINUE
      CALL FLUSH(9)
*
      RETURN
*
      END
