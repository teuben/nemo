      SUBROUTINE OUTPUT
*
*
*       Output and data save.
*       ---------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/ECHAIN/  ECH
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
      REAL*8  X1(3,4),V1(3,4),UI(4),VI(4),XREL2(3),VREL2(3)
      REAL*4 XS(3,NMAX),VS(3,NMAX),BODYS(NMAX),AS(20)
      REAL*4 RHO1(NMAX),PHI1(NMAX)
      REAL*4  XJ(3,6),VJ(3,6),BODYJ(6)
      REAL*8 R2,RHO
      COMMON/WORK1/ R2(NMAX),RHO(NMAX)
      LOGICAL FIRST,SECOND
      SAVE FIRST,SECOND
      DATA FIRST,SECOND /.TRUE.,.TRUE./
*
*
*       Obtain energy error in case routine ADJUST not called recently.
      IF (TIME.GE.TADJ.OR.TIME.LE.0.0D0) GO TO 10
*
*       Predict X & XDOT for all particles (except unperturbed pairs).
      CALL XVPRED(IFIRST,NTOT)
*
*       Obtain the total energy at current time (resolve all KS pairs).
      CALL ENERGY
*
*       Include KS, triple & quad, mergers, mass loss, collisions & chain.
      ETOT = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + EMDOT + ECOLL
     &                                                         + ECDOT
      IF (NCH.GT.0) THEN
          ETOT = ETOT + ECH
      END IF
*
*       Update energies and form the relative error (divide by ZKIN or ETOT).
      BE(2) = BE(3)
      BE(3) = ETOT
      DE = BE(3) - BE(2)
      DETOT = DETOT + DE
      DE = DE/MAX(ZKIN,ABS(ETOT))
*       Save sum of relative energy error for main output and accumulate DE.
      ERROR = ERROR + DE
      ERRTOT = ERRTOT + DE
*
*       Find density centre & core radius (Casertano & Hut, Ap.J. 298, 80).
      IF (N.GT.10.AND.KZ(29).EQ.0) THEN
          CALL CORE
      END IF
*
*       Check optional sorting of Lagrangian radii & half-mass radius.
      IF (KZ(7).GT.0) THEN
          CALL LAGR(RDENS)
      END IF
*
*       Initialize diagnostic variables.
   10 NP = 0
      IUNP = 0
      AMIN = 100.0
      MULT = 0
      SUM = 0.0
*
*       Find smallest semi-major axis and count unperturbed KS pairs.
      DO 20 IPAIR = 1,NPAIRS
          NP = NP + LIST(1,2*IPAIR-1)
          SEMI = -0.5*BODY(N+IPAIR)/H(IPAIR)
          IF (SEMI.GT.0.0) AMIN = MIN(AMIN,SEMI)
          IF (LIST(1,2*IPAIR-1).EQ.0) IUNP = IUNP + 1
          IF (NAME(N+IPAIR).LT.-2*NZERO) MULT = MULT + 1
          SUM = SUM + BODY(N+IPAIR)**2
   20 CONTINUE
*
*       Count number of single stars (any NSUB contributes one).
      NS = 0
      DO 30 I = IFIRST,N
          IF (BODY(I).GT.0.0D0) NS = NS + 1
          SUM = SUM + BODY(I)**2
   30 CONTINUE
      NS = NS - NSUB
      NEFF = ZMASS**2/SUM
*
*       Set density centre displacement.
      RD = SQRT(RDENS(1)**2 + RDENS(2)**2 + RDENS(3)**2)
*
*       Check print frequency indicator & optional model counter.
      NPRINT = NPRINT + 1
      IF (NPRINT.GT.NFIX.OR.TIME.LE.0.0) THEN
          NPRINT = 1
          IF (KZ(3).GT.0) MODEL = MODEL + 1
      END IF
*
*       Form binary & merger energy ratios.
      EB = EBIN/(ZKIN - POT)
      IF (EB.LT.0.0) EB = EBIN/BE(3)
      EM = EMERGE/(ZKIN - POT)
      CALL JACOBI(NTESC)
*
*       Print main output diagnostics.
      I6 = TSCALE*TTOT
*
      WRITE (6,40)  TTOT, N, NPAIRS, NMERGE, MULT, NS, NSTEPI, NSTEPU,
     &              ERROR, BE(3)
   40 FORMAT (//,' T =',F7.1,'  N =',I6,'  KS =',I5,'  NM =',I3,
     &           '  MM =',I2,'  NS =',I6,'  NSTEPS =',I11,I12,
     &           '  DE =',F10.6,'  E =',F10.6)
*
      IF (KZ(21).GT.0) THEN
          CALL CPUTIM(TCOMP)
          DMIN1 = MIN(DMIN1, DMIN2, DMIN3, DMIN4, DMINC)
          WRITE (6,45)  NRUN, MODEL, TCOMP, DMIN1, DMIN2, AMIN, RMAX,
     &                  NBLOCK, NIRECT, NURECT, NEFF
   45     FORMAT (/,' RUN =',I3,'  M# =',I3,'  CPU =',F7.1,
     &              '  DMIN =',1P,2E8.1,'  AMIN =',E8.1,'  RMAX =',E8.1,
     &              '  NBLOCK =',I10,'  NIRECT =',I2,'  NURECT =',I2, 
     &              '  NEFF =',I6)
      END IF
*
      WRITE (6,50)
   50 FORMAT (/,'    <R>  RTIDE  RDENS   RC     NC   MC   RHOD   RHOM',
     &                                        '    UN  NPT  RCM    VCM',
     &                '        AZ     EB/E   EM/E   TCR      T6  NTESC')
*
      WRITE (6,55)  RSCALE, RTIDE, RD, RC, NC, ZMC, RHOD, RHOM, IUNP,
     &              NP, CMR(4), CMRDOT(4), AZ, EB, EM, TCR, I6, NTESC
   55 FORMAT (' #1',F5.2,F6.1,F7.2,F7.3,I5,F7.3,F6.1,F7.1,2I5,F7.3,
     &                                        F8.4,F10.6,2F7.3,F6.2,2I7)
*
      WRITE (6,65)
   65 FORMAT (/,6X,'NKSTRY  NKSREG  NKSHYP     NKSPER  NPRECT  NKSMOD',
     &             '   NTTRY  NTRIP  NQUAD  NCHAIN  NMERG  NEWHI',
     &             '  NSTEPC    NBCALL    NTPERT    NWARN  NHI')
      WRITE (6,70)  NKSTRY, NKSREG,  NKSHYP, NKSPER, NPRECT, NKSMOD,
     &              NTTRY, NTRIP, NQUAD, NCHAIN, NMERG, NEWHI, NSTEPC,
     &              NBCALL, NTPERT, NWARN, NHI
   70 FORMAT (' #2',I9,2I8,I11,3I8,2I7,I8,2I7,I8,2I10,I9,I5)
*
*       Check output for mass loss or tidal capture.
      IF (KZ(19).GT.0.OR.KZ(27).GT.0) THEN
          CALL EVENTS
      END IF
*
*       Obtain half-mass radii for two groups (NAME <= NZERO/5 & > NZERO/5).
      IF (KZ(7).EQ.6.AND.TTOT.GE.TCRIT.AND.BODY1.GT.2.0*BODYM.AND.
     &    KZ(5).NE.3.AND.KZ(5).NE.4) THEN
          CALL LAGR2(RDENS)
      END IF
*
*       Include diagnostics about cluster orbit in general external field.
      IF (KZ(14).EQ.3) THEN
          GZ = RG(1)*VG(2) - RG(2)*VG(1)
          SX = RBAR/1000.0
          WRITE (6,80)  NTAIL, NSTAIL, (RG(K)*SX,K=1,3),
     &                                 (VG(K)*VSTAR,K=1,3),
     &                  GZ, ETIDE
   80     FORMAT (/,5X,'CLUSTER ORBIT    NT NST RG VG JZ ET ',
     &                              I5,I9,3F7.2,2X,3F7.1,1P,E16.8,E10.2)
      END IF
*
*       Reset minimum encounter distances & maximum apocentre separation.
      DMIN2 = 100.0
      DMIN3 = 100.0
      DMIN4 = 100.0
      DMINC = 100.0
      RMAX = 0.0
*
*       Check integer overflows (2^{32} or 2.1 billion).
      IF (NSTEPI.GT.2000000000) THEN
          NSTEPI = 0
          NIRECT = NIRECT + 1
      END IF
      IF (NSTEPU.GT.2000000000) THEN
          NSTEPU = 0
          NURECT = NURECT + 1
      END IF
      IF (NBLOCK.GT.2000000000) THEN
          NBLOCK = 0
      END IF
      IF (NSTEPS.GT.2000000000) THEN
          NSTEPS = 0
      END IF
      IF (NBCALL.GT.2000000000) THEN
          NBCALL = 0
      END IF
      IF (NSTAIL.GT.2000000000) THEN
          NSTAIL = 0
      END IF
*
*       Ensure NLIST does not become large for block-step version.
      IF (TIME.LE.TBLOCK) TLIST = 0.0
*
*       Exit if error exceeds restart tolerance (TIME < TADJ means no CHECK).
      IF (ABS(ERROR).GT.5.0*QE.AND.TIME.LT.TADJ) GO TO 100
*
*       Check optional analysis & output of KS binaries.
      IF (KZ(8).GT.0.OR.NPAIRS.GT.0) THEN
          CALL BINOUT
      END IF
*
*       Include optional diagnostics of block-steps.
      IF (KZ(33).GT.0) THEN
          CALL LEVELS
      END IF
*
*       Check optional output of additional cluster parameters.
      IF (KZ(21).GT.1) THEN
          CALL CENTRE
      END IF
*
*       Include optional output for global cluster analysis.
      IF (KZ(14).GT.2.OR.KZ(21).GT.3) THEN
          CALL GLOBAL
      END IF
*
*       Check optional output of single bodies & binaries.
      CALL BODIES
*
*       See whether to write data bank of binary diagnostics on unit 9 or 4.
      IF (KZ(8).GE.2.AND.NPAIRS.GT.0.AND.TIME.GE.TPLOT) THEN
          CALL BINDAT
          IF (KZ(8).GE.2) THEN
              CALL DEGEN(1,NPAIRS,0)
              CALL HIDAT
          END IF
      END IF
      IF (KZ(8).GT.4) CALL PERIOD
*
*       Check optional diagnostics of evolving stars.
      IF (KZ(25).GT.0.AND.KZ(19).GE.3.AND.TIME.GE.TPLOT.AND.
     &    BODY1.GT.2.0*BODYM) THEN
          CALL HRPLOT
*
*       Reset GRAPE and send all particles (temporary measure).
*
*         CALL GPWIPE(GPID,TIME)
*         CALL GPSEND
*
      END IF
*
*       Check optional writing of data on unit 3 (frequency NFIX). 
      IF (KZ(3).EQ.0.OR.NPRINT.NE.1) GO TO 100
      IF (KZ(3).GT.1) GO TO 99
*
      AS(1) = TTOT
      AS(2) = FLOAT(NPAIRS)
      AS(3) = RBAR
      AS(4) = ZMBAR
      AS(5) = RTIDE
      AS(6) = TIDAL(4)
      AS(7) = RDENS(1)
      AS(8) = RDENS(2)
      AS(9) = RDENS(3)
      AS(10) = TTOT/TCR
      AS(11) = TSCALE
      AS(12) = VSTAR
      AS(13) = RC
      AS(14) = NC
      AS(15) = TURN
      AS(16) = RHOM
      AS(17) = RHOD
      AS(18) = RSCALE
      AS(19) = 0.0
      AS(20) = DMIN1
*
*       Convert masses, coordinates, velocities, RHO & PHI to REAL*4.
      DO 90 I = 1,NTOT
          BODYS(I) = BODY(I)
          RHO1(I) = RHO(I)
          PHI1(I) = PHI(I)
          DO 85 K = 1,3
              XS(K,I) = X(K,I)
              VS(K,I) = XDOT(K,I)
   85     CONTINUE
   90 CONTINUE
*
*       Copy density at c.m. for any KS components.
      DO 92 I = 1,NPAIRS
          RHO(2*I-1) = RHO(N+I)
          RHO(2*I) = RHO(N+I)
   92 CONTINUE
*
*       Replace any ghosts by actual M, R & V (including 2 binaries).
      DO 95 JPAIR = 1,NPAIRS
          J2 = 2*JPAIR
          J1 = J2 - 1
          ICM = N + JPAIR
*       Determine merger & ghost index for negative c.m. name.

          IF (NAME(ICM).LT.0.AND.BODY(ICM).GT.0.0) THEN
              CALL FINDJ(J1,J,IM)
*       Note: J is ghost index and IM is merger index.
              IF (J.LE.0) GO TO 95
              BODYS(J1) = CM(1,IM)
              BODYS(J) = CM(2,IM)
              ZMB = CM(1,IM) + CM(2,IM)
*       Form global coordinates and velocities from c.m. with XREL & VREL.
              DO K = 1,3
                  X1(K,1) = X(K,J1) + CM(2,IM)*XREL(K,IM)/ZMB
                  X1(K,2) = X(K,J1) - CM(1,IM)*XREL(K,IM)/ZMB
                  V1(K,1) = XDOT(K,J1) + CM(2,IM)*VREL(K,IM)/ZMB
                  V1(K,2) = XDOT(K,J1) - CM(1,IM)*VREL(K,IM)/ZMB
*
                  XS(K,J1) = X1(K,1)
                  XS(K,J)  = X1(K,2)
                  VS(K,J1) = V1(K,1)
                  VS(K,J)  = V1(K,2)
              END DO
*       Look for ghosts of possible second (i.e. outer) merged binary.
              IF (NAME(J).GT.NZERO) THEN
                  ICM2 = 0
                  DO  JJ = N+1,NTOT
                      IF (NAME(JJ).EQ.NAME(J)) ICM2 = JJ
                  END DO
*       Treat the second binary using inactive KS variables.
                  IF (ICM2.GT.0) THEN
                      IPAIR = ICM2 - N
                      I1 = 2*IPAIR - 1
                      I2 = I1 + 1
                      BODYS(I1) = CM(3,IM)
                      BODYS(I2) = CM(4,IM)
*       Copy KS variables to local scalars.
                      DO K = 1,4
                          UI(K) = U(K,IPAIR)
                          VI(K) = UDOT(K,IPAIR)
                      END DO
*       Transform to physical variables and multiply by 4 (momentum formula).
                      CALL KSPHYS(UI,VI,XREL2,VREL2)
                      ZM = CM(3,IM) + CM(4,IM)
                      DO K = 1,3
                          VREL2(K) = 4.0*VREL2(K)
                          X1(K,3) = X(K,J2) + CM(4,IM)*XREL2(K)/ZM
                          X1(K,4) = X(K,J2) - CM(3,IM)*XREL2(K)/ZM
                          V1(K,3) = XDOT(K,J2) + CM(4,IM)*VREL2(K)/ZM
                          V1(K,4) = XDOT(K,J2) - CM(3,IM)*VREL2(K)/ZM
*
                          XS(K,I1) = X1(K,3)
                          XS(K,I2)  = X1(K,4)
                          VS(K,I1) = V1(K,3)
                          VS(K,I2)  = V1(K,4)
                          XS(K,ICM2) = X(K,J2)
                          VS(K,ICM2) = XDOT(K,J2)
                      END DO
                  END IF
              END IF
          END IF
   95 CONTINUE
*
*       Check modification for chain regularization (case NAME(ICM) = 0).
      IF (NCH.GT.0) THEN
          CALL CHDATA(XJ,VJ,BODYJ)
          DO 98 L = 1,NCH
*       Copy global address from common JLIST (set in CHDATA).
              J = JLIST(L)
              BODYS(J) = BODYJ(L)
              DO 97 K = 1,3
                  XS(K,J) = XJ(K,L)
                  VS(K,J) = VJ(K,L)
   97         CONTINUE
   98     CONTINUE
      END IF
*
*       Split into WRITE (6) NTOT & WRITE (3) ..  if disc instead of tape.
      IF (FIRST) THEN
          OPEN (UNIT=3,STATUS='NEW',FORM='UNFORMATTED',FILE='OUT3')
          FIRST = .FALSE.
      END IF
      NK = 20
      WRITE (3)  NTOT, MODEL, NRUN, NK
      WRITE (3)  (AS(K),K=1,NK), (BODYS(J),J=1,NTOT),
     &           ((XS(K,J),K=1,3),J=1,NTOT), ((VS(K,J),K=1,3),J=1,NTOT),
     &           (RHO1(J),J=1,NTOT),(PHI1(J),J=1,NTOT),
     &           (NAME(J),J=1,NTOT),(KSTAR(J),J=1,NTOT)
      CALL FLUSH(3)
*     CLOSE (UNIT=3)
*
*       Include all stars in same file (KZ(3) > 1; astrophysical units). 
   99 IF (KZ(3).GT.1.AND.NTAIL.GT.0) THEN
          IF (SECOND) THEN
              OPEN (UNIT=34,STATUS='NEW',FORM='FORMATTED',FILE='OUT34')
              SECOND = .FALSE.
          END IF
          NP = 0
*       Copy cluster members with respect to density centre.
          DO 120 I = IFIRST,NTOT
              IF (BODY(I).LE.0.0D0) GO TO 120
              NP = NP + 1
              DO 118 K = 1,3
                  XS(K,NP) = (X(K,I) - RDENS(K))*RBAR
                  VS(K,NP) = XDOT(K,I)*VSTAR
  118         CONTINUE
              BODYS(NP) = BODY(I)*SMU
  120     CONTINUE
          N1 = NP
*       Add tidal tail in the same frame.
          DO 130 I = ITAIL0,NTTOT
              NP = NP + 1
              DO 125 K = 1,3
                  XS(K,NP) = (X(K,I) - RG(K) - RDENS(K))*RBAR
                  VS(K,NP) = (XDOT(K,I) - VG(K))*VSTAR
  125         CONTINUE
              BODYS(NP) = BODY(I)*SMU
  130     CONTINUE
          WRITE (34,140)  NP, N1, (TIME+TOFF)*TSCALE
  140     FORMAT (' ',2I6,F8.1)
          DO 150 I = 1,NP
              WRITE (34,145) (XS(K,I),K=1,3), (VS(K,I),K=1,3), BODYS(I)
  145         FORMAT (' ',3F10.3,3F8.1,F7.2)
  150     CONTINUE
          CALL FLUSH(34)
      END IF
*
*       Update output interval and initialize the corresponding error.
  100 TPRINT = TPRINT + DELTAT
      ERROR = 0.0D0
*
      RETURN
*
      END
