      SUBROUTINE BINOUT
*
*
*       Binary analysis & output.
*       -------------------------
*
      INCLUDE 'common4.h'
      COMMON/ECHAIN/  ECH
      INTEGER  JD(15),JEB(15),JE(10)
*
*
*       Define semi-major axis & binding energy of hard binary (= KT).
      A0 = 1.0/FLOAT(NZERO)
      EB0 = -0.25/FLOAT(NZERO)
*       Initialize counters & variables.
      DO 10 J = 1,15
          JD(J) = 0
          JEB(J) = 0
   10 CONTINUE
      DO 20 J = 1,10
          JE(J) = 0
   20 CONTINUE
      EBMAX = 0.0
      DISP = 0.0
      EMAX = 0.0
      JOR = 0
      JEX = 0
      JC = 0
      JLAST = 0
      KLAST = 0
      E(1) = 0.0
      E(2) = 0.0
      NPOP(1) = 0
      NPOP(2) = 0
      NPOP(3) = N - 2*NPAIRS
*
*       Obtain relevant distributions of KS binaries.
      DO 50 J = 1,NPAIRS
          I = N + J
          KST = KSTAR(I)
          IF (H(J).GE.0.0.OR.BODY(I).EQ.0.0) GO TO 50
          J1 = 2*J - 1
          J2 = 2*J
*
*       Count original & exchanged binaries and core members.
          IF (IABS(NAME(J1) - NAME(J2)).EQ.1) THEN
              JOR = JOR + 1
          ELSE
              IF (MIN(NAME(J1),NAME(J2)).LE.2*NBIN0) JEX = JEX + 1
          END IF
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          IF (RI2.LT.RC**2) JC = JC + 1
*
*       Adopt logarithmic distribution of semi-major axis (factor 2).
          SEMI = -0.5D0*BODY(I)/H(J)
          IF (SEMI.LT.A0) THEN
              K = 2 + LOG10(A0/SEMI)/LOG10(2.0)
              K = MIN(K,10)
              JLAST = MAX(JLAST,K)
              JD(K) = JD(K) + 1
          ELSE
              JD(1) = JD(1) + 1
          END IF
*
*       Form eccentricity dispersion & histogram.
          ECC2 = (1.0 - R(J)/SEMI)**2 + TDOT2(J)**2/(BODY(I)*SEMI)
          DISP = DISP + ECC2
          ECC = SQRT(ECC2)
          IF (ECC.GT.EMAX) EMAX = ECC
          K = 1 + 10.0*ECC
          IF (K.LE.10) JE(K) = JE(K) + 1
*
*       Set up logarithmic distribution of binding energy (factor 2).
          EB = -0.5*BODY(J1)*BODY(J2)/SEMI
          IF (EB.LT.EB0) THEN
              K = 2 + LOG10(EB/EB0)/LOG10(2.0)
              K = MIN(K,14)
              KLAST = MAX(KLAST,K)
              JEB(K) = JEB(K) + 1
              EBMAX = MIN(EBMAX,EB)
          ELSE
              JEB(1) = JEB(1) + 1
          END IF
*
*       Define flag to distinguish primordial or non-primordial binary.
          IP = 1
          IF (LIST(2,J2).EQ.0) IP = 2
*
*       Sum the respective energies & populations.
          E(IP) = E(IP) + EB
          NPOP(IP) = NPOP(IP) + 1
*
*       Perform consistency check in case of over-writing.
          IF (LIST(2,J2).NE.-1.AND.LIST(2,J2).NE.0) THEN
              WRITE (6,30)  J, LIST(1,J1), LIST(1,J2), LIST(2,J2),
     &                      NAME(N+J)
   30         FORMAT (/,' FLAG PROBLEM!   PAIR NP NP2 FLAG NAME ',5I6)
          END IF
*
*       Produce special diagnostics for new binaries (unit 8).
          IF (IP.EQ.2.AND.H(J).LT.0.0) THEN
              VR = ((X(1,I) - RDENS(1))*XDOT(1,I) +
     &              (X(2,I) - RDENS(2))*XDOT(2,I) +
     &              (X(3,I) - RDENS(3))*XDOT(3,I))/SQRT(RI2)
              K = 0
              IF (NAME(I).LT.0) K = -1
              GX = GAMMA(J)*(SEMI*(1.0 + ECC)/R(J))**3
              WRITE (8,35)  TTOT, NAME(J1), NAME(J2), LIST(2,J2), K,
     &                      BODY(J1), BODY(J2), EB, SEMI, ECC, GX,
     &                      SQRT(RI2), VR
   35         FORMAT (' T =',F7.1,'  NAME = ',2I6,2I3,'  M =',2F8.4,
     &                         '  EB =',F10.5,'  A =',F8.5,'  E =',F5.2,
     &                        '  GX =',F6.3,'  RI =',F6.2,'  VR =',F4.1)
              CALL FLUSH(8)
          END IF
*
*       Rectify KSTAR for inconsistent non-circular and circular binaries.
          IF (ECC.GT.0.01.AND.KSTAR(I).GE.10) THEN
              WRITE (6,40)  J, NAME(J1), NAME(J2), KSTAR(J1), KSTAR(J2),
     &                      KSTAR(I), LIST(1,J1), ECC, SEMI, GAMMA(J)
   40         FORMAT (' NON-CIRCULAR:    KS NAM K* NP E A G ',
     &                                   I4,2I6,4I4,F7.3,1P,2E9.1)
              KSTAR(I) = 0
          ELSE IF (ECC.LT.0.002.AND.KSTAR(I).LE.0.AND.
     &             KZ(19).GT.0) THEN
              WRITE (6,45)  J, NAME(J1), NAME(J2), KSTAR(J1), KSTAR(J2),
     &                      KSTAR(I), LIST(1,J1), ECC, SEMI, GAMMA(J)
   45         FORMAT (' CIRCULARIZED:    KS NAM K* NP E A G ',
     &                                   I4,2I6,4I4,F7.3,1P,2E9.1)
              KSTAR(I) = 10
              NCIRC = NCIRC + 1
              IF (KZ(34).GT.0) THEN
                  CALL TRFLOW(J,DTR)
                  TEV(I) = TIME + DTR
              END IF
          END IF
*
   50 CONTINUE
*
*       Set fractional gain of binding energy (initial energy = -0.25).
      IF (TIME.LE.0.0D0) THEN
          JOR = NBIN0
          DB = 0.0
      ELSE
          E(9) = EMERGE
          DB = -4.0*(E(1) + E(2) + E(5) + E(7) + E(9) + ECOLL + EMDOT
     &                                                + EKICK - EBIN0)
      END IF
*
*       Print eccentricity & energy distributions.
      IF (NPAIRS.EQ.0) GO TO 70 
      WRITE (6,60)  JOR, JEX, DB, SBCOLL, BBCOLL, CHCOLL, JC,
     &              (JD(J),J=1,JLAST)
   60 FORMAT (/,' OR =',I4,'  EX =',I3,'  DB =',F7.3,'  SB =',F8.4,
     &          '  BB =',F8.4,'  CH =',F8.4,'  NC =',I3,
     &          '  N(A) =',14I4)
*
      IF (DISP.GT.0.0D0) DISP = SQRT(DISP/FLOAT(NPAIRS))
      WRITE (6,65)  DISP, EMAX, (NPOP(J),J=1,9), (JEB(K),K=1,KLAST)
   65 FORMAT (' <E> =',F5.2,'  EMAX =',F6.3,'  NPOP =',I4,I3,2I6,I4,4I3,
     &                                                 '  EB/KT =',14I4)
*
*       Form the basic internal energy (binaries & single particles).
   70 ETOT = 0.0
      E(3) = ZKIN - POT + ETIDE
      DO 80 J = 1,3
          ETOT = ETOT + E(J)
   80 CONTINUE
*
*       Sum all separate energy components to monitor conservation.
      ETOT = ETOT + ESUB + EMERGE + EMDOT + ECDOT + ECOLL
*
*       Include energy of chain if NCH > 0.
      IF (NCH.GT.0) THEN
          ETOT = ETOT + ECH
      END IF
*
      WRITE (6,90)  (E(J),J=1,10), ETOT
   90 FORMAT (' ENERGIES   ',10F9.5,'  ETOT =',F10.6)
*
      RETURN
*
      END
