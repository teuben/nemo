      SUBROUTINE REGINT(I,XI,YI,ZI,XIDOT)
*
*
*       Regular integration.
*       --------------------
*
      INCLUDE 'common2.h'
      REAL*8  XI,YI,ZI,FIRR(3),FREG(3)
      REAL*4  XIDOT(3),F1DOT(4),F2DOT(4),F3DOT(4)
*
*
*       Set neighbour number, time-step factors & central distance.
      NNB0 = LIST(1,I)
      DT = TIME - T1(I)
      DT1 = TIME - T2(I)
      DTR = TIME - T0R(I)
      DTRIN = 1.0/DTR
      IRSKIP = 0
      RI2 = (XI - RDENS(1))**2 + (YI - RDENS(2))**2 + (ZI - RDENS(3))**2
*
*       Obtain irregular & regular force and determine current neighbours.
      RS2 = RS(I)**2
*       Take volume between inner and outer radius equal to basic sphere.
    1 RCRIT2 = 1.59*RS2
*       Set radial velocity factor for the outer shell.
      VRFAC = -0.1*RS2*DTRIN
*       Start count at 2 and subtract 1 at the end to avoid ILIST(NNB+1).
      NNB = 1
*
*       Initialize the force components.
      DO 5 K = 1,3
          FIRR(K) = 0.0
          FREG(K) = 0.0
    5 CONTINUE
*
*       Sum over all the particles.
      DO 10 J = 1,N
          A1 = X(1,J) - XI
          A2 = X(2,J) - YI
          A3 = X(3,J) - ZI
*       Predicted coordinates avoids spurious force differences.
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
*       First see whether the distance exceeds the outer shell radius.
          IF (RIJ2.GT.RCRIT2) GO TO 8
*
          IF (RIJ2.GT.RS2) THEN
              A4 = XDOT(1,J) - XIDOT(1)
              A5 = XDOT(2,J) - XIDOT(2)
              A6 = XDOT(3,J) - XIDOT(3)
              A7 = A1*A4 + A2*A5 + A3*A6
*       Accept member if maximum penetration factor exceeds 8 per cent.
              IF (A7.GT.VRFAC) GO TO 8
          ELSE
              IF (J.EQ.I) GO TO 10
          END IF
*
*       Increase neighbour counter and obtain current irregular force.
          NNB = NNB + 1
          ILIST(NNB) = J
          RIJ2 = RIJ2 + EPS2
          A5 = BODY(J)/(RIJ2*SQRT(RIJ2))
          FIRR(1) = FIRR(1) + A1*A5
          FIRR(2) = FIRR(2) + A2*A5
          FIRR(3) = FIRR(3) + A3*A5
          GO TO 10
*
*       Obtain the regular force (FIRR & FREG summed in double precision).
    8     RIJ2 = RIJ2 + EPS2
          A5 = BODY(J)/(RIJ2*SQRT(RIJ2))
          FREG(1) = FREG(1) + A1*A5
          FREG(2) = FREG(2) + A2*A5
          FREG(3) = FREG(3) + A3*A5
   10 CONTINUE
*
*       See whether an external regular force should be added.
      IF (KZ(15).GT.0) THEN
*       Choose between Plummer model & logarithmic potential.
          IF (KZ(15).EQ.1) THEN
              CALL XTRNL1(I,1,FREG)
          ELSE
              CALL XTRNL2(I,1,FREG)
          END IF
      END IF
*
      NNB = NNB - 1
      IF (NNB.EQ.0) THEN
*       Double the neighbour sphere and try again unless RS > 50*RSCALE.
          IF (RS(I).GT.50.0*RSCALE) THEN
              IRSKIP = 1
*       Assume small mass at centre for rare case of no neighbours.
              A6 = 0.01*BODYM/(RI2*SQRT(RI2))
              DO 25 K = 1,3
                  FIRR(K) = -A6*X(K,I)
   25         CONTINUE
              LIST(1,I) = 0
              GO TO 50
          ELSE
              RS2 = 1.59*RS2
          END IF
*
          RS(I) = SQRT(RS2)
          NBVOID = NBVOID + 1
          IF (RS(I).GT.10.0*RSCALE) IRSKIP = 1
          GO TO 1
      END IF
*
*       Check maximum neighbour number.
      IF (NNB.LE.NNBMAX) GO TO 40
*
*       Reduce search radius by cube root of conservative volume factor.
   30 NNB2 = ZNBMAX
      A1 = ZNBMAX/FLOAT(NNB)
      IF (RS(I).GT.5.0*RSCALE) THEN
          A1 = MIN(5.0*A1,0.9)
          IRSKIP = 1
      END IF
      RS2 = RS2*A1**0.66667
      RCRIT2 = 1.59*RS2
      RS(I) = SQRT(RS2)
      NNB1 = 0
*
      DO 35 L = 1,NNB
          J = ILIST(L+1)
          IF (L + NNB2.GT.NNB + NNB1) GO TO 32
*       Sum of neighbours (NNB1) and those left (NNB+1-L) set to NNB2.
          A1 = X(1,J) - XI
          A2 = X(2,J) - YI
          A3 = X(3,J) - ZI
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          IF (RIJ2.GT.RCRIT2) GO TO 34
          IF (RIJ2.GT.RS2) THEN
              A4 = XDOT(1,J) - XIDOT(1)
              A5 = XDOT(2,J) - XIDOT(2)
              A6 = XDOT(3,J) - XIDOT(3)
              A7 = A1*A4 + A2*A5 + A3*A6
              IF (A7.GT.VRFAC) GO TO 34
          END IF
   32     NNB1 = NNB1 + 1
          JLIST(NNB1+1) = J
          GO TO 35
*
*       Subtract neighbour force included above and add to regular force.
   34     RIJ2 = RIJ2 + EPS2
          A5 = BODY(J)/(RIJ2*SQRT(RIJ2))
          FIRR(1) = FIRR(1) - A1*A5
          FIRR(2) = FIRR(2) - A2*A5
          FIRR(3) = FIRR(3) - A3*A5
          FREG(1) = FREG(1) + A1*A5
          FREG(2) = FREG(2) + A2*A5
          FREG(3) = FREG(3) + A3*A5
   35 CONTINUE
*
*       Copy reduced membership to current list.
      DO 38 L = 2,NNB1+1
          ILIST(L) = JLIST(L)
   38 CONTINUE
      NNB = NNB1
      NBFULL = NBFULL + 1
*       See whether to reduce NNB further.
      IF (NNB.GT.NNBMAX) GO TO 30
*
*       Stabilize NNB between ZNBMIN & ZNBMAX by square root of contrast.
   40 A3 = ALPHA*SQRT(FLOAT(NNB)*RS(I))/RS2
      A3 = MIN(A3,ZNBMAX) 
      A4 = MAX(A3,ZNBMIN)/FLOAT(NNB)
*       Include inertial factor to prevent resonance oscillations in RS.
      IF ((A3 - FLOAT(NNB0))*(A3 - FLOAT(NNB)).LT.0.0) A4 = SQRT(A4)
*
*       Check option for unique density centre.
      IF (KZ(10).EQ.0) THEN
*       Modify volume ratio by radial velocity factor outside the core.
          IF (RI2.GT.RC2) THEN
              RIDOT = (XI - RDENS(1))*XIDOT(1) +
     &                (YI - RDENS(2))*XIDOT(2) +
     &                (ZI - RDENS(3))*XIDOT(3)
              A4 = A4*(1.0 + RIDOT*DTR/RI2)
          END IF
*       Impose a time-step dependent maximum change below 0.01*TCR.
          A5 = MIN(25.0*DTR/TCR,0.25)
      ELSE
          A5 = 0.25
      END IF
*
*       Restrict volume ratio by inertial factor of 25 per cent either way.
      A4 = A4 - 1.0
      IF (A4.GT.A5) THEN
          A4 = A5
      ELSE
          A4 = MAX(A4,-A5)
      END IF
*
*       Modify neighbour sphere radius by volume factor.
      IF (IRSKIP.EQ.0) THEN
          IF (RS(I).GT.50.0*RSCALE) GO TO 50
          A3 = ONE3*A4
          A1 = 1.0 + A3 - A3*A3
*       Second-order cube root expansion (maximum volume error < 0.3 %).
          IF (RS(I).GT.5.0*RSCALE) A1 = SQRT(A1)
*       Inertial factor for distant particle avoids oscillations in RS.
          RS(I) = A1*RS(I)
      END IF
*
*       Calculate radial velocity with respect to at most 3 neighbours.
      IF (NNB.LE.3) THEN
          A1 = 2.0*RS(I)
*
          DO 45 L = 1,NNB
              J = ILIST(L+1)
              RIJ = SQRT((XI - X(1,J))**2 + (YI - X(2,J))**2 +
     &                                      (ZI - X(3,J))**2)
              RSDOT = ((XI - X(1,J))*(XIDOT(1) - XDOT(1,J)) +
     &                 (YI - X(2,J))*(XIDOT(2) - XDOT(2,J)) +
     &                 (ZI - X(3,J))*(XIDOT(3) - XDOT(3,J)))/RIJ
*       Find smallest neighbour distance assuming constant regular step.
              A1 = MIN(A1,RIJ + RSDOT*STEPR(I))
   45     CONTINUE
*
*       Increase neighbour sphere if all members are leaving inner region.
          RS(I) = MAX(A1,1.1*RS(I))
      END IF
*
*       Check minimum neighbour sphere since last output (diagnostic only).
      RSMIN = MIN(RSMIN,RS(I))
*
*       Find loss or gain of neighbours at the same time.
   50 NBLOSS = 0
      NBGAIN = 0
*
*       Check case of zero old or new membership (skip if both are zero).
      IF (NNB0.EQ.0) THEN
          IF (NNB.EQ.0) GO TO 70
          LIST(2,I) = 0
      END IF
*
      L = 2
      LG = 2
*       Set termination value in ILIST(NNB+2) and save last list member.
      ILIST(NNB+2) = N + 1
      ILIST(1) = LIST(NNB0+1,I)
*
*       Compare old and new list members in locations L & LG.
   56 IF (LIST(L,I).EQ.ILIST(LG)) GO TO 58
*
*       Now check whether inequality means gain or loss.
      IF (LIST(L,I).GE.ILIST(LG)) THEN
          NBGAIN = NBGAIN + 1
          JLIST(NNB0+NBGAIN) = ILIST(LG)
*       Number of neighbour losses can at most be NNB0.
          L = L - 1
*       The same location will be searched again after increasing L below.
      ELSE
          NBLOSS = NBLOSS + 1
          J = LIST(L,I)
          JLIST(NBLOSS) = J
          LG = LG - 1
      END IF
*
*       See whether the last location has been checked.
   58 IF (L.LE.NNB0) THEN
          L = L + 1
          LG = LG + 1
*       Last value of second search index is NNB + 2 which holds N + 1.
          GO TO 56
      ELSE IF (LG.LE.NNB) THEN
          LG = LG + 1
          LIST(L,I) = N + 1
*       Last location of list holds termination value (saved in ILIST(1)).
          GO TO 56
      END IF
*
*       Set time intervals for corrector and update all regular force times.
   70 DT1R = TIME - T1R(I)
      DT2R = TIME - T2R(I) 
      T1PR = T0R(I) - T1R(I)
      T2PR = T0R(I) - T2R(I)
      T3PR = T0R(I) - T3R(I)
      DT3R = TIME - T3R(I)
      DT06 = 0.6*DTR
      DT12 = ONE12*DTR
      S2 = T1PR*T2PR
      S3 = S2*T3PR
      S4 = S2 + T3PR*(T1PR + T2PR)
      S5 = T1PR + T2PR + T3PR
      S6 = (((0.6666667*DT + S5)*DT06 + S4)*DT12 + ONE6*S3)*DT
      S7 = ((0.2*DTR + 0.25*S5)*DTR + ONE3*S4)*DTR + 0.5*S3
      T3R(I) = T2R(I)
      T2R(I) = T1R(I)
      T1R(I) = T0R(I)
      T0R(I) = TIME
*       Current regular step is reduced to nearest irregular step.
      A2 = 1.0/DT1R
      A3 = 1.0/DT2R
      A4 = DTR*DTR/DT3R
*
*       Evaluate the regular divided differences and save F, FI, FR & FDOT.
      DO 75 K = 1,3
*       Subtract change of neighbour force to get actual first difference.
          D1RK = (FREG(K) - (FI(K,I) - FIRR(K)) - FR(K,I))*DTRIN
          D2RK = (D1RK - D1R(K,I))*A2
          D3RK = (D2RK - D2R(K,I))*A3
          F4DOTK = (D3RK - D3R(K,I))*A4
          F(K,I) = 0.5*(FREG(K) + FIRR(K))
*       One-half the total force used for fast coordinate prediction.
          FI(K,I) = FIRR(K)
          FR(K,I) = FREG(K)
*       Current force components include possible change of neighbours.
          D1R(K,I) = D1RK
          D2R(K,I) = D2RK
          D3R(K,I) = D3RK
*       Convert from differences to sixth force derivative at current time.
          FDOT(K,I) = ONE6*((D3RK*DT1R + D2RK)*DTR + D1RK +
     &                             (D3(K,I)*DT1 + D2(K,I))*DT + D1(K,I))
*       Include semi-iteration for the regular force.
          X(K,I) = F4DOTK*S6 + X(K,I)
          X0(K,I) = X(K,I)
          X0DOT(K,I) = F4DOTK*S7 + X0DOT(K,I)
          XDOT(K,I) = X0DOT(K,I)
   75 CONTINUE
*
*       Correct force polynomials due to neighbour changes.
      CALL FPCORR(I,NBLOSS,NBGAIN,NNB,XI,YI,ZI,XIDOT,DTR,DT1R)
*
*       Obtain new regular integration step using composite expression.
*       STEPR = (ETAR*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
      SR = DTR + DT1R
  100 DO 102 K = 1,3
          F1DOT(K) = D2R(K,I)*DTR + D1R(K,I)
          F2DOT(K) = D3R(K,I)*SR + D2R(K,I)
          F3DOT(K) = D3R(K,I)
  102 CONTINUE
*
      FR2 = FREG(1)**2 + FREG(2)**2 + FREG(3)**2
      F1DOT(4) = F1DOT(1)**2 + F1DOT(2)**2 + F1DOT(3)**2
      F2DOT(4) = 4.0*(F2DOT(1)**2 + F2DOT(2)**2 + F2DOT(3)**2)
      F3DOT(4) = 36.0*(F3DOT(1)**2 + F3DOT(2)**2 + F3DOT(3)**2)
*
*       Form new step by relative criterion (extra SQRT for large F3DOT).
      IF (F3DOT(4).LT.1.0E+20) THEN
          A1 = (SQRT(FR2*F2DOT(4)) + F1DOT(4))/
     &                              (SQRT(F1DOT(4)*F3DOT(4)) + F2DOT(4))
      ELSE
          A1 = (SQRT(FR2*F2DOT(4)) + F1DOT(4))/
     &                        (SQRT(F1DOT(4))*SQRT(F3DOT(4)) + F2DOT(4))
      END IF
*
*       Restrict increase of regular time-step by stability factor 1.2.
      STEPR(I) = MIN(SQRT(ETAR*A1),1.2*STEPR(I))
      NSTEPR = NSTEPR + 1
*
      RETURN
*
      END
