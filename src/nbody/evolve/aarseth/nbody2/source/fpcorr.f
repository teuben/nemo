      SUBROUTINE FPCORR(I,NBLOSS,NBGAIN,NNB,XI,YI,ZI,XIDOT,DTR,DT1R)
*
*
*       Force polynomial corrections.
*       -----------------------------
*
      INCLUDE 'common2.h'
      REAL*8  XI,YI,ZI
      REAL*4  XIDOT(3),SAVE1(3),SAVE2(3),SAVE3(3),A(12),
     &        F1DOT(3),F2DOT(3),F3DOT(4)
*
*
*       See whether there has been a change of neighbours.
      NBFLUX = NBLOSS + NBGAIN
      IF (NBFLUX.EQ.0) GO TO 60
*
*       Initialize the derivative corrections.
      DO 10 K = 1,3
          SAVE1(K) = 0.0
          SAVE2(K) = 0.0
          SAVE3(K) = 0.0
   10 CONTINUE
*
*       Form compact list of NBLOSS & NBGAIN.
      IF (NBGAIN.GT.0) THEN
          NNB0 = LIST(1,I)
          DO 15 L = 1,NBGAIN
              JLIST(NBLOSS+L) = JLIST(NNB0+L)
   15     CONTINUE
      END IF
*
*       Accumulate derivative corrections.
      L = 1
   20 J = JLIST(L)
      S = TIME - T0(J)
      S3 = 3.0*S
*
*       Predict velocity and force of body #J to order FDOT.
      DO 25 K = 1,3
          A(K+3) = (FDOT(K,J)*S3 + 2.0*F(K,J))*S + X0DOT(K,J) - XIDOT(K)
          A(K+6) = 2.0*(FDOT(K,J)*S3 + F(K,J) - F(K,I))
          A(K+9) = 6.0*(FDOT(K,J) - FDOT(K,I))
   25 CONTINUE
*
      A(1) = X(1,J) - XI
      A(2) = X(2,J) - YI
      A(3) = X(3,J) - ZI
*
      A13 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
      A14 = BODY(J)*A13*SQRT(A13)
      A15 = (A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A13
      A16 = A15*A15
      A17 = 3.0*A15
      A18 = 6.0*A15
      A19 = 9.0*A15
      A20 = (A(4)*A(4) + A(5)*A(5) + A(6)*A(6) + A(1)*A(7) + A(2)*A(8)
     &                                            + A(3)*A(9))*A13 + A16
      A21 = 9.0*A20
      A20 = 3.0*A20
      A22 = (9.0*(A(4)*A(7) + A(5)*A(8) + A(6)*A(9)) + 3.0*(A(1)*A(10)
     &             + A(2)*A(11) + A(3)*A(12)))*A13 + A17*(A20 - 4.0*A16)
*
      DO 30 K = 1,3
          F1DOTK = A(K+3) - A17*A(K)
          F2DOT(K) = (A(K+6) - A18*F1DOTK - A20*A(K))*A14
          F3DOT(K) = (A(K+9) - A21*F1DOTK - A22*A(K))*A14 - A19*F2DOT(K)
          F1DOT(K) = F1DOTK*A14
   30 CONTINUE
*
*       Change the sign for NBLOSS contributions.
      IF (L.LE.NBLOSS) THEN
          DO 35 K = 1,3
              F1DOT(K) = -F1DOT(K)
              F2DOT(K) = -F2DOT(K)
              F3DOT(K) = -F3DOT(K)
   35     CONTINUE
      END IF
*
*       Include derivative corrections from losses & gains.
      DO 40 K = 1,3
          SAVE1(K) = SAVE1(K) + F1DOT(K)
          SAVE2(K) = SAVE2(K) + F2DOT(K)
          SAVE3(K) = SAVE3(K) + F3DOT(K)
   40 CONTINUE
*
      L = L + 1
      IF (L.LE.NBFLUX) GO TO 20
*
      NBCORR = NBCORR + 1
*       Copy new list of neighbours from working space to COMMON.
      LIST(1,I) = NNB
      DO 45 L = 2,NNB+1
          LIST(L,I) = ILIST(L)
   45 CONTINUE
*
*       Perform corrections to irregular and regular force differences.
      DT = TIME - T1(I)
      DT1 = TIME - T2(I)
      S = DT + DT1
      SR = DTR + DT1R
      DO 50 K = 1,3
          F1DOTK = (D3(K,I)*DT1 + D2(K,I))*DT + D1(K,I) + SAVE1(K)
          F2DOTK = D3(K,I)*S + D2(K,I) + 0.5D0*SAVE2(K)
          F3DOTK = D3(K,I) + ONE6*SAVE3(K)
          D1(K,I) = (F3DOTK*DT - F2DOTK)*DT + F1DOTK
          D2(K,I) = F2DOTK - F3DOTK*S
          D3(K,I) = F3DOTK
          F1DOTK = (D3R(K,I)*DT1R + D2R(K,I))*DTR + D1R(K,I) - SAVE1(K)
          F2DOTK = D3R(K,I)*SR + D2R(K,I) - 0.5D0*SAVE2(K)
          F3DOTK = D3R(K,I) - ONE6*SAVE3(K)
          D1R(K,I) = (F3DOTK*DTR - F2DOTK)*DTR + F1DOTK
          D2R(K,I) = F2DOTK - F3DOTK*SR
          D3R(K,I) = F3DOTK
   50 CONTINUE
*
   60 RETURN
*
      END
