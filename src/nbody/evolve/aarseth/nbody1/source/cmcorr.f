      SUBROUTINE CMCORR
*
*
*       Center of mass corrections.
*       ---------------------------
*
      INCLUDE 'common1.h'
      REAL*8  DX0,DV0
*
*
*       Initialize centre of mass variables.
      DO 10 K = 1,3
          CMR(K) = 0.0
          CMRDOT(K) = 0.0
   10 CONTINUE
*
*       Form c.m. coordinate & velocity displacements.
      DO 20 I = 1,N
          DO 15 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   15     CONTINUE
   20 CONTINUE
*
      DO 30 K = 1,3
          CMR(K) = CMR(K)/ZMASS
          CMRDOT(K) = CMRDOT(K)/ZMASS
   30 CONTINUE
*
*       Apply c.m. corrections and accumulate energy changes.
      ERRV = 0.0
      DO 40 I = 1,N
          DO 35 K = 1,3
              VI2 = XDOT(K,I)**2
              X(K,I) = X(K,I) - CMR(K)
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)
              ERRV = ERRV + BODY(I)*(XDOT(K,I)**2 - VI2)
   35     CONTINUE
   40 CONTINUE
*
*       Adjust the total energy to new kinetic energy.
      BE(3) = BE(3) + 0.5*ERRV
*
*       Redetermine X0 & X0DOT consistently with current corrected X & XDOT.
      DO 60 I = 1,N
          DT = TIME - T0(I)
          A1 = 0.05*DT
          A2 = 0.25*DT
          A3 = (T0(I) - T1(I)) + (T0(I) - T2(I))
          DO 55 K = 1,3
              F2DOTK = D3(K,I)*A3 + D2(K,I)
              F3DOTK = D3(K,I)
              DV0 = (((F3DOTK*A2 + ONE3*F2DOTK)*DT +
     &                                3.0*FDOT(K,I))*DT + 2.0*F(K,I))*DT
              X0DOT(K,I) = XDOT(K,I) - DV0
              DX0 = ((((F3DOTK*A1 + ONE12*F2DOTK)*DT + FDOT(K,I))*DT +
     &                                       F(K,I))*DT + X0DOT(K,I))*DT
              X0(K,I) = X(K,I) - DX0
   55     CONTINUE
   60 CONTINUE
*
      RETURN
*
      END
