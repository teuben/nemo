      SUBROUTINE XVPRED(I1,I2)
*
*
*       Prediction of coordinates & velocities.
*       ---------------------------------------
*
      INCLUDE 'common2.h'
*
*
*       Predict current coordinates & velocities to order F3DOT.
      DO 10 I = I1,I2
          DT = TIME - T0(I)
          A1 = 0.05*DT
          A2 = 0.25*DT
          A3 = (T0(I) - T1(I)) + (T0(I) - T2(I))
          A4 = (T0(I) - T0R(I)) + (T0(I) - T1R(I)) + (T0(I) - T2R(I))
          DO 5 K = 1,3
              F2DOTK = D3R(K,I)*A4 + D2R(K,I) + D3(K,I)*A3 + D2(K,I)
              F3DOTK = D3R(K,I) + D3(K,I)
              X(K,I) = ((((F3DOTK*A1 + ONE12*F2DOTK)*DT +
     &                     FDOT(K,I))*DT + F(K,I))*DT + X0DOT(K,I))*DT +
     &                                                           X0(K,I)
              XDOT(K,I) = (((F3DOTK*A2 + ONE3*F2DOTK)*DT +
     &                              3.0*FDOT(K,I))*DT + 2.0*F(K,I))*DT +
     &                                                        X0DOT(K,I)
    5     CONTINUE
   10 CONTINUE
*
      RETURN
*
      END
