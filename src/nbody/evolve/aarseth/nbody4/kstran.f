      SUBROUTINE KSTRAN(IPAIR,I1,I,BODYIN,RI,UI,UIDOT,XI,VI)
*
*
*       KS transformation.
*       ------------------
*
      INCLUDE 'common4.h'
      REAL*8  UI(4),UIDOT(4),XI(6),VI(6),RDOT(3),A1(3,4)
*
*
*       Form relative coordinates from explicit KS transformation.
      Q1 = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
      Q2 = UI(1)*UI(2) - UI(3)*UI(4)
      Q3 = UI(1)*UI(3) + UI(2)*UI(4)
      Q2 = Q2 + Q2
      Q3 = Q3 + Q3
*
*       Assign global coordinates of regularized components.
      A2 = BODY(I1+1)*BODYIN
      XI(1) = X(1,I) + A2*Q1
      XI(2) = X(2,I) + A2*Q2
      XI(3) = X(3,I) + A2*Q3
      XI(4) = XI(1) - Q1
      XI(5) = XI(2) - Q2
      XI(6) = XI(3) - Q3
*
*       Set current transformation matrix.
      CALL MATRIX(UI,A1)
*
*       Obtain relative velocities from KS transformation.
      RINV = 2.0D0/RI
      DO 30 L = 1,3
          RDOT(L) = 0.0D0
          DO 25 K = 1,4
              RDOT(L) = RDOT(L) + A1(L,K)*UIDOT(K)
   25     CONTINUE
          RDOT(L) = RDOT(L)*RINV
   30 CONTINUE
*
*       Set global velocities of KS components.
      VI(1) = XDOT(1,I) + A2*RDOT(1)
      VI(2) = XDOT(2,I) + A2*RDOT(2)
      VI(3) = XDOT(3,I) + A2*RDOT(3)
      VI(4) = VI(1) - RDOT(1)
      VI(5) = VI(2) - RDOT(2)
      VI(6) = VI(3) - RDOT(3)
*
      RETURN
*
      END
