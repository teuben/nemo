      SUBROUTINE REMOVE(I,KCASE)
*
*
*       Removal of particle.
*       --------------------
*
      INCLUDE 'common2.h'
      REAL*4  A(6)
*
*
*       Remove escaper or merger component (KCASE = 1, 2).
      IF (KCASE.EQ.2) GO TO 10
*
*       Correct force & first derivative of neighbours (only for escape).
      NNB1 = LIST(1,I) + 1
      DO 5 L = 2,NNB1
          J = LIST(L,I)
          RIJ2 = 0.0
          RIJDOT = 0.0
*
          DO 1 K = 1,3
              A(K) = X(K,I) - X(K,J)
              A(K+3) = XDOT(K,I) - XDOT(K,J)
              RIJ2 = RIJ2 + A(K)**2
              RIJDOT = RIJDOT + A(K)*A(K+3)
    1     CONTINUE
*
          RIJ2 = RIJ2 + EPS2
          A8 = BODY(I)/(RIJ2*SQRT(RIJ2))
          A9 = 3.0*RIJDOT/RIJ2
*
          DO 2 K = 1,3
              DF = A(K)*A8
              DD = (A(K+3) - A9*A(K))*A8
              F(K,J) = F(K,J) - 0.5*DF
              FI(K,J) = FI(K,J) - DF
              FDOT(K,J) = FDOT(K,J) - ONE6*DD
              D1(K,J) = D1(K,J) - DD
    2     CONTINUE
    5 CONTINUE
*
*       Skip table updates for last body.
   10 IF (I.GT.N) GO TO 60
*
*       Move up all COMMON variables (escaper or merged body).
      DO 15 J = I,N
          J1 = J + 1
          DO 12 K = 1,3
              X(K,J) = X(K,J1)
              X0(K,J) = X0(K,J1)
              X0DOT(K,J) = X0DOT(K,J1)
              XDOT(K,J) = XDOT(K,J1)
              F(K,J) = F(K,J1)
              FDOT(K,J) = FDOT(K,J1)
              FI(K,J) = FI(K,J1)
              D1(K,J) = D1(K,J1)
              D2(K,J) = D2(K,J1)
              D3(K,J) = D3(K,J1)
              FR(K,J) = FR(K,J1)
              D1R(K,J) = D1R(K,J1)
              D2R(K,J) = D2R(K,J1)
              D3R(K,J) = D3R(K,J1)
   12     CONTINUE
*
          BODY(J) = BODY(J1)
          RS(J) = RS(J1)
          NAME(J) = NAME(J1)
          STEP(J) = STEP(J1)
          STEPR(J) = STEPR(J1)
          T0(J) = T0(J1)
          T1(J) = T1(J1)
          T2(J) = T2(J1)
          T3(J) = T3(J1)
          T0R(J) = T0R(J1)
          T1R(J) = T1R(J1)
          T2R(J) = T2R(J1)
          T3R(J) = T3R(J1)
*       Transfer unmodified neighbour list.
          NNB = LIST(1,J1) + 1
          DO 14 L = 1,NNB
              LIST(L,J) = LIST(L,J1)
   14     CONTINUE
   15 CONTINUE
*
*       Delete body #I from neighbour lists and reduce higher locations.
   60 DO 150 J = 1,N
          NNB = LIST(1,J)
          L = 2
   70     IF (LIST(L,J).NE.I) GO TO 130
*
*       Move up the remaining list members and reduce membership by one.
          DO 80 K = L,NNB
              LIST(K,J) = LIST(K+1,J)
   80     CONTINUE
          LIST(1,J) = LIST(1,J) - 1
*
*       Reduce both time-steps in order to minimize the error effect.
          A1 = T0(J) + STEP(J)
          STEP(J) = STEP(J) - 0.5*(A1 - TIME)
          STEPR(J) = 0.5*STEPR(J)
*
*       Add body #J to time-step list unless already a member. 
          IF (A1.GT.TLIST.AND.T0(J) + STEP(J).LT.TLIST) THEN
              NNB1 = NLIST(1) + 1
              NLIST(NNB1+1) = J
              NLIST(1) = NNB1
          END IF
*
          IF (LIST(1,J).EQ.0) THEN
*       Add a spurious body as neighbour if list only contains body #I.
              K = N
              IF (K.EQ.J) K = N - 1
              LIST(1,J) = 1
              LIST(2,J) = K
*
*       Initialize neighbour force & higher differences.
              DO 120 K = 1,3
                  FI(K,J) = 0.0
                  D1(K,J) = 0.0
                  D2(K,J) = 0.0
                  D3(K,J) = 0.0
  120         CONTINUE
              GO TO 150
          END IF
*
*       Reduce higher particle locations by one.
  130     IF (LIST(L,J).GT.I) LIST(L,J) = LIST(L,J) - 1
          L = L + 1
          IF (L.LE.LIST(1,J) + 1) GO TO 70
  150 CONTINUE
*
*       Modify the time-step list due to removal of body #I.
      NNB = NLIST(1)
      L = 2
  160 IF (NLIST(L).EQ.I) THEN
*       Move up the subsequent list members and reduce membership by one.
          DO 170 K = L,NNB
              NLIST(K) = NLIST(K+1)
  170     CONTINUE
          NLIST(1) = NLIST(1) - 1
      END IF
*
*       Reduce higher particle locations by one.
      IF (NLIST(L).GT.I) NLIST(L) = NLIST(L) - 1
      L = L + 1
      IF (L.LE.NLIST(1) + 1) GO TO 160
*
      RETURN
*
      END
