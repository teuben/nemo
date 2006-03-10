      SUBROUTINE KCPERT(I,I1,FIRR,FD)
*
*
*       Perturbation on KS components due to chain.
*       -------------------------------------------
*
      INCLUDE 'common4.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX),LISTCM(LMAX)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FD(3),FP(3),FPD(3),FPS(3),FDS(3)
*
*
*       Obtain force & derivative on each KS component (saving the first).
      J1 = I1
      DO 10 KCOMP = 1,2
          DO 5 K = 1,3
              FP(K) = 0.0D0
              FPD(K) = 0.0D0
              XI(K) = X(K,J1)
              XIDOT(K) = XDOT(K,J1)
    5     CONTINUE
          CALL FCHAIN(J1,XI,XIDOT,FP,FPD)
          IF (KCOMP.EQ.1) THEN
              DO 8 K = 1,3
                  FPS(K) = FP(K)
                  FDS(K) = FPD(K)
    8         CONTINUE
          END IF
          J1 = J1 + 1
   10 CONTINUE
*
*       Add mass-weighted contributions to irregular force & derivative.
      BODYIN = 1.0/BODY(I)
      DO 15 K = 1,3
          FIRR(K) = FIRR(K) + (BODY(I1)*FPS(K) + BODY(J1)*FP(K))*BODYIN
          FD(K) = FD(K) + (BODY(I1)*FDS(K) + BODY(J1)*FPD(K))*BODYIN
   15 CONTINUE
*
      RETURN
*
      END
