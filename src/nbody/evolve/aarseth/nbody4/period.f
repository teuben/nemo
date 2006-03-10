      SUBROUTINE PERIOD
*
*
*       Distribution of periods & semi-major axes.
*       ------------------------------------------
*
      INCLUDE 'common4.h'
      INTEGER  NP(30),NA(30)
*
*
      DO 1 K = 1,30
          NP(K) = 0
          NA(K) = 0
    1 CONTINUE
*
      NB = 0
      A0 = 4.0/FLOAT(NZERO)
      P0 = DAYS*A0*SQRT(A0/(2.0*BODYM))
*
      II = 0
*       Form histograms for all KS pairs.
      DO 10 JP = 1,NPAIRS
          IF (H(JP).GT.0.0.OR.BODY(N+JP).LE.0.0D0) GO TO 10
          SEMI = -0.5*BODY(N+JP)/H(JP)
          P = DAYS*SEMI*SQRT(SEMI/BODY(N+JP))
          IF (SEMI.LT.A0) THEN
              K = 2 + LOG10(A0/SEMI)/LOG10(2.0)
          ELSE
              K = 1
          END IF
          IF (P.LT.P0) THEN
              L = 2 + LOG10(P0/P)/LOG10(2.0)
          ELSE
              L = 1
          END IF
          II = II + 1
          NA(K) = NA(K) + 1
          NP(L) = NP(L) + 1
   10 CONTINUE
*
*       Search original binary components.
      DO 30 I = IFIRST,N
          IF (NAME(I).GT.2*NBIN0) GO TO 30
          DO 20 J = I+1,N
              IF (IABS(NAME(I) - NAME(J)).GT.1) GO TO 20
              RIJ2 = 0.0
              VIJ2 = 0.0
              DO 15 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                  VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
   15         CONTINUE
              SEMI = 2.0/SQRT(RIJ2) - VIJ2/(BODY(I) + BODY(J))
              IF (SEMI.LE.0.0) GO TO 30
              SEMI = 1.0/SEMI
              P = DAYS*SEMI*SQRT(SEMI/(BODY(I) + BODY(J)))
              IF (SEMI.LT.A0) THEN
                  K = 2 + LOG10(A0/SEMI)/LOG10(2.0)
              ELSE
                  K = 1
              END IF
              IF (P.LT.P0) THEN
                  L = 2 + LOG10(P0/P)/LOG10(2.0)
              ELSE
                  L = 1
              END IF
              II = II + 1
              NA(K) = NA(K) + 1
              NP(L) = NP(L) + 1
   20     CONTINUE
   30 CONTINUE
*
      WRITE(47,*)TTOT,N,NPAIRS,II
      WRITE(47,*)A0,P0
      DO 40 I = 1,30
          WRITE(47,*)I,NA(I),NP(I)
   40 CONTINUE
      CALL FLUSH(47)
*
      RETURN
*
      END
