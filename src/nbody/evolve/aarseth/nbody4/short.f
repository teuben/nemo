      SUBROUTINE SHORT(LENGTH,NXTLST)
*
*
*       Short time-step list.
*       ---------------------
*
      INCLUDE 'common4.h'
      INTEGER  NXTLST(100)
*
*
*       Decide between modifying list or increasing membership.
    1 NNB = LSHORT(1)
      IF (NNB.GE.8) THEN
          K1 = NNB - 7
          K2 = NNB - (K1 - 1)
          K10 = K1
          IF(NNB.EQ.8) K1 = 2
*       Reduce membership below 8 by removing old entries.
          DO 5 K = K1,K2
              LSHORT(K) = LSHORT(K+K10)
    5     CONTINUE
          NNB = NNB - K10
*       Add new KS candidates and c.m. particles unless already present.
          DO 10 K = 1,LENGTH
              J = NXTLST(K)
              DO 8 L = 2,NNB+1
                  IF (LSHORT(L).EQ.J) GO TO 10
    8         CONTINUE
              IF (STEP(J).LT.SMIN) THEN
                  NNB = NNB + 1
                  LSHORT(NNB+1) = J
              END IF
   10     CONTINUE
          LSHORT(1) = NNB
      ELSE
*       Check existing members before adding current small step particles.
          DO 20 K = 1,LENGTH
              J = NXTLST(K)
              DO 15 L = 2,NNB+1
                  IF (LSHORT(L).EQ.J) GO TO 20
   15         CONTINUE
              IF (STEP(J).LT.SMIN) THEN
                  NNB = NNB + 1
                  LSHORT(NNB+1) = J
*                 LSHORT(1) = NNB
              END IF
   20     CONTINUE
          LSHORT(1) = NNB
      END IF
*
      RETURN
*
      END
