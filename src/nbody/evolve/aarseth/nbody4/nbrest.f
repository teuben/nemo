      SUBROUTINE NBREST(ICM,NSYS,NP)
*
*
*       Restore ghosts in perturber lists.
*       ----------------------------------
*
      INCLUDE 'common4.h'
*
*
*       Examine all KS perturber lists.
      DO 100 JPAIR = 1,NP
          I = 2*JPAIR - 1
          NNB1 = LIST(1,I) + 1
          IF (NNB1.LE.1) GO TO 100
*
*       First see whether body #ICM is a perturber.
          DO 20 L = 2,NNB1
              IF (LIST(L,I).EQ.ICM) THEN
                  GO TO 30
              ELSE
                  IF (LIST(L,I).GT.ICM) GO TO 100
              END IF
   20     CONTINUE
*
*       Add other members of the subsystem.
   30     DO 60 K = 1,NSYS
              J = JLIST(K)
              IF (J.EQ.ICM) GO TO 60
*
*       Skip addition if body #J is already a member.
              DO 40 L = 2,NNB1
                  IF (LIST(L,I).EQ.J) GO TO 60
   40         CONTINUE
*
*       Move members down until correct location identified.
              DO 50 L = NNB1,1,-1
                  IF (LIST(L,I).GT.J.AND.L.GT.1) THEN
                      LIST(L+1,I) = LIST(L,I)
                  ELSE
*       Place body #J in vacated location and increase membership.
                      LIST(L+1,I) = J
                      LIST(1,I) = LIST(1,I) + 1
                      NNB1 = NNB1 + 1
                      GO TO 60
                  END IF
   50         CONTINUE
   60     CONTINUE
  100 CONTINUE
*
      RETURN
*
      END
