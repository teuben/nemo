      SUBROUTINE UPDATE(IPAIR)
*
*
*       List modifications after KS termination.
*       ----------------------------------------
*
      INCLUDE 'common4.h'
*
*
*       Adjust perturber lists to new sequence.
      ICM = N + IPAIR
      IF (NPAIRS.EQ.0) GO TO 90
*
*       Replace c.m. by components and reduce subsequent members by one.
      DO 80 JPAIR = 1,NPAIRS
          J = 2*JPAIR - 1
          NNB = LIST(1,J)
          DO 76 L = NNB+1,2,-1
              IF (LIST(L,J).EQ.ICM) THEN
*       Move up all subsequent list members by one and reduce membership.
                  DO 72 K = L,NNB
                      LIST(K,J) = LIST(K+1,J)
   72             CONTINUE
                  LIST(1,J) = LIST(1,J) - 1
                  KCOMP = 2
                  IF (NNB.GT.LMAX-2) KCOMP = 1
                  DO 74 K = NNB,2,-1
                      IF (LIST(K,J).GT.ICOMP) THEN
*       Move members down by two (or one) to make room for components.
                          LIST(K+KCOMP,J) = LIST(K,J)
                      ELSE
*       Replace old c.m. by components and increase membership.
                          LIST(K+1,J) = ICOMP
                          IF (KCOMP.EQ.2) LIST(K+2,J) = JCOMP
                          LIST(1,J) = LIST(1,J) + KCOMP
                          GO TO 80
                      END IF
   74             CONTINUE
*       Include special case of c.m. being the last or only perturber.
                  LIST(2,J) = ICOMP
                  LIST(3,J) = JCOMP
                  LIST(1,J) = NNB + 1 
              END IF
   76     CONTINUE
   80 CONTINUE
*
*       Reduce recent c.m. members by one.
      DO 85 JPAIR = 1,NPAIRS
          J = 2*JPAIR - 1
          NNB1 = LIST(1,J) + 1
          DO 82 L = 2,NNB1
              IF (LIST(L,J).GT.ICM) LIST(L,J) = LIST(L,J) - 1
   82     CONTINUE
   85 CONTINUE
*
*       Copy flag index of disrupted pair (set in KSTERM).
   90 IFLAG = JLIST(3)
*       Add primordial pairs to LISTD (skip new KS pairs).
      IF (IFLAG.EQ.0) GO TO 110
*
*       Check list of disrupted component names.
      NNB = LISTD(1) - 1
      KCOMP = 0
      DO 100 K = 2,NNB+1,2
          IF (LISTD(K).EQ.JLIST(1).AND.LISTD(K+1).EQ.JLIST(2)) KCOMP = 1
  100 CONTINUE
*
*       Include both components if not already members.
      IF (KCOMP.EQ.0) THEN
          IF (NNB.GT.46) THEN
              DO 102 K = 2,NNB
                  LISTD(K) = LISTD(K+2)
  102         CONTINUE
              NNB = NNB - 2
          END IF
*       Add most recent names at the end (limit is 10 pairs).
          LISTD(NNB+3) = JLIST(1)
          LISTD(NNB+4) = JLIST(2)
          LISTD(1) = NNB + 3
      END IF
*
*       Remove first component of old KS pair and c.m. body from NLIST.
  110 I = 2*IPAIR - 1
      CALL NLMOD(I,-1)
      CALL NLMOD(ICM,-1)
*
*       Rename any older single KS components and more recent c.m. bodies.
      IF (IPAIR.LE.NPAIRS) THEN
          NNB1 = NLIST(1) + 1
          DO 120 L = 2,NNB1
*       Reduce index of any subsequent c.m. and first components by 1 & 2.
              IF (NLIST(L).GT.ICM) NLIST(L) = NLIST(L) - 1
              IF (NLIST(L).LE.JCOMP.AND.NLIST(L).GT.2*IPAIR - 1)
     &                                           NLIST(L) = NLIST(L) - 2
  120     CONTINUE
      END IF
*
      RETURN
*
      END
