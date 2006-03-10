      SUBROUTINE RENAME
*
*
*       Renaming of list arrays.
*       ------------------------
*
      INCLUDE 'common4.h'
*
*
*       Remove components from perturber lists and replace by new c.m.
      IF (NPAIRS.EQ.1) GO TO 70
*
      DO 20 JPAIR = 1,NPAIRS-1
*       Only consider first list for each KS pair (note ICOMP < JCOMP).
          J = 2*JPAIR - 1
*       Skip modification if first perturber comes after second component.
          IF (LIST(2,J).GT.JCOMP) GO TO 20
          KCOMP = 0
          NNB = LIST(1,J)
          L = 1
   10     L = L + 1
          IF (LIST(L,J).EQ.ICOMP.OR.LIST(L,J).EQ.JCOMP) THEN
              KCOMP = KCOMP + 1
*       Remove component (unless last) and reduce membership.
              DO 15 K = L,NNB
                  LIST(K,J) = LIST(K+1,J)
   15         CONTINUE
              NNB = NNB - 1
*       Consider same location again unless last member.
              L = L - 1
          END IF
          IF (L.LE.NNB) GO TO 10
          IF (KCOMP.GT.0) THEN
*       Place c.m. last to ensure sequential ordering and set membership.
              LIST(NNB+2,J) = NTOT
              LIST(1,J) = NNB + 1
          END IF
   20 CONTINUE
*
*       Rename any exchanged particles using sequential ordering.
      DO 50 JPAIR = 1,NPAIRS-1
          J = 2*JPAIR - 1
          NNB1 = LIST(1,J) + 1
          KCOMP = 0
*       Check whether two first perturbers occupy new KS locations.
          NNB2 = MIN(NNB1,3)
          DO 30 L = 2,NNB2
              IF (LIST(L,J).EQ.2*NPAIRS-1) THEN
                  KCOMP = KCOMP + 1
                  JPERT(KCOMP) = ICOMP
*       Include case of ICOMP occupying new location of JCOMP (one exchange).
                  IF (ICOMP.EQ.2*NPAIRS) THEN
                      JPERT(KCOMP) = JCOMP
                  END IF
              ELSE IF (LIST(L,J).EQ.2*NPAIRS) THEN
                  KCOMP = KCOMP + 1
                  JPERT(KCOMP) = JCOMP
              END IF
   30     CONTINUE
*
*       Move members up to create new location for exchanged component.
          DO 40 K = 1,KCOMP
              DO 35 L = 2,NNB1
                  IF (L.LT.NNB1.AND.LIST(L+1,J).LT.JPERT(K)) THEN
                      LIST(L,J) = LIST(L+1,J)
                  ELSE
*       Place switched component in correct sequential location and quit.
                      LIST(L,J) = JPERT(K)
                      GO TO 40
                  END IF
   35         CONTINUE
   40     CONTINUE
   50 CONTINUE
*
*       Check removal of close encounter components from NLIST.
   70 CALL NLMOD(ICOMP,-1)
      CALL NLMOD(JCOMP,-1)
*
*       Rename exchanged particles (two stages required for special case).
      DO 80 KCOMP = 1,2
          NNB1 = NLIST(1) + 1
          DO 75 L = 2,NNB1
              IF (NLIST(L).EQ.2*NPAIRS + KCOMP - 2) THEN
                  IF (KCOMP.EQ.1) NLIST(L) = ICOMP
                  IF (KCOMP.EQ.2) NLIST(L) = JCOMP
              END IF
   75     CONTINUE
   80 CONTINUE
*
      RETURN
*
      END
