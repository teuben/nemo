      
      SUBROUTINE NLMOD(I,KMOD)
*
*
*       Time-step list modification.
*       ----------------------------
*
      INCLUDE 'common4.h'
*
*
*       Distinguish between addition or removal (KMOD > 0 or < 0).
      NL1 = NLIST(1) + 1
      IF (KMOD.GT.0) THEN
*
*       First see whether body #I is already a member.
          DO 10 L = 2,NL1
              IF (NLIST(L).EQ.I) GO TO 50
   10     CONTINUE
*
*       Add body #I at the end and increase membership.
          NLIST(NL1+1) = I
          NLIST(1) = NLIST(1) + 1
      ELSE
*
*       Search all members since NLIST is not always sequential.
          DO 30 L = 2,NL1
              IF (NLIST(L).EQ.I) THEN
*       Move up subsequent list members and reduce membership.
                  DO 20 K = L,NL1-1
                      NLIST(K) = NLIST(K+1)
   20             CONTINUE
                  NLIST(1) = NLIST(1) - 1
*
*       Ensure membership > 0.
                  IF (NLIST(1).EQ.0) THEN
                      NLIST(1) = 1
                      NLIST(2) = 1
                  END IF
              END IF
*
*       Check reduction of higher locations by one (KMOD = -2: from ESCAPE).
              IF (KMOD.LE.-2.AND.NLIST(L).GT.I) THEN
                  NLIST(L) = NLIST(L) - 1
              END IF
   30     CONTINUE
      END IF
*
   50 RETURN
*
      END
