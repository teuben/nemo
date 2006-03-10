      SUBROUTINE NCLIST(I,RS2,RI2,RS2I)
*
*
*       Neighbour list for core particles.
*       ----------------------------------
*
      INCLUDE 'common4.h'
      DATA NBPRV /0/
      SAVE NBPRV
*
*
*       Obtain neighbour list inside RS2 for body #I.
      IF (RI2.LT.0.01.AND.RS2.GT.0.02*RI2) THEN
*         WRITE (7,5)  TTOT,I,SQRT(RS2),SQRT(RI2)
*   5     FORMAT (' NCLIST REDUCE    T I RS RI ',F8.1,I6,1P,2E10.2)
          RS2 = MIN(0.06*RI2,RC**2)
      END IF
      IF (NBPRV.GT.20)  RS2 = 0.5*RS2
   10 CALL NBLIST(I,RS2,NNB)
      RS2I = RS2
      ILIST(1) = NNB
*
*     RS20 = RS2
*       Ensure at least six members.
      IF (NNB.LT.6) THEN
          FAC = 1.0 + 0.1*(6.0 - FLOAT(NNB))
          RS2 = FAC*RS2
*         RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
*    &                                   (X(3,I) - RDENS(3))**2
*     WRITE (33,22)  I,NNB,FAC,SQRT(RS20),SQRT(RS2),SQRT(RI2)
*     WRITE (6,22)  I,NNB,FAC,SQRT(RS20),SQRT(RS2),SQRT(RI2)
*  22 FORMAT (' REPEAT:    I NNB FAC RS0 RS R ',I6,I5,F6.2,2F7.3,F7.3)
          GO TO 10
      END IF
*
*       Modify next neighbour radius according to previous membership.
      IF (NNB.LT.10) THEN
          FAC = 1.0 + 0.1*(10.0 - FLOAT(NNB))
          FAC = MIN(FAC,1.2D0)
          RS2 = FAC*RS2
      ELSE
          FAC = 1.0 - 0.1*(FLOAT(NNB) - 9.0)
          FAC = MAX(FAC,0.5D0)
          RS2 = FAC*RS2
      END IF
      NBPRV = NNB
*
      RETURN
*
      END
