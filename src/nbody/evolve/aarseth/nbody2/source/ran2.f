      FUNCTION RAN2(IDUM)
*
*
*       Random number generator (Press p. 195).
*       ---------------------------------------
*
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1./M)
      COMMON/RAND2/  IY,IFF,IR(97) 
*     DATA  IFF /0/
*
*
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
          IFF = 1
          IDUM = MOD(IC-IDUM,M)
          DO 11 J = 1,97
              IDUM = MOD(IA*IDUM+IC,M)
              IR(J) = IDUM
   11     CONTINUE
          IDUM = MOD(IA*IDUM+IC,M)
          IY = IDUM
      END IF
      J = 1 + (97*IY)/M
      IF (J.GT.97.OR.J.LT.1) WRITE (6,12)  J, IDUM
   12 FORMAT (/,'  TROUBLES IN RAN2   J IDUM ',2I12)
      IY = IR(J)
      RAN2 = IY*RM
      IDUM = MOD(IA*IDUM+IC,M)
      IR(J) = IDUM
*
      RETURN
*
      END 
