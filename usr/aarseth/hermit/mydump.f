      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
      PARAMETER  (NMAX=10,NA=61*NMAX+20,NP=2*NMAX+84)
      REAL*4  A,E,P
      INTEGER  IR
*
*
      COMMON/NBODY/  A(NA)
      COMMON/PARAMS/ E(60)
      COMMON/RAND2/  IR(99)
      COMMON/BLOCKS/ P(NP)
*
*
*       Open unit #J by reading dummy and rewinding.
      REWIND J
      READ (J,ERR=10,END=10)  DUMMY
   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
          READ (J)   A, E, IR, P
      ELSE
          WRITE (J)  A, E, IR, P
          END FILE J
          CLOSE (UNIT=J)
      END IF
*
      RETURN
*
      END
