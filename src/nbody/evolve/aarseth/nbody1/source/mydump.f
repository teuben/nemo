      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
      INCLUDE 'params.h'
      PARAMETER  (NA=43*NMAX,NB=2*NMAX+24)
*       --------------------------------------------
*       Alternative declaration for full *8 version.
*     PARAMETER  (NA=66*NMAX,NB=2*NMAX+24)
*       --------------------------------------------
      REAL*4  A,C
      INTEGER IB
*
      COMMON/NBODY/  A(NA)
      COMMON/NAMES/  IB(NB)
      COMMON/PARAMS/ C(43)
*       Alternative declaration for full *8 version.
*     COMMON/PARAMS/ C(84)
*
*
*       Open unit #J by reading dummy and rewinding.
      REWIND J
      READ (J,ERR=10,END=10)  DUMMY
   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
          READ (J)   A, IB, C
      ELSE
          WRITE (J)  A, IB, C
          END FILE J
          CLOSE (UNIT=J)
      END IF
*
      RETURN
*
      END
