      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
      INCLUDE 'params.h'
      PARAMETER  (NA=64*NMAX,NB=(LMAX+2)*NMAX+29)
*       ----------------------------------------------
*       Alternative declaration for full *8 version.
*     PARAMETER  (NA=108*NMAX,NB=(LMAX+2)*NMAX+29)
*       ----------------------------------------------
      REAL*4  A,D
      INTEGER IB,IC
*
      COMMON/NBODY/  A(NA)
      COMMON/NAMES/  IB(NB)
      COMMON/COUNTS/ IC(13)
      COMMON/PARAMS/ D(66)
*       Alternative declaration for full *8 version.
*     COMMON/PARAMS/ D(130)
*
*
*       Open unit #J by reading dummy and rewinding.
      REWIND J
      READ (J,ERR=10,END=10)  DUMMY
   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
          READ (J)   A, IB, IC, D
      ELSE
          WRITE (J)  A, IB, IC, D
          END FILE J
          CLOSE (UNIT=J)
      END IF
*
      RETURN
*
      END
