      SUBROUTINE FRAME 
*
*
*       Movie frame output.
*       -------------------
*
      INCLUDE 'common1.h'
*
*
*       Specify colour type (unless individual value).
      IC = 1
*
*       Begin each frame by the current particle number.
      WRITE (7,1)  N
    1 FORMAT (I5)
*
*       Produce data file of coordinates & colour type.
      DO 10 I = 1,N
          WRITE (7,5)  (X(K,I),K=1,3),IC
    5     FORMAT (3F7.2,I2)
   10 CONTINUE
*
*       Update time for next frame output.
      TFRAME = TFRAME + DELTAF
*
      RETURN
*
      END
