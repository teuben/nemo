*
      SUBROUTINE ZERO
*
*
*       Initialization of global scalars. 
*       ---------------------------------
*
      INCLUDE 'common1.h'
*
*
*       Initialize parameters & counters and set useful constants.
      TIME = 0.0D0
      TNEXT = 0.0
      TLIST = 0.0
      TFRAME = 0.0
      CPUTOT = 0.0
      ERRTOT = 0.0
      NPRINT = 0
      MODEL = 0
      NDUMP = 0
      NTIMER = 0
      NSTEPI = 0
*
      ONE3 = 1.0/3.0D0
      ONE6 = 1.0/6.0D0
      ONE9 = 1.0/9.0D0
      ONE12 = 1.0/12.0D0
*
      RETURN
*
      END
