      SUBROUTINE PIPES(NXTLEN)
*
*
*       Diagnostics of pipe activity on GRAPE.
*       ---------------------------------------
*
      INCLUDE 'common4.h'
*
*
*       Count number of blocks and sum block size.
      IPIPE(1) = IPIPE(1) + 1
      IPIPE(2) = IPIPE(2) + NXTLEN
*
*       Include special count of all pipes and just 1 or 2.
      IF (NXTLEN.GT.NPIPE) THEN
          IPIPE(3) = IPIPE(3) + 1
      ELSE IF (NXTLEN.EQ.1) THEN
          IPIPE(4) = IPIPE(4) + 1
      ELSE IF (NXTLEN.EQ.2) THEN
          IPIPE(5) = IPIPE(5) + 1
      END IF
*
*       Form histogram of active pipes.
      FAC = FLOAT(NXTLEN)/FLOAT(NPIPE) + 1.0E-08
      IFAC = NPIPE*(FAC - INT(FAC))
      IF (IFAC.EQ.0) IFAC = NPIPE
      IFAC = MIN(IFAC,15)
      IPIPE(IFAC+5) = IPIPE(IFAC+5) + 1
*
      RETURN
*
      END
