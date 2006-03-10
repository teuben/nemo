      SUBROUTINE JACOBI(NTESC)
*
*
*       Jacobi escape criterion.
*       ------------------------
*
      INCLUDE 'common4.h'
*
*
*       Specify escape energy (tidal field or isolated system).
      IF (KZ(14).GT.0) THEN
          ECRIT = -1.5*(TIDAL(1)*ZMASS**2)**0.333
      ELSE
          ECRIT = 0.0
      END IF
*
*       Count all escapers.
      NTESC = 0
      DO 50 I = IFIRST,NTOT
*         RI2 = 0.0
          VI2 = 0.0
          DO 30 K = 1,3
*             RI2 = RI2 + (X(K,I) - RDENS(K))**2
              VI2 = VI2 + XDOT(K,I)**2
   30     CONTINUE
          EI = 0.5*VI2 + PHI(I)
          IF (EI.GT.ECRIT) NTESC = NTESC + 1
   50 CONTINUE
*
*     NS = 0
*     DO 55 I = IFIRST,N
*         IF (BODY(I).GT.0.0D0) NS = NS + 1
*  55 CONTINUE
*     NS = NS - NSUB
*     FRAC = FLOAT(NTESC)/FLOAT(NTOT-IFIRST+1)
*     WRITE (88,60) TTOT, ZMASS, FRAC, NS, NTESC
*  60 FORMAT (' ',F10.3,2F9.4,2I6)
      WRITE (88,65)  TPHYS, N-NPAIRS, NPAIRS, NTESC
   65 FORMAT (' ',F8.1,3I6)
      CALL FLUSH(88)
*
      RETURN
*
      END
