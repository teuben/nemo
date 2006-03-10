      SUBROUTINE GHOST(J)
*
*
*       Initialization of ghost particle.
*       ---------------------------------
*
      INCLUDE 'common4.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX),LISTCM(LMAX)
*
*
*       Form ghost and initialize integration variables for body #J.
      BODY(J) = 0.0D0
      T0(J) = 1.0E+06
      DO 10 K = 1,3
          X0DOT(K,J) = 0.0D0
          XDOT(K,J) = 0.0D0
          F(K,J) = 0.0D0
          FDOT(K,J) = 0.0D0
          D2(K,J) = 0.0D0
          D3(K,J) = 0.0D0
   10 CONTINUE
*
*       Set large X0 & X to avoid perturber selection (no escape).
      X0(1,J) = 1.0D+06
      X(1,J) = 1.0D+06
*
*       Remove ghost from KS perturber lists containing body #ICH.
      JLIST(1) = J
      CALL NBREM(ICH,1,NPAIRS)
*
      RETURN
*
      END
