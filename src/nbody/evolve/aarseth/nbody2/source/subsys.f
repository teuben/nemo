      SUBROUTINE SUBSYS
*
*
*       Initial subsystems.
*       -------------------
*
      INCLUDE 'common2.h'
*
*
*       Read c.m. displacement & eccentricity of relative motion.
      READ (5,*)  XCM, ECC
*
*       Set relative velocity for two clumps of ZMASS at distance 2*XCM.
      VREL = SQRT(ZMASS*(1.0 - ECC)/XCM)
*
*       Copy second system from the first and separate the two clumps.
      DO 10 I = 1,N
          J = I + N
          BODY(J) = BODY(I)
          DO 5 K = 1,3
              X(K,J) = X(K,I)
              XDOT(K,J) = XDOT(K,I)
    5     CONTINUE
*       Displace each clump by XCM and introduce apocentre y-velocity.
          X(1,I) = X(1,I) + XCM
          X(1,J) = X(1,J) - XCM
          XDOT(2,I) = XDOT(2,I) + 0.5*VREL
          XDOT(2,J) = XDOT(2,J) - 0.5*VREL
   10 CONTINUE
*
*       Set new total particle number (NB! Do not exceed NMAX).
      N = 2*N
      
      WRITE (6,20)  XCM, ECC, VREL
   20 FORMAT (/,12X,'SUBSYSTEMS:    XCM =',F5.2,'  ECC =',F5.2,
     &                                          '  VREL =',F5.2)
*
      RETURN
*
      END
