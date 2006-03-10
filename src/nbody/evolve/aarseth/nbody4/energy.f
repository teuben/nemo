      SUBROUTINE ENERGY
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'common4.h'
*
*
*       Sum the total energy of regularized pairs.
      EBIN = 0.0D0
      DO 10 IPAIR = 1,NPAIRS
*       Skip pairs with zero mass of c.m. particle (merged binary ghost).
          IF (BODY(N+IPAIR).GT.0.0D0) THEN
*       Predict coordinates, velocities & binding energy.
              CALL RESOLV(IPAIR,1)
              EBIN = EBIN + BODY(2*IPAIR-1)*BODY(2*IPAIR)*HT/
     &                                                     BODY(N+IPAIR)
          END IF
   10 CONTINUE
*
*       Calculate the potential energy on GRAPE (or host).
      IF(TIME.LE.0.D0)THEN
         CALL POTN0(SUMP)
      ELSE
         CALL POTN2(SUMP)
      ENDIF
      POT = SUMP
*
*       Sum the kinetic energy (include c.m. bodies but not components).
      ZKIN = 0.0D0
      DO 40 I = IFIRST,NTOT
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
   40 CONTINUE
      ZKIN = 0.5D0*ZKIN
*
*       Obtain the tidal potential if linearized external field is present. 
      IF (KZ(14).EQ.1.OR.KZ(14).EQ.2) THEN
          CALL XTRNLV(IFIRST,NTOT)
      END IF
*
*       Check differential potential energy due to chain subsystem.
      IF (NCH.GT.0) THEN
          CALL CHPOT(DP)
          POT = POT + DP
      END IF
*
*       Energy = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL + ECH.
*
      RETURN
*
      END
