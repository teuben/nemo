      SUBROUTINE BINARY(IESC)
*
*
*       Binary determination.
*       ---------------------
*
      IMPLICIT  REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  Q(8),P(8),R,R1,R2,ENERGY,M(3),X(3,3),XDOT(3,3),
     &               RCOLL,ERROR,C11,C12,C19,C20,C24,C25,NSTEPS,NAME(3)
      COMMON/TBIN0/  RGRAV,EBIN0,ECC0,DB,I10,I20
*
*
*       Set index of second binary component & escaper.
      IMIN = 1
      IF (R2.LT.R1) IMIN = 2
      I = 3 - IMIN
*
*       Evaluate orbital elements.
      RDOT = 0.0D0
      VREL2 = 0.0D0
      RI2 = 0.0
      RIDOT = 0.0
      DO 10 K = 1,3
          RDOT = RDOT + (X(K,3) - X(K,IMIN))*(XDOT(K,3) - XDOT(K,IMIN))
          VREL2 = VREL2 + (XDOT(K,3) - XDOT(K,IMIN))**2
          RI2 = RI2 + X(K,I)**2
          RIDOT = RIDOT + X(K,I)*XDOT(K,I)
   10 CONTINUE
*
      RB = MIN(R1,R2)
      MB = M(3) + M(IMIN)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      E = SQRT((1.0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      EB = -0.5D0*M(3)*M(IMIN)/SEMI
*
*       Save initial elements on first call.
      IF (IESC.EQ.0) THEN
          EBIN0 = MIN(EB,0.0)
          ECC0 = E
          I10 = NAME(IMIN)
          I20 = NAME(3)
          IESC = I
          GO TO 30
      END IF
*
*       Set distance & radial velocity of body #I with respect to binary.
      RI = SQRT(RI2)
      ZMT = M(1) + M(2) + M(3)
      FAC = ZMT/MB
      RIDOT = FAC*RIDOT/RI
      RI = FAC*RI
*
*       Evaluate the escape criterion due to Standish (Celes. Mech. 4, 44).
      RATIO = RGRAV/(MB*RI)
      VC2 = 2.0*ZMT*(1.0/RI + M(3)*M(IMIN)*RATIO**2/(RI - RGRAV))
      DV2 = RIDOT**2 - VC2
      DB = (EB - EBIN0)/(0.5*ENERGY)
      ZMU = M(I)*MB/ZMT
      EFAC = (ENERGY - 2.0*EB)/ZMU
*
*       Express terminal velocity in units of mean binary velocity.
      IF (EFAC.GT.0.0) THEN
          VESC = SQRT(EFAC*SEMI/MB)
      ELSE
          VESC = RIDOT
      END IF
*
*       Check for diagnostic output of dominant binary.
      IF (EB.LT.0.5*ENERGY) THEN
          EBE = EB/(0.5*ENERGY)
*       Binding energy scaled by total energy.
*
          WRITE (6,20)  NAME(IMIN), NAME(3), R1, R2, SEMI, E, EB,
     &                  EBE, VESC, DV2, DB
   20     FORMAT  (/,5X,'BINARY ',2I3,'  R1 =',F5.1,'  R2 =',F5.1,
     &                  '  A =',F7.3,'  E =',F6.3,'  EB =',F8.2,
     &                  '  EB/E =',F7.3,'  VESC =',F5.2,'  DV2 =',F6.2,
     &                  '  DB =',F6.2)
      END IF
*
      IESC = I
   30 RETURN
*
      END
