      SUBROUTINE SHOCK
*
*
*       Disk shocking.
*       --------------
*
      INCLUDE 'common4.h'
      REAL*8  DV(3)
*
*
*       Obtain the velocity change due to tidal shocking.
      ZK1 = 0.0
      ZK2 = 0.0
      zkr1 = 0.0
      zkr2 = 0.0
      DO 10 I = IFIRST,NTOT
          CALL IMPULS(I,DV)
*       Set T0 and X0 for GRAPE predictions and add velocity perturbation.
          T0(I) = TIME
          rpa=(x(1,i)-rdens(1))**2+(x(2,i)-rdens(2))**2+
     &     (x(3,i)-rdens(3))**2
          DO 5 K = 1,3
             if(rpa.lt.rtide**2) then
                zkr1=zkr1+body(i)*xdot(k,i)**2
              end if
              ZK1 = ZK1 + BODY(I)*XDOT(K,I)**2
              XDOT(K,I) = XDOT(K,I) + DV(K)
              X0DOT(K,I) = XDOT(K,I)
              X0(K,I) = X(K,I)
              ZK2 = ZK2 + BODY(I)*XDOT(K,I)**2
             if(rpa.lt.rtide**2) then
                zkr2=zkr2+body(i)*xdot(k,i)**2
              end if
    5     CONTINUE
   10 CONTINUE
*
*       Update the total energy.
      BE(3) = BE(3) + 0.5*(ZK2 - ZK1)
      WRITE (6,15)  TPHYS, NSHOCK, 0.5*(ZK2 - ZK1),
     &  0.5*(zkr2-zkr1),BE(3)
   15 FORMAT (' SHOCK:    TPHYS NSH DE E   ',F8.1,I5,3F12.6)
*
*       Predict all coordinates and new velocities of KS pairs.
      DO 20 J = 1,NPAIRS
          CALL KSRES2(J,J1,J2,0.0D0)
   20 CONTINUE
*
*       Initialize all force polynomials and time-steps.
      DO 30 I = IFIRST,NTOT
*       Specify indicator -2 for skipping GPSEND for all I <= N except first.
          IF (I.GT.IFIRST.AND.I.LE.N) THEN
              IPHASE = -2
          ELSE
              IPHASE = 0
          END IF
          CALL FPOLYI(I)
   30 CONTINUE
*
*       Apply c.m. corrections and update the total energy.
      CALL CMCORR
*
*       Update time for next disk shocking and increase counter.
      TSHOCK = TSHOCK + DTSHOCK
      NSHOCK = NSHOCK + 1
*
*       Specify IPHASE < 0 to ensure new time-step sorting.
      IPHASE = -1
*
      RETURN
*
      END
