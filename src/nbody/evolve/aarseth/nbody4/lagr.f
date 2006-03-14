      SUBROUTINE LAGR(C)
*
*
*       Lagrangian radii.
*       -----------------
*
      INCLUDE 'common4.h'
      REAL*8 R2,RHO,MRT,RT
      COMMON/WORK1/ R2(NMAX),RHO(NMAX)
      PARAMETER (LX=11)
      REAL*8 C(3),FLAGR(LX),RLAGR(LX),RM(LX),DENS(LX),VR(LX),AVM(LX)
*     DATA FLAGR/-1.9,-1.7,-1.5,-1.3,-1.1,-.9,-.7,-.5,-.3,-.1/
*     DATA FLAGR/0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,
*    &           0.75,0.9/
      DATA FLAGR/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.625,0.75,0.9/
*
*
*       Set square radii of single particles & c.m. bodies.
      NP = 0
      DO 10 I = IFIRST,NTOT
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
   10 CONTINUE
*
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)

*       Determine Lagrangian radii for specified mass fractions.
      ZM = 0.0
      DO 20 IL = 1,LX
          ZM1 = ZM
          ZM = 0.0
          ZMH = FLAGR(IL)*ZMASS
          VM2 = 0.0
          IF (IL.GT.1) THEN
              R1 = SQRT(R2(I))
              IPREV = I
          ELSE
              R1 = 0.0
              IPREV = 0
          END IF
          I = 0
   15     I = I + 1
          IM = JLIST(I)
          ZM = ZM + BODY(IM)
*       Sum mass-weighted square velocity in current shell.
          IF (ZM.GT.ZM1) THEN
              VM2 = VM2 + BODY(IM)*(XDOT(1,IM)**2 + XDOT(2,IM)**2 +
     &                                              XDOT(3,IM)**2)
          END IF
          IF (ZM.LT.ZMH) GO TO 15
          RLAGR(IL) = SQRT(R2(I))
*       Form mean square velocity and density (print outer shell radii).
          DM = ZM - ZM1
          IF (ABS(DM).LT.1.0D-10) THEN
              VR(IL) = 0.0
          ELSE
              VR(IL) = SQRT(VM2/DM)
          END IF
          DV = 2.0*TWOPI/3.0*(R2(I)**1.5 - R1**3)
          DENS(IL) = DM/DV
          IF (DENS(IL).LE.0.0D0) DENS(IL) = 1.0
          RM(IL) = 0.5*(R1 + RLAGR(IL))
*       Obtain the average mass (units of the mean value).
          IF (I.GT.IPREV) THEN
              AVM(IL) = DM/(FLOAT(I - IPREV)*BODYM)
          ELSE
              AVM(IL) = 0.0
          END IF
   20 CONTINUE
*
*       Obtain half-mass radius separately.
      ZM = 0.0
      ZMH = 0.5*ZMASS
      I = 0
   25 I = I + 1
      IM = JLIST(I)
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZMH) GO TO 25
*
*       Replace approximate half-mass radius by actual value.
      RSCALE = SQRT(R2(I))
*
*       Determine tidal radius and mass if kz(14) = 2.
      if (kz(14).eq.2.and.tidal(1).ne.0.) then
      I = np
      mrt = zmass
   30 continue
      if (mrt.gt.0.) then
         rt = (mrt/tidal(1))**(1./3.)
      else
         rt = 0.0
      endif
      IF (sqrt(r2(i)).gt.rt.and.rt.gt.0.and.i.gt.1) then
         IM = JLIST(I)
         Mrt = Mrt - BODY(IM)
         i = i - 1
         GO TO 30
      endif
      endif

      IF ((KZ(7).EQ.2.OR.KZ(7).EQ.4).AND.TTOT.GE.TPRINT) THEN
          IF (KZ(14).EQ.2) WRITE (6,*)  tphys,mrt,rt
      END IF
*
      IF (KZ(7).GE.3.AND.TTOT.GE.TPRINT) THEN
          IF (KZ(14).EQ.2) WRITE (12,*)  tphys,mrt,rt
      END IF
*
*       Check output options (line printer or unit 7 or both).
      IF (KZ(7).EQ.2.OR.KZ(7).EQ.4.AND.TTOT.GE.TPRINT) THEN
          WRITE (6,40)  (LOG10(RLAGR(K)),K=1,LX)
   40     FORMAT (/,3X,'LAGR:  ',13F7.3)
      END IF
*
      IF (KZ(7).GE.3.AND.TTOT.GE.TPRINT) THEN
          WRITE (12,50)  TTOT, (LOG10(RLAGR(K)),K=1,LX)
   50     FORMAT ('  LAGR:  ',F7.1,13F7.3)
          CALL FLUSH(12)
      END IF
*
      IF (KZ(7).GE.5.AND.TTOT.GE.TPRINT) THEN
          WRITE (26,60)  TTOT, (DENS(K),K=1,LX)
   60     FORMAT ('  DENSITY (T =',F6.1,'): ',1P,13E10.2)
          WRITE (26,65)  (RM(K),K=1,LX)
   65     FORMAT ('  DISTANCE:  ',1P,13E10.2)
          WRITE (27,70)  TTOT, (VR(K),K=1,LX)
   70     FORMAT ('  VELOCITY (T =',F6.1,'): ',1P,13E10.2)
          WRITE (27,65)  (RM(K),K=1,LX)
          WRITE (36,75)  TTOT, (AVM(K),K=1,LX)
   75     FORMAT ('  AVERAGE MASS (T =',F6.1,'): ',13F7.3)
          WRITE (36,65)  (RM(K),K=1,LX)
      END IF
*
      RETURN
*
      END
