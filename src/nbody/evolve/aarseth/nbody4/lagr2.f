      SUBROUTINE LAGR2(C)
*
*
*       Mass distribution for two mass groups.
*       --------------------------------------
*
      INCLUDE 'common4.h'
      REAL*8 R2,RHO,MRT,RT
      COMMON/WORK1/ R2(NMAX),RHO(NMAX)
      PARAMETER (LX=11)
      REAL*8 C(3),FLAGR(LX),RLAGR(LX),RM(LX),DENS(LX),VR(LX)
*     DATA FLAGR/-1.9,-1.7,-1.5,-1.3,-1.1,-.9,-.7,-.5,-.3,-.1/
*     DATA FLAGR/0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,
*    &           0.75,0.9/
      DATA FLAGR/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.625,0.75,0.9/
*
*
*       Set square radii of single particles & c.m. bodies (NAME <= NZERO/5).
      NAM1 = NZERO/5
*       Apply the mass test in same proportion for two Plummer spheres.
      IF (KZ(5).EQ.2) THEN
          NAM1 = (NZERO - N1)/5
          NAM2 = N1 + (NZERO - N1)/5
      END IF
      ITER = 0
*
    1 NP = 0
      ZM1 = 0.0
      DO 10 I = IFIRST,NTOT
          IF (KZ(5).EQ.1) THEN
              IF (NAME(I).GT.NAM1.AND.NAME(I).LE.NZERO) GO TO 10
          ELSE IF (KZ(5).EQ.2) THEN
              IF ((NAME(I).GT.NAM1.AND.NAME(I).LE.N1).OR.
     &            (NAME(I).GT.NAM2.AND.NAME(I).LE.NZERO)) GO TO 10
          END IF
          ZM1 = ZM1 + BODY(I)
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
   10 CONTINUE
*
*       Improve choice of NAM1 using relative deviation (max 10 tries).
      DM = (ZM1 - 0.5*ZMASS)/ZMASS
      IF (ABS(DM).GT.0.001.AND.ITER.LE.10) THEN
          ITER = ITER + 1
          IF (ABS(DM).GT.0.05) DM = DM*0.05/ABS(DM)
          NM = DM*FLOAT(N-NPAIRS)
          IF (NM.NE.0) THEN
              NAM1 = NAM1 - NM
              IF (KZ(5).EQ.2) THEN
                  NAM2 = NAM2 - NM*(NZERO - N1)/(5*NAM2)
              END IF
              GO TO 1
          END IF
      END IF
*
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)

*       Determine Lagrangian radii for specified mass fractions.
      DO 20 IL = 1,LX
          ZM = 0.0
          ZMH = FLAGR(IL)*ZM1
          I = 0
   15     I = I + 1
          IM = JLIST(I)
          ZM = ZM + BODY(IM)
          IF (ZM.LT.ZMH) GO TO 15
          RLAGR(IL) = SQRT(R2(I))
   20 CONTINUE
*
*       Obtain half-mass radius separately.
      ZM = 0.0
      ZMH = 0.5*ZM1
      I = 0
   25 I = I + 1
      IM = JLIST(I)
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZMH) GO TO 25

      RH1 = SQRT(R2(I))
      NP1 = NP
*
      WRITE (31,30)  TIME+TOFF, (LOG10(RLAGR(K)),K=1,LX)
   30 FORMAT (' ',F8.1,13F7.3)
      CALL FLUSH(31)
*
*       Treat the low-mass particles in the same way (exclude cm bodies).
      NP = 0
      ZM2 = 0.0
      DO 50 I = IFIRST,N
          IF (KZ(5).EQ.1) THEN
              IF (NAME(I).LE.NAM1) GO TO 50
          ELSE IF (KZ(5).EQ.2) THEN
              IF (NAME(I).LE.NAM1.OR.
     &           (NAME(I).GT.N1.AND.NAME(I).LE.NAM2)) GO TO 50
          END IF
          ZM2 = ZM2 + BODY(I)
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
   50 CONTINUE
*
*       Sort square distances with respect to the centre C.
      CALL SORT1(NP,R2,JLIST)

      DO 60 IL = 1,LX
          ZM = 0.0
          ZMH = FLAGR(IL)*ZM2
          I = 0
   55     I = I + 1
          IM = JLIST(I)
          ZM = ZM + BODY(IM)
          IF (ZM.LT.ZMH) GO TO 55
          RLAGR(IL) = SQRT(R2(I))
   60 CONTINUE
*
*       Obtain half-mass radius separately.
      ZM = 0.0
      ZMH = 0.5*ZM2
      I = 0
   65 I = I + 1
      IM = JLIST(I)
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZMH) GO TO 65
*
      RH2 = SQRT(R2(I))
      NP2 = NP
*
      WRITE (32,30)  TIME+TOFF, (LOG10(RLAGR(K)),K=1,LX)
      CALL FLUSH(32)
*
      WRITE (6,70)  NP1, NP2, ZM1, ZM2, RH1, RH2
   70 FORMAT(/,' TWO MASS GROUPS:    NM1 NM2 M1 M2 RM1 RM2 ',
     &                               2I7,1X,4F7.3)
*
      RETURN
*
      END
