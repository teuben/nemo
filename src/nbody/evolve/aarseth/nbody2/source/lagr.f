      SUBROUTINE LAGR(C)
*
*
*       Lagrangian radii.
*       -----------------
*
      INCLUDE 'common2.h'
      COMMON/WORK1/  R2(NMAX)
      REAL*4  C(3),FLAGR(10),RLAGR(10)
*       Lagrangian radii at 1,2,5,7.5,10,15,20,50,75,85 % of total mass.
      DATA FLAGR  /1.E-2,2.E-2,5.E-2,7.5E-2,1.E-1,1.5E-1,2.E-1,
     &             5.E-1,7.5E-1,8.5E-1/
*
*
*       Set square radii of all particles.
      NP = 0
      DO 10 I = 1,N
          NP = NP + 1
          R2(NP) = (X(1,I) - C(1))**2 + (X(2,I) - C(2))**2 +
     &                                  (X(3,I) - C(3))**2
          JLIST(NP) = I
   10 CONTINUE
*
*       Sort square distances with respect to the centre C.
      CALL SORT2(NP,R2,JLIST)
*
*         Determine Lagrangian radii for specified mass fractions.
      DO 20 IL = 1,10
          ZM = 0.0
          ZMH = FLAGR(IL)*ZMASS
          I = 0
   15     I = I + 1
          IM = JLIST(I)
          ZM = ZM + BODY(IM)
          IF (ZM.LT.ZMH) GO TO 15
          RLAGR(IL) = SQRT(R2(I))
   20 CONTINUE
*
*       Determine half-mass radius separately.
      ZM = 0.0
      ZMH = 5.0E-01*ZMASS
      I = 0
   30 I = I + 1
      IM = JLIST(I)
      ZM = ZM + BODY(IM)
      IF (ZM.LT.ZMH) GO TO 30
*
*       Replace approximate half-mass radius by actual value.
      RSCALE = SQRT(R2(I))
*
*       Check output options (line printer or unit 7 or both).
      IF (KZ(7).EQ.2.OR.KZ(7).EQ.4) THEN
          WRITE (6,40)  (LOG10(RLAGR(K)),K=1,10)
   40     FORMAT (/,3X,'LAGR:  ',10F7.3)
      END IF
*
      IF (KZ(7).GE.3) THEN
          WRITE (7)  TIME, (LOG10(RLAGR(K)),K=1,10)
      END IF
*
      RETURN
*
      END
