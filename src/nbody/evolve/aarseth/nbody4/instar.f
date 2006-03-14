      SUBROUTINE INSTAR
*
*
*       Initialization of stars.
*       ------------------------
*
      INCLUDE 'common4.h'
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN
      REAL*8 M0,M1,RM,LUM,AGE,MC,RCC,DLOGM
      REAL*8 MENV,RENV,K2,K3
      PARAMETER(K3=0.21D0)
      REAL*8 JSPIN,OSPIN
      REAL*8 VROTF
      EXTERNAL VROTF
*
*
*       Initialize mass loss variables & counters.
      TPHYS = 0.d0
      ZMRG = 0.d0
      ZMHE = 0.d0
      ZMRS = 0.d0
      ZMNH = 0.d0
      ZMWD = 0.d0
      ZMSN = 0.d0
      ZMBH = 0.d0
      ZMDOT = 0.d0
      ZMSY = 0.d0
      NMDOT = 0
      NROCHE = 0
      IQCOLL = 0
      NCHA = 0
      NSP = 0
      NHG = 0
      NRG = 0
      NHE = 0
      NRS = 0
      NNH = 0
      NWD = 0
      NSN = 0
      NTZ = 0
      NBS = 0
      NBR = 0
      NAS = 0
      NBH = 0
      NRO = 0
      NDD = 0
      NSPIR = 0
      NCIRC = 0
      NSLP = 0
      NCONT = 0
      NCOAL = 0
      NSTAB = 0
      NEINT = 0
      NEMOD = 0
      NHYP = 0
      NGB = 0
      NMS = N
      IBLUE = 0
*
      TMDOT = 1.0d+10
*
*     Set evolution parameters which depend solely on metallicity.
*
      CALL zcnsts(ZMET,ZPARS)
      WRITE(6,9)ZPARS(11),ZPARS(12),ZMET
 9    FORMAT(//,12X,'ABUNDANCES: X =',F6.3,' Y =',F6.3,' Z =',F6.3)
*
*       Calculate scale factor for spin angular momentum.
      SPNFAC = ZMBAR*SU**2/(1.0D+06*TSTAR)
*
      EPOCH1 = EPOCH0
      DO 10 I = 1,N
*
*       Obtain stellar parameters at current epoch (NB! dummy test).
          IF(KZ(22).LT.0)THEN
             READ(12,*)M1,KW,M0,EPOCH1,OSPIN
             IF(KW.LE.1.AND.(M1.LT.0.99*M0.OR.M1.GT.1.01*M0))THEN
                WRITE(6,*)' INSTAR ERROR M1 NE M0 FOR MS ',I,M1,M0
                STOP
             ENDIF
          ELSE
             M1 = BODY(I)*ZMBAR
             M0 = M1
             KW = 1
             IF(M0.LE.0.01D0) KW = 10
          ENDIF
          MC = 0.D0
          AGE = TIME*TSTAR - EPOCH1
          CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
          CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                RM,LUM,KW,MC,RCC,MENV,RENV,K2)
*
*       Assign initial spin angular momentum. 
          IF(KZ(22).LT.0)THEN
             JSPIN = (K2*(M1-MC)*RM**2 + K3*MC*RCC**2)*OSPIN
          ELSE
             OSPIN = 45.35d0*VROTF(M1)/RM
             JSPIN = K2*M1*RM**2*OSPIN
             IF(KW.GT.1) JSPIN = 0.D0
          ENDIF
*
*       Convert from solar radii to scaled units (assume Sun = 1/215 AU).
          IF (KZ(27).EQ.-2) THEN
              RADIUS(I) = 0.d0
              ZLMSTY(I) = 0.d0
              SPIN(I) = 0.d0
          ELSE
              RADIUS(I) = RM/SU
              ZLMSTY(I) = LUM
              SPIN(I) = JSPIN/SPNFAC
          END IF
*
*       Initialize the stellar classification type (KW = 0 - 10).
          KSTAR(I) = KW
*
*       Save the initial mass of all the stars in sequential order.
          BODY(I) = M1/ZMBAR
          BODY0(I) = M0/ZMBAR
*
*       Set initial look-up time.
          EPOCH(I) = TIME*TSTAR - AGE
          TEV0(I) = TIME
          CALL TRDOT(I,DTM,RM)
          TEV(I) = DTM
*
*       Determine the time for next stellar evolution check.
          IF (TEV(I).LT.TMDOT) THEN
              TMDOT = TEV(I)
          END IF
*
*         WRITE(13,*)M1,KW,M0,EPOCH(I)
   10 CONTINUE
*
*       Include option for Chernoff-Weinberg mass loss time-scale.
      IF (KZ(19).GT.4) THEN
          TMDOT = 1.0d+10
          DO 12 I = 1,N
              IF (KZ(13).GT.0) THEN
                  M0 = BODY(I)*ZMBAR
                  TM = 0.d0
                  CALL INTEV(M0,TM,DLOGM)
              ELSE
*       Set TEV for all masses in the case of variable time scale.
                  M0 = BODY(I)*ZMBAR
                  TM = 0.d0
                  CALL INTEV(M0,TM,DLOGM)
              END IF
              TEV(I) = TM/TSTAR
              IF (TEV(I).LT.TMDOT) THEN
                  TMDOT = TEV(I)
              END IF
   12     CONTINUE
          IF (KZ(13).NE.0) THEN
              M0 = BODY1*ZMBAR
              TM = 0.d0
              CALL INTEV(M0,TM,DLOGM)
              TEV10 = TM
          END IF
      END IF
*
*        Ensure binary components will be updated at the same time.
      DO 15 I = 1,NBIN0
         I1 = 2*I - 1
         I2 = I1 + 1
         TEV(I1) = MIN(TEV(I1),TEV(I2))
         TEV(I2) = TEV(I1)
   15 CONTINUE
*
*       Define first quantized step < 100 yrs (minimum interval for MDOT).
      DT = 1.0d-04/TSCALE
      CALL STEPK(DT,DTN)
      IF (DTN*TSCALE.LT.100.0) DTN = 2.0*DTN
      STEPX = DTN
*
*       Initialize stellar collision matrix.
*
      ktype(0,0) = 1
      do 20 , j = 1,6
         ktype(0,j) = j
         ktype(1,j) = j
 20   continue
      ktype(0,7) = 4
      ktype(1,7) = 4
      do 25 , j = 8,12
         if(j.ne.10)then
            ktype(0,j) = 6
         else
            ktype(0,j) = 3
         endif
         ktype(1,j) = ktype(0,j)
 25   continue
      ktype(2,2) = 3
      do 30 , i = 3,14
         ktype(i,i) = i
 30   continue
      ktype(5,5) = 4
      ktype(7,7) = 1
      ktype(10,10) = 15
      ktype(13,13) = 14
      do 35 , i = 2,5
         do 40 j = i+1,12
            ktype(i,j) = 4
 40      continue
 35   continue
      ktype(2,3) = 3
      ktype(2,6) = 5
      ktype(2,10) = 3
      ktype(2,11) = 5
      ktype(2,12) = 5
      ktype(3,6) = 5
      ktype(3,10) = 3
      ktype(3,11) = 5
      ktype(3,12) = 5
      ktype(6,7) = 4
      ktype(6,8) = 6
      ktype(6,9) = 6
      ktype(6,10) = 5
      ktype(6,11) = 6
      ktype(6,12) = 6
      ktype(7,8) = 8
      ktype(7,9) = 9
      ktype(7,10) = 7
      ktype(7,11) = 9
      ktype(7,12) = 9
      ktype(8,9) = 9
      ktype(8,10) = 7
      ktype(8,11) = 9
      ktype(8,12) = 9
      ktype(9,10) = 7
      ktype(9,11) = 9
      ktype(9,12) = 9
      ktype(10,11) = 9
      ktype(10,12) = 9
      ktype(11,12) = 12
      do 45 , i = 0,12
         ktype(i,13) = 13
         ktype(i,14) = 14
 45   continue
      ktype(13,14) = 14
*
* Increase common-envelope cases by 100.
      do 50 , i = 0,9
         do 55 , j = i,14
            if(i.le.1.or.i.eq.7)then
               if(j.ge.2.and.j.le.9.and.j.ne.7)then
                  ktype(i,j) = ktype(i,j) + 100
               endif
            else
               ktype(i,j) = ktype(i,j) + 100
            endif
 55      continue
 50   continue
*
*       Assign the remaining values by symmetry.
      do 60 , i = 1,14
         do 65 , j = 0,i-1
            ktype(i,j) = ktype(j,i)
 65      continue
 60   continue
*
      WRITE (6,75)  KTYPE
   75 FORMAT (/,11X,' KTYPE: ',15I4,14(/,19X,15I4))
*
      RETURN
      END
