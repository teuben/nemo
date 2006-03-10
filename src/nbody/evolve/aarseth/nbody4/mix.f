***
      SUBROUTINE MIX(J1,J2,DM)
*
*     Author : J. R. Hurley
*     Date :   7th July 1998
*
*       Evolution parameters for mixed star.
*       ------------------------------------
*
      INCLUDE 'common4.h'
*
      REAL*8 TSCLS(20),LUMS(10),GB(10),TMS1,TMS2,TMS3,TN
      REAL*8 M01,M02,M03,M1,M2,M3,AGE1,AGE2,AGE3,MC3,LUM,RM,RCC
      REAL*8 MENV,RENV,K2E
      REAL*8 MCH,MXNS
      PARAMETER(MCH=1.44D0,MXNS = 3.0d0)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*
*       Define global indices with body #I1 being most evolved.
      IF(KSTAR(J1).GE.KSTAR(J2))THEN
          I1 = J1
          I2 = J2
      ELSE
          I1 = J2
          I2 = J1
      END IF
*
*       Specify case index for collision treatment.
      K1 = KSTAR(I1)
      K2 = KSTAR(I2)
      ICASE = KTYPE(K1,K2)
      IF(ICASE.GT.100) WRITE(38,*)' MIX ERROR ICASE>100 ',ICASE,K1,K2
*
*       Record name and stellar type of the colliding stars.
      JC = NCOLL + 1
      ITYPE(JC,1) = NAME(I1)
      ITYPE(JC,2) = NAME(I2)
      ITYPE(JC,3) = KSTAR(I1)
      ITYPE(JC,4) = KSTAR(I2)
*
*       Set physical time and initialize mass loss & time increment.
      TTOT = TIME + TOFF
      TPHYS = TTOT*TSTAR
      DM = 0.d0
      DT = 0.d0
      RS1 = RADIUS(I1)*SU
      RS2 = RADIUS(I2)*SU
*
*       Evolve the stars to the current time unless they have been
*       evolved further by recent Roche interaction.
      TEV1 = MAX(TIME,TEV0(I1))
*
*       Determine evolution time scales for first star.
      M01 = BODY0(I1)*ZMBAR
      M1 = BODY(I1)*ZMBAR
      AGE1 = TEV1*TSTAR - EPOCH(I1)
      CALL star(K1,M01,M1,TMS1,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Obtain time scales for second star.
      M02 = BODY0(I2)*ZMBAR
      M2 = BODY(I2)*ZMBAR
      AGE2 = TEV1*TSTAR - EPOCH(I2)
      CALL star(K2,M02,M2,TMS2,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Check for planetary systems - defined as HeWDs and low-mass WDs!
      IF(K1.EQ.10.AND.M1.LT.0.05)THEN
         ICASE = K2
         IF(K2.LE.1)THEN
            ICASE = 1
            AGE1 = 0.D0
         ENDIF
      ELSEIF(K1.GE.11.AND.M1.LT.0.5.AND.ICASE.EQ.6)THEN
         ICASE = 9
      ENDIF
      IF(K2.EQ.10.AND.M2.LT.0.05)THEN
         ICASE = K1
         IF(K1.LE.1)THEN
            ICASE = 1
            AGE2 = 0.D0
         ENDIF
      ENDIF
*
      WRITE(38,67)K1,M01,M1
 67   FORMAT(' MIX OLD *1:',I4,2F7.2)
      WRITE(38,68)K2,M02,M2
 68   FORMAT(' MIX OLD *2:',I4,2F7.2)
*
*       Specify total mass.
      M3 = M1 + M2
      M03 = M01 + M02
      MC3 = 0.d0
      KW = ICASE
      AGE3 = 0.d0
      TMS3 = 0.d0
*
*       Evaluate apparent age and other parameters.
*
      IF(ICASE.EQ.1)THEN
*       Specify new age based on complete mixing.
         IF(K1.EQ.7) KW = 7
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = 0.1d0*TMS3*(AGE1*M1/TMS1 + AGE2*M2/TMS2)/M3
      ELSEIF(ICASE.EQ.3.OR.ICASE.EQ.6.OR.ICASE.EQ.9)THEN
         MC3 = M1
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
      ELSEIF(ICASE.EQ.4)THEN
         MC3 = M1
         AGE3 = AGE1/TMS1
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
      ELSEIF(ICASE.EQ.7)THEN
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = TMS3*(AGE2*M2/TMS2)/M3
      ELSEIF(ICASE.LE.12)THEN
*       Ensure that a new WD has the initial mass set correctly.
         M03 = M3
         IF(ICASE.LT.12.AND.M3.GE.MCH)THEN
            DM = M3
            M3 = 0.D0
            KW = 15
            RM = 0.D0
            WRITE (6,2)  NAME(I1), NAME(I2), K1, K2, KW, DM
    2       FORMAT (' NEW SN:    NAM K* M3 ',2I6,3I4,F6.2)
         ENDIF
         DT = 1.0d+04
      ELSEIF(ICASE.EQ.13.OR.ICASE.EQ.14)THEN
*       Set unstable Thorne-Zytkow object with fast mass loss of envelope
*       unless the less evolved star is a WD, NS or BH.
         IF(K2.LT.10)THEN
            M03 = M1
            M3 = M1
            DM = M2
            NTZ = NTZ + 1
            WRITE (6,5)  NAME(I1), NAME(I2)
    5       FORMAT (' NEW TZ    NM ',2I6)
         ELSEIF(ICASE.EQ.13)THEN
            NGB = NGB + 1
            IF(M3.GT.MXNS) KW = 14
         ENDIF
         DT = 1.0D+04
      ELSEIF(ICASE.EQ.15)THEN
         DM = M3
         M3 = 0.D0
      ELSEIF(ICASE.GT.100)THEN
*       Common envelope case which should only be used after COMENV.
         KW = K1
         AGE3 = AGE1
         M3 = M1
         M03 = M01
      ELSE
*       This should not be reached.
        WRITE(6,*)' ERROR MIX: ICASE NOT CAUGHT!!!'
        WRITE(6,*)' K1 K2 -> K3 ',K2,K2,KW
        CALL gpfree
        STOP
      ENDIF
*
      WRITE(38,69)KW,M03,M3
 69   FORMAT(' MIX NEW *3:',I4,2F7.2)
*
*       Determine consistent stellar type and specify mass loss.
      IF(KW.LE.14)THEN
         KW0 = KW 
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         CALL hrdiag(M03,AGE3,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW,MC3,RCC,MENV,RENV,K2E)
         if(kw.ne.kw0)then
            write(38,*)' MIX TYPE CHANGE ',kw0,kw,kstar(i1),icase
         endif
      ENDIF
      KSTAR(I1) = KW
      RADIUS(I1) = RM/SU
      ZLMSTY(I1) = LUM
*       Update initial mass, epoch & evolution times.
      BODY0(I1) = M03/ZMBAR
      EPOCH(I1) = TEV1*TSTAR - AGE3
      TEV(I1) = TEV1 + DT
      TEV0(I1) = TEV1
      DM = DM/ZMBAR
*
*       Record final type of mixed star.
      ITYPE(JC,5) = KSTAR(I1)
*
*       Check for blue straggler formation (TM < TPHYS & KSTAR <= 1).
      IF(TMS3.LT.TPHYS.AND.KSTAR(I1).LE.1)THEN
          WRITE (6,8)  NAME(I1), M3, TMS3, TPHYS, AGE3
    8     FORMAT (' NEW BS (MIX):    NAM M3 TM TP AGE ',
     &                               I6,F6.2,3F8.1)
      ENDIF
*
      WRITE (6,10)  IQCOLL, K1, K2, KSTAR(I1), M1, M2, M3, RS1, RS2,
     &              RADIUS(I1)*SU, DM*ZMBAR
   10 FORMAT (' NEW STAR:    IQ K1 K2 K* M1 M2 M3 R1 R2 R* DM ',
     &                       4I3,3F6.1,3F7.1,F5.1)
*
*       Open unit #13 the first time.
      IF(FIRST)THEN
          OPEN (UNIT=13,STATUS='NEW',FORM='FORMATTED',FILE='COLL')
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          IF(NCOLL.EQ.0)THEN
              WRITE (13,20)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &                       NBIN0, NZERO
   20         FORMAT (/,6X,'MODEL:    RBAR =',F5.1,'  <M> =',F6.2,
     &                     '  M1 =',F6.1,'  TSCALE =',F6.2,
     &                     '  NB =',I4,'  N0 =',I6,//)
          ENDIF
*
          WRITE (13,25)
   25     FORMAT ('   TIME  NAME  NAME  K1  K2  KC  M1   M2   MC',
     &            '   DM    R1     R2    r/Rc   R     P',/)
      ENDIF
*
*       Form central distance (scaled by RC) and period in days.
      RI2 = 0.d0
      RIJ2 = 0.d0
      VIJ2 = 0.d0
      DO 30 K = 1,3
          RI2 = RI2 + (X(K,I1) - RDENS(K))**2
          RIJ2 = RIJ2 + (X(K,I1) - X(K,I2))**2
          VIJ2 = VIJ2 + (XDOT(K,I1) - XDOT(K,I2))**2
   30 CONTINUE
      RI = SQRT(RI2)/RC
      RIJ = SQRT(RIJ2)
      ZMB = BODY(I1) + BODY(I2)
      SEMI = 2.d0/RIJ - VIJ2/ZMB
      SEMI = 1.d0/SEMI
      TK = DAYS*SEMI*SQRT(ABS(SEMI)/ZMB)
*
*       Accumulate collision diagnostics on unit #13.
      WRITE (13,35)  TTOT, (ITYPE(JC,K),K=1,5), M1, M2, M3, DM*ZMBAR,
     &               RS1, RS2, RI, RIJ*SU, TK
   35 FORMAT (1X,F6.1,2I6,3I4,4F5.1,2F9.2,F6.2,F9.2,1P,E9.1)
      CALL FLUSH(13)
*
*       Re-define indices of colliding bodies with J1 as new c.m.
      J1 = I1
      J2 = I2
*
      IF(KSTAR(I1).GT.12)THEN
          WRITE (15,40)  K1, K2, KSTAR(I1), BODY(I1)*ZMBAR,
     &                   BODY(I2)*ZMBAR, M3
   40     FORMAT (' MIX:    K1 K2 K* M1 M2 M3 ',3I4,3F7.2)
          CALL FLUSH(15)
      ENDIF
*
      RETURN
      END
***
