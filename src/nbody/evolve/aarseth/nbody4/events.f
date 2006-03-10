      SUBROUTINE EVENTS
*
*
*       Mass loss or tidal interaction events.
*       --------------------------------------
*
      INCLUDE 'common4.h'
      INTEGER  NTYPE(17)
*
      TTOT = TIME + TOFF
      TPHYS = TTOT*TSTAR
*
*       Determine turnoff mass at current cluster age (cf routine STAR).
      IF (TIME.LE.0.0D0) THEN
          TURN = BODY1*ZMBAR
          TURN0 = TURN
      ELSEIF(KZ(19).GE.3)THEN
          TURN0 = TURN
          T6 = TPHYS - EPOCH0
          CALL MTURN(TURN,T6,ZPARS)
      END IF
*
*       Check counter for stellar evolution events.
      IF (NMDOT.GT.0) THEN
          DO 5 J = 1,16
              NTYPE(J) = 0
    5     CONTINUE
*
          KM = 1
          RSTAR = 0.0
          DO 10 J = 1,N
              KW = KSTAR(J) + 1
              KW = MIN(KW,16)
              KW = MAX(KW,1)
              NTYPE(KW) = NTYPE(KW) + 1
              KM = MAX(KM,KW)
              RJ = RADIUS(J)
              RSTAR = MAX(RJ,RSTAR)
   10     CONTINUE
          NMS = NTYPE(1) + NTYPE(2)
*
          WRITE (6,15)
   15    FORMAT(/,5X,'NMDOT   NRG  NHE   NRS  NNH  NWD  NSN  NBH  TURN',
     &             '   ZMRG   ZMHE   ZMRS   ZMNH   ZMWD   ZMSN   ZMDOT',
     &             '  NTYPE')
         WRITE (6,20)  NMDOT, NRG, NHE, NRS, NNH, NWD, NSN, NBH, TURN,
     &                 ZMRG, ZMHE, ZMRS, ZMNH, ZMWD, ZMSN, ZMDOT,
     &                 (NTYPE(J),J=1,KM)
   20     FORMAT (' #4',I8,2I5,I6,4I5,F6.2,6F7.1,F8.1,I7,I6,14I4)
      END IF
*
*       Check output for tidal capture & collisions or Roche overflow.
      IF (NDISS + NCOLL + NRO.GT.0) THEN
*
*       Include check for BS (mass exceeding 1.1 of current turnoff mass).
          NBS2 = NBS
          DO 50 I = 1,IFIRST-1
              ZM = BODY(I)*ZMBAR
              IF(KSTAR(I).LE.1.AND.ZM.GT.1.02*TURN.AND.ZM.LT.TURN0)THEN
                  NBS = NBS + 1
                  AGE = TIME*TSTAR - EPOCH(I)
                  WRITE (6,55)  NAME(I), NBS, ZM, TURN, TPHYS, AGE
   55             FORMAT (' NEW BS:    NAM NBS M TURN TP AGE ',
     &                                 I6,I4,2F6.2,2F8.1)
              END IF
   50     CONTINUE
          WRITE (6,60)
   60     FORMAT (/,5X,'NDISS   NSPIR  NTD  NCOLL  NCOAL  NCONT  NBS',
     &               '  NHYP  NCH  NSP  NCIRC  NSLP  NROCHE  NRO  NAS',
     &               '  NBR  NDD  ECOLL  EMDOT  ECDOT  EKICK')
          WRITE (6,65)  NDISS, NSPIR, NTIDE, NCOLL, NCOAL, NCONT, NBS,
     &                  NHYP, NCHA, NSP, NCIRC, NSLP, NROCHE, NRO, NAS,
     &                  NBR, NDD, ECOLL, EMDOT, ECDOT, EKICK
   65     FORMAT (' #5',I7,I9,I4,3I7,I5,I6,2I5,I7,I6,I8,4I5,4F7.2)
          IF(NBS.NE.NBS2.AND.KZ(25).GE.2)THEN
             WRITE(98,70)TPHYS,N,NBS,TURN
 70          FORMAT(F9.1,I8,I4,2F6.2,' EVENTS')
             CALL FLUSH(99)
          ENDIF
      END IF
*
*       Check Roche look-up times in case of orbital shrinkage.
      IF(KZ(19).GT.5.AND.KZ(34).GT.0)THEN
      DO 80 IPAIR = 1,NPAIRS
          I = N + IPAIR
          IF (BODY(I).GT.0.0D0.AND.NAME(I).GT.0.AND.
     &        KSTAR(I).GT.0) THEN
              CALL TRFLOW(IPAIR,DTR)
              TEV(I) = MIN(TIME+DTR,TEV(I))
          END IF
   80 CONTINUE
      ENDIF
*
      RETURN
*
      END
