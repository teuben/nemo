      SUBROUTINE IMF(BODY10,BODYN0)
*
*
*       Mass function for binaries & single stars.
*       ------------------------------------------
*
* KZ(20) = 
*
*     0,1 Salpeter power-law with exponent ALPHA.
*       2 Scalo (see Eggleton, Fitchett & Tout, Ap.J. 347, 998)
*       3 Eggleton book
*       4 KTG 93
*       5 KTG93 and binary masses from KTG91
*       6 Pavel99 extending to Brown-dwarf regime
*
*       ------------------------------------------
      INCLUDE 'common4.h'
      INTEGER NPLNT,ISPLIT
      REAL*4 RAN2
      REAL*8 MPLNT,BCM0(2*KMAX)
      REAL*8 IMFBD
      EXTERNAL IMFBD
      REAL*8 BCM,RHO
      COMMON/WORK1/ BCM(NMAX),RHO(NMAX)
      DATA G1,G2,G3,G4 /0.19,1.55,0.05,0.6/
      DATA G5,G6,G7,G8 /0.75,0.04,0.25,1.04/
*
*
*       Generate initial mass function.
      ZMASS = 0.D0
      KDUM = IDUM1
*
*       Determine whether binary mass is given by the mass of 1 or 2 stars 
*       (set by the sign of NBIN0). 
*       Note: if -10 < RATIO <= 0 in BINPOP then ISPLIT must equal 2.
      ISPLIT = 1
*       Note added 27/8/04: ignore Jarrod's NBIN0 < 0 and take ISPLIT = 2.
      ISPLIT = 2
      IF(NBIN0.LT.0)THEN
         ISPLIT = 2
         NBIN0 = -1*NBIN0
      ENDIF
*
*       Setup the possibility of including planets. 
      NPLNT = MIN(0,NBIN0)
      MPLNT = 0.001D0
      NBIN0 = NBIN0 - NPLNT
*
      IF(KZ(20).EQ.0)THEN
         ALPHA1 = ALPHA - 1.0
         FM1 = 1.0/BODY10**ALPHA1
         FMN = (FM1 - 1.0/BODYN0**ALPHA1)/(FLOAT(N) - 1.0)
         CONST = 1.0/ALPHA1
*
*       Assign individual masses sequentially.
         DO 10 I = 1,N
            FMI = FM1 - FLOAT(I - 1)*FMN
            BODY(I) = 1.0/FMI**CONST
            ZMASS = ZMASS + BODY(I)
 10      CONTINUE
         WRITE(6,12)ALPHA,BODY10,BODYN0
 12      FORMAT(/,12X,'STANDARD IMF    ALPHA =',F5.2,
     &                '  BODY1 =',F5.1,'  BODYN =',F5.2)
         GOTO 90
      ELSEIF(KZ(20).EQ.1)THEN
         ALPHA1 = ALPHA - 1.0
         FM1 = 1.0/BODY10**ALPHA1
         FMN = (FM1 - 1.0/BODYN0**ALPHA1)
         CONST = 1.0/ALPHA1
*
*       Assign individual masses.
         DO 15 I = 1,N
            XX = RAN2(KDUM)
            FMI = FM1 - XX*FMN
            BODY(I) = 1.D0/FMI**CONST
            IF(I.LE.NBIN0)THEN
               BCM0(2*I-1) = BODY(I)
               IF(ISPLIT.GT.1)THEN
                  XX = RAN2(KDUM)
                  FMI = FM1 - XX*FMN
                  BCM0(2*I) = 1.D0/FMI**CONST
                  BODY(I) = BODY(I) + BCM0(2*I)
               ELSE
                  BCM0(2*I) = 0.D0
               ENDIF
            ENDIF
            ZMASS = ZMASS + BODY(I)
 15      CONTINUE
      ELSEIF(KZ(20).EQ.2)THEN
         ITER = 1
         BODYM = BODY10
         XM = 0.998D0
*
*       Find initial value of XM for upper & lower mass by iteration.
 20      Y0 = (1.D0 - XM)**0.75 + 0.032D0*(1.D0 - XM)**0.25
         Y2 = (0.75D0 + 0.008D0/SQRT(1.D0 - XM))/
     &                          ((1.D0 - XM) + 0.032D0*SQRT(1.D0 - XM))
         Y1 = BODYM - 0.19D0*XM/Y0
         YPR = -0.19D0*(1.D0 + XM*Y2)/Y0
*
*       Try next guess of Newton-Raphson iteration.
         XM = XM - Y1/YPR
         IF(XM.GT.1.0) XM = 0.99999D0
         IF(XM.LT.0.0) XM = 0.001D0
         BODYS = 0.19D0*XM/Y0
         IF(ABS(BODYS - BODYM).GT.1.0D-06*BODYM) GOTO 20
*
*       Save upper value of XM and perform second iteration.
         IF(ITER.EQ.1)THEN
            X1 = XM
            BODYM = BODYN0
            XM = 0.4D0
            ITER = 2
            GOTO 20
         ENDIF
*
*       Assign individual masses.
         DX = (X1 - XM)
         DO 25 I = 1,N
            XX = RAN2(KDUM)
            XI = X1 - XX*DX
            ZM0 = (1.D0 - XI)**0.75 + 0.032D0*(1.D0 - XI)**0.25
            BODY(I) = 0.19D0*XI/ZM0
            IF(I.LE.NBIN0)THEN
               BCM0(2*I-1) = BODY(I)
               IF(ISPLIT.GT.1)THEN
                  XX = RAN2(KDUM)
                  XI = X1 - XX*DX
                  ZM0 = (1.D0 - XI)**0.75 + 0.032D0*(1.0 - XI)**0.25
                  BCM0(2*I) = 0.19D0*XI/ZM0
                  BODY(I) = BODY(I) + BCM0(2*I)
               ELSE
                  BCM0(2*I) = 0.D0
               ENDIF
            ENDIF
            ZMASS = ZMASS + BODY(I)
 25      CONTINUE
      ENDIF
*
      IF(KZ(20).LE.2) GOTO 40
*
*       Option for including extra massive stars.  
      IEXTR1 = 0
      IEXTR2 = -1
*
      DO 30 I = 1,N
 35      XX = RAN2(KDUM)
*
         IF(KZ(20).EQ.3)THEN
            ZM = 0.3D0*XX/(1.D0 - XX)**0.55
         ELSEIF(KZ(20).EQ.4.OR.KZ(20).EQ.5)THEN
            XX1 = 1.0 - XX
            IF(KZ(20).EQ.5.AND.I.LE.NBIN0)THEN
               ZM = 0.33D0*((1.D0/(XX1**G5 + G6*XX1**G7)) - (XX1**2/G8))
               ZM = 0.5D0*ZM
            ELSE
               ZM = 0.08D0 + (G1*XX**G2 + G3*XX**G4)/XX1**0.58
               IF(ZM.LT.0.15.AND.IEXTR1.LE.IEXTR2)THEN
                  IEXTR1 = IEXTR1 + 1
                  ZM = 20.D0 + FLOAT(IEXTR1-1)*20.D0/FLOAT(IEXTR2)
               ENDIF
            ENDIF
         ELSEIF(KZ(20).EQ.6)THEN
            ZM = IMFBD(XX,BODYN0,BODY10)
         ENDIF
*
*       See whether the mass falls within the specified range.
         IF(ZM.GE.BODYN0.AND.ZM.LE.BODY10)THEN
            IF(KZ(20).EQ.5.AND.I.LE.NBIN0) ZM = 2.D0*ZM
            BODY(I) = ZM
            IF(ISPLIT.GT.1.AND.I.LE.NBIN0.AND.KZ(20).NE.5)THEN
 37            XX = RAN2(KDUM)
               IF(KZ(20).EQ.3)THEN
                  ZM = 0.3D0*XX/(1.D0 - XX)**0.55
               ELSEIF(KZ(20).EQ.4)THEN
                  XX1 = 1.0 - XX
                  ZM = 0.08D0 + (G1*XX**G2 + G3*XX**G4)/XX1**0.58
               ELSEIF(KZ(20).EQ.6)THEN
                  ZM = IMFBD(XX,BODYN0,BODY10)
               ENDIF
               IF(ZM.LT.BODYN0.OR.ZM.GT.BODY10) GOTO 37
               BCM0(2*I-1) = MAX(BODY(I),ZM)
               BCM0(2*I) = MIN(BODY(I),ZM)
               BODY(I) = BODY(I) + ZM
            ELSEIF(I.LE.NBIN0)THEN
               BCM0(2*I-1) = ZM
               BCM0(2*I) = 0.D0
            ENDIF
            ZMASS = ZMASS + BODY(I)
         ELSE
            GOTO 35
         ENDIF
 30   CONTINUE
*
*       Place any planets around the first NPLNT single stars. 
 40   CONTINUE
      IF(NPLNT.GT.0)THEN
         DO 42 I = NBIN0+1,NBIN0+NPLNT
            BCM0(2*I-1) = BODY(I)
            BCM0(2*I) = MPLNT
            BODY(I) = BODY(I) + MPLNT
            ZMASS = ZMASS + MPLNT
 42      CONTINUE
         NBIN0 = NBIN0 + NPLNT
      ENDIF
*
*       See whether to skip mass sorting for binaries.
      IF(NBIN0.EQ.0) GOTO 50
*
*       Merge binary components in temporary variable for sorting.
      DO 44 I = 1,NBIN0
          BCM(I) = BODY(I)
          JLIST(I) = I
 44   CONTINUE
*
*       Sort total binary masses in increasing order.
      IF(NBIN0.GT.1)THEN
         CALL SORT1(NBIN0,BCM,JLIST)
      ENDIF
*
*       Copy the masses of binaries to COMMON in decreasing order.
      ZMB = 0.D0
      DO 46 I = 1,NBIN0
          BODY(I) = BCM(NBIN0-I+1)
          ZMB = ZMB + BODY(I)
          JB = JLIST(NBIN0-I+1)
          BODY0(2*I-1) = BCM0(2*JB-1)/ZMASS
          BODY0(2*I) = BCM0(2*JB)/ZMASS
 46   CONTINUE
*
      WRITE(6,48)NBIN0,BODY(1),BODY(NBIN0),ZMB/FLOAT(NBIN0)
 48   FORMAT(//,12X,'BINARY STAR IMF:    NB =',I6,
     &              '  RANGE =',1P,2E10.2,'  <MB> =',E9.2)
*
 50   IF(N.LE.NBIN0) GOTO 90
*
*       Move the single stars up to form compact array of N members.
      NS = 0
      DO 60 L = 1,N-NBIN0
         NS = NS + 1
         BCM(NS) = BODY(NBIN0+L)
         JLIST(NS) = NBIN0 + L
 60   CONTINUE
*
*       Sort masses of single stars in increasing order.
      CALL SORT1(NS,BCM,JLIST)
*
*       Copy the masses of single stars to COMMON in decreasing order.
      ZMS = 0.D0
      DO 70 I = 1,NS
         BODY(N-I+1) = BCM(I)
         ZMS = ZMS + BCM(I)
 70   CONTINUE
*
      WRITE(6,80)N-NBIN0,BODY(NBIN0+1),BODY(N),ZMS/FLOAT(N-NBIN0)
 80   FORMAT(/,12X,'SINGLE STAR IMF:    NS =',I6,'  RANGE =',1P,2E10.2,
     &           '  <MS> =',E9.2)
*
   90 RETURN
*
      END
