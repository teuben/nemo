      SUBROUTINE GLOBAL
*
*
*       Global cluster diagnostics.
*       ---------------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      INTEGER MF(30)
      CHARACTER *110 LINE
      CHARACTER *10 TEMP
      INTEGER START
      LOGICAL IM
*
*
*       Determine global anisotropy ratio (eq. (19) of FH).
      VT2 = 0.0
      VR2 = 0.0
      DO 10 I = IFIRST,NTOT
          RI2 = 0.0
          VI2 = 0.0
          RD = 0.0
          DO 5 K = 1,3
              RI2 = RI2 + (X(K,I) - RDENS(K))**2
              VI2 = VI2 + XDOT(K,I)**2
              RD = RD + (X(K,I) - RDENS(K))*XDOT(K,I)
    5     CONTINUE
          RD = RD/SQRT(RI2)
          VT2 = VT2 + (VI2 - RD**2)
          VR2 = VR2 + RD**2
   10 CONTINUE
      ANIS = 1.0 - 0.5*VT2/VR2
*
*       Include mean stellar mass in solar units.
      ZMB = ZMASS*ZMBAR/FLOAT(N)
*
      WRITE (9,30)  TTOT, N, ZMASS, ZMB, ZMDOT, RSCALE, RTIDE, ANIS,
     &              POT-ETIDE, ZKIN, ETIDE
   30 FORMAT (' GLOBAL ',F7.1,I7,F7.2,F6.2,F7.1,2F6.2,F7.3,3F8.4)
*
*       Initialize logical variable and IMF histogram.
      IM = .TRUE.
      DO 40 I = 1,30
         MF(I) = 0
   40 CONTINUE
*
*       Form the IMF in logarithmic bins (count any ghosts separately).
      IMERGE = 0
      DO 50 I = 1,N
*       Note that KSTAR(I) & RADIUS(I) are available for every particle.
          IF (BODY(I).GT.0.0D0) THEN
              IBIN = 10.0*LOG10(ZMBAR*BODY(I)) + 11
              IF (I.LT.IFIRST) THEN
*       Set KS index and check corresponding c.m. NAME for merger.
                  JP = KVEC(I)
                  IF (NAME(N+JP).LT.0) THEN
*       Treat each KS component separately using sequential merger index.
                      IF (I.EQ.2*JP-1) THEN
                          IMERGE = IMERGE + 1
*       Copy the two original KS components directly from merger table.
                          DO 45 K = 1,2
                              BODYK = CM(K,IMERGE)
                              IBIN = 10.0*LOG10(ZMBAR*BODYK) + 11
                              MF(IBIN) = MF(IBIN) + 1
   45                     CONTINUE
                          GO TO 50
*       See whether outer component is single or KS (initialized in MERGE).
                      ELSE IF (CM(3,IMERGE).GT.0.0D0) THEN
*       Copy the two original outer KS components directly from merger table.
                          DO 48 K = 3,4
                              BODYK = CM(K,IMERGE)
                              IBIN = 10.0*LOG10(ZMBAR*BODYK) + 11
                              MF(IBIN) = MF(IBIN) + 1
   48                     CONTINUE
                          GO TO 50
                      END IF
                  END IF
              END IF
          ELSE
              IBIN = 30
          END IF
          IF (IBIN.GE.1.AND.IBIN.LE.30) THEN
              MF(IBIN) = MF(IBIN) + 1
          ELSE
              IF (IM) THEN
                 WRITE (6,*) ' ERROR IN MF; ZMBAR, BODY(I), IBIN:',
     &                         ZMBAR, BODY(I), IBIN
                 IM = .FALSE.
              END IF
          END IF
   50 CONTINUE
*
*       Prepare the IMF output for one line (Douglas Heggie 1/11/95).
      START = 1
      WRITE (LINE,'(110X)')
      DO I = 1,30
         IF (MF(I).GT.0) THEN
            LENGTH = ALOG10(FLOAT(MF(I))) + 1
         ELSE
            LENGTH = 1
         ENDIF
         WRITE (TEMP,'(I10)') MF(I)
         IF (START+LENGTH-1.GT.110) THEN
            WRITE (6,*) 'OVERFLOW IN IMF OUTPUT'
            GOTO 60
         ELSE
            LINE(START:START+LENGTH-1) = TEMP(10-LENGTH+1:)
            START = START+LENGTH
            LINE(START:START) = ' '
            START = START + 1
         ENDIF
      ENDDO
   60 CONTINUE
      WRITE (21,'(F9.1,1X,110A)') TTOT,LINE
*
      RETURN
*
      END
