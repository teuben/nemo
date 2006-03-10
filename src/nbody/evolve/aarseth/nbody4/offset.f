      SUBROUTINE OFFSET(DTOFF)
*
*
*       Offset of global times.
*       -----------------------
*
      INCLUDE 'common4.h'
      COMMON/BINARY/  BM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
*
*
*       Update the global offset time.
    1 TOFF = TOFF + DTOFF
*
*       Reduce all individual times and epochs by offset interval.
      DO 10 I = 1,NTOT
          T0(I) = T0(I) - DTOFF
          TEV(I) = TEV(I) - DTOFF
          TEV0(I) = TEV0(I) - DTOFF
          EPOCH(I) = EPOCH(I) - DTOFF*TSTAR
   10 CONTINUE
*
*       Update TNEXT for new scheduling at start of INTGRT.
      DO 15 I = IFIRST,NTOT
          TNEXT(I) = T0(I) + STEP(I)
   15 CONTINUE
*
*       Set new global times.
      TIME = TIME - DTOFF
      TADJ = TADJ - DTOFF
      TPRINT = TPRINT - DTOFF
      TPREV = TPREV - DTOFF
      TBLIST = TBLIST - DTOFF
      IF (KZ(19).GT.2) THEN
          TPLOT = TPLOT - DTOFF
          TMDOT = TMDOT - DTOFF
      END IF
      DO 30 I = 1,NMERGE
          TMDIS(I) = TMDIS(I) - DTOFF
   30 CONTINUE
      DO 40 I = 1,NSUB
          TS(I) = TS(I) - DTOFF
          T0S(I) = T0S(I) - DTOFF
   40 CONTINUE
*
*       Check tidal tail members (note TNEXT also needs doing).
      IF (NTAIL.GT.0) THEN
          DO 50 I = ITAIL0,NTTOT
              T0(I) = T0(I) - DTOFF
              TNEXT(I) = T0(I) + STEP(I)
   50     CONTINUE
      END IF
*
*       See whether more reductions are needed.
      IF (TIME.GE.DTOFF) GO TO 1
*
*       Activate control indicator for new scheduling.
      IPHASE = -1
*
      RETURN
*
      END
