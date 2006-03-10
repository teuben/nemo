      SUBROUTINE INTIDE
*
*
*       Input & scaling for tidal dissipation.
*       --------------------------------------
*
      INCLUDE 'common4.h'
*
*
*       Read parameters for tidal capture simulation.
      READ (5,*)  RSTAR, IMS, IEV, RMS, REV
*
*       Convert radii from S.U. to internal units.
      RSTAR = RSTAR/SU
*
      WRITE (6,10)  RSTAR, IMS, IEV, RMS, REV
   10 FORMAT (/,12X,'OLD TIDAL:    RSTAR =',1PE8.1,'  IMS =',0P,I5,
     &              '  IEV =',I4,'  RMS/R* =',F6.2,'  REV/R* =',F6.1)
*
*       Assign individual radii for main-sequence and evolved stars.
      DO 30 I = 1,N
*       Adopt a primitive scheme in case of no stellar evolution.
          IF (I.LE.IMS) THEN
              RADIUS(I) = RMS*RSTAR
          ELSE
              RADIUS(I) = REV*RSTAR
          END IF
*       Set a rough estimate of the luminosity.
          ZLMSTY(I) = (BODY(I)*ZMBAR)**4
   30 CONTINUE
*
*     WRITE (6,40)  BODY(1), BODY(N), RADIUS(1), RADIUS(N)
*  40 FORMAT (/,12X,'SCALED RADII:    M(1) =',F8.4,'  M(N) =',F8.4,
*    &                                '  R(1) =',1PE8.1,'  R(N) =',E8.1)
*
      RETURN
*
      END
