      SUBROUTINE INPUT
*
*
*       Parameter input.
*       ----------------
*
      INCLUDE 'common2.h'
*
*
*       Make a formal call to define input parameters & counters.
      CALL DEFINE
*
*       Read & print the main input parameters.
      READ (5,*)  N, NFIX, NRAND, NNBMAX, NRUN
      READ (5,*)  ETAI, ETAR, RS0, DELTAT, TCRIT, QE, EPS
      READ (5,*)  (KZ(J),J=1,20)
*
      WRITE (6,10)
   10 FORMAT (/////,12X,'N  NFIX  NRAND  NNBMAX  NRUN')
      WRITE (6,12)  N, NFIX, NRAND, NNBMAX, NRUN
   12 FORMAT (/,I13,I6,I7,I8,I6)
      WRITE (6,15)
   15 FORMAT (//,12X,'ETAI      ETAR      RS0       DELTAT    TCRIT',
     &                                             '     QE        EPS')
      WRITE (6,20)  ETAI, ETAR, RS0, DELTAT, TCRIT, QE, EPS
   20 FORMAT (/,9X,1P,9E10.1)
      WRITE (6,24)  (J,J=1,20)
   24 FORMAT (//,12X,'OPTIONS   ',20I4)
      WRITE (6,26)  (KZ(J),J=1,20)
   26 FORMAT (/,22X,20I4)
*
*       Check optional input of Plummer model or logarithmic potential.
      IF (KZ(15).GT.0) THEN
          IF (KZ(15).EQ.1) THEN
              READ (5,*)  XTPAR(1), XTPAR(2)
              WRITE (6,30)  XTPAR(1), XTPAR(2)
   30         FORMAT (/,12X,'PLUMMER MODEL:   MASS =',F7.1,
     &                                                 '  SCALE =',F6.1)
              XTPAR(2) = XTPAR(2)**2
              CGAS = 0.0
          ELSE
              READ (5,*)  ZMGAS, RGAS
              CGAS = ZMGAS/RGAS
              WRITE (6,35)  ZMGAS, RGAS, CGAS
   35         FORMAT (/,12X,'LOGARITHMIC POTENTIAL:   MGAS =',F8.1,
     &                                  '  RGAS =',F6.1,'  CGAS =',F6.2)
              XTPAR(1) = 0.0
          END IF
      END IF
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
*       Set square softening and save initial time-step parameter.
      EPS2 = EPS**2
      ETA0 = ETAI
*       Copy random number skip for routine DATA.
      JLIST(1) = NRAND
*
*       Define neighbour membership range & initialize minimum radius.
      ZNBMIN = 0.2*FLOAT(NNBMAX)
      ZNBMAX = 0.9*FLOAT(NNBMAX)
      RSMIN = RS0
      RC = RS0
*       Temporary save of initial neighbour sphere for routine NBLIST.
*
      RETURN
*
      END
