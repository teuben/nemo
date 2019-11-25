      SUBROUTINE INPUT
*
*
*       Parameter input.
*       ----------------
*
      INCLUDE 'common1.h'
      EXTERNAL VERIFY
*
*
*       Make a formal call to define input parameters & counters.
      CALL DEFINE
*
*       Read & print the main input parameters.
      READ (5,*)  N, NFIX, NRAND, NRUN
      READ (5,*)  ETA, DELTAT, TCRIT, QE, EPS
      READ (5,*)  (KZ(J),J=1,15)
*
      WRITE (6,10)
   10 FORMAT (////,12X,'N  NRAND    ETA   DELTAT   TCRIT   QE',
     &                                                    '        EPS')
      WRITE (6,20)  N, NRAND, ETA, DELTAT, TCRIT, QE, EPS
   20 FORMAT (/,8X,I5,I7,F8.3,F8.1,F8.1,1P,2E10.1)
      WRITE (6,30)  (J,J=1,15)
   30 FORMAT (//,12X,'OPTIONS   ',15I4)
      WRITE (6,40)  (KZ(J),J=1,15)
   40 FORMAT (/,22X,15I4)
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
*       Set square softening parameter.
      EPS2 = EPS**2
*       Copy random number skip for routine DATA.
      NLIST(1) = NRAND
      ETA0 = ETA
*
      RETURN
*
      END
