      SUBROUTINE UNITS
*
*
*       Initialization of units & scaling factors.
*       ------------------------------------------
*
      INCLUDE 'common4.h'
*
      
*       Define GM & PC in cgs units and AU in pc.
      GM = 6.67D-08*1.989D+33
      PC = 3.0857D+18
      AU = 2.0627E+05
*
*       Form scaling factors for binary periods A*SQRT(A/M) to yrs and days.
      YRS = (RBAR*AU)**1.5/SQRT(ZMBAR)
      DAYS = 365.24*YRS
*
*       Specify conversion factors for lengths to solar radii & AU.
      SU = PC/(AU*6.96D+10)*RBAR*AU
      RAU = RBAR*AU
*
*       Copy solar mass scaling to new variable (M = BODY*<M>).
      SMU = ZMBAR
*
*       Form time scale in seconds and velocity scale in km/sec.
      TSTAR = SQRT(PC/GM)*PC
      VSTAR = 1.0D-05*SQRT(GM/PC)
*
*       Convert time scale from units of seconds to million years.
      TSTAR = TSTAR/(3.15D+07*1.0D+06)
*
*       Ensure ZMBAR & RBAR > 0 (=0: assume <M>/Sun = 1, RBAR = 1 pc).
      IF (ZMBAR.LE.0.0D0) ZMBAR = FLOAT(N)/ZMASS
      IF (RBAR.LE.0.0D0) RBAR = 1.0
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      TSTAR = TSTAR*SQRT(RBAR**3/(ZMASS*ZMBAR))
      VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
*
*       Copy TSTAR to secondary time-scale factor.
*     TSTAR = 2.62d0*TSTAR
      TSCALE = TSTAR
*
*       Physical scaling: X, M, V, T from RBAR*X, ZMBAR*M, VSTAR*V, TSTAR*T.
      WRITE (6,10)  RBAR, ZMBAR, VSTAR, TSTAR, BODYM*ZMBAR, SU
   10 FORMAT (/,12X,'PHYSICAL SCALING:    R* =',F4.1,'  M* =',F8.1,
     &              '  V* =',F6.3,'  T* =',F6.3,'  <M> =',F5.2,
     &              '  SU =',1P,E8.1,/)
*
      RETURN
*
      END
