      SUBROUTINE HOTSYS
*
*
*       Hot initial system.
*       -------------------
*
      INCLUDE 'common4.h'
*
*
*       Determine the rms velocity from current kinetic energy.
      ZK = 0.0D0
      DO 10 I = 1,N
          DO 5 K = 1,3
              ZK = ZK + BODY(I)*XDOT(K,I)**2
    5     CONTINUE
   10 CONTINUE
      VRMS = SQRT(ZK/ZMASS)
*
*       Read central velocity dispersion and form scaling factor.
      READ (5,*)  SIGMA0
      VSCALE = SIGMA0/(VSTAR*VRMS)
*
*       Scale the velocities to central velocity dispersion of SIGMA0 km/sec.
      DO 20 I = 1,N
          DO 15 K = 1,3
              XDOT(K,I) = VSCALE*XDOT(K,I)
   15     CONTINUE
   20 CONTINUE
*
*       Rescale crossing time, output times & termination time.
      RATIO = SIGMA0/VSTAR
      TCR = TCR/RATIO
      TCR0 = TCR
      DTADJ = DTADJ/RATIO
      DELTAT = DELTAT/RATIO
      TCRIT = TCRIT/RATIO
*
      WRITE (6,30)  SIGMA0, VRMS, VSCALE
   30 FORMAT (/,12X,'HOT SYSTEM:   SIGMA0 =',F5.1,'  <V> =',F6.3,
     &                                              '  VSCALE =',F6.3,/)
*
      RETURN
*
      END
