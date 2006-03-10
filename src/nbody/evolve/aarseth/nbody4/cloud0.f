      SUBROUTINE CLOUD0
*
*
*       Cloud initialization.
*       ---------------------
*
      INCLUDE 'common4.h'
      COMMON/CLOUDS/  XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &                CLM(MCL),CLMDOT(MCL),VCL,SIGMA,RB2,RB3,PCL2,
     &                TCL,STEPCL,TBIG,DTBIG,NCL,NEWCL
*
*
*       Initialize cloud variables.
      NCL = 0
      TCL = 0.0D0
      NEWCL = 0
      PCL2 = 0.0
      TBIG = 1.0E+10
*
*       Read the cloud parameters.
      READ (5,*)  NCL, RB2, VCL, SIGMA, DTBIG
      READ (5,*)  (CLM(J), RCL2(J),J=1,NCL)
      WRITE (6,5)  NCL, RB2, VCL, SIGMA
    5 FORMAT (/,12X,'NCL =',I3,'  RB =',F6.1,'  <V> =',F5.1,
     &                                       '  SIGMA =',F5.1)
      RBAR1 = RBAR
      IF (RBAR.EQ.0.0) RBAR1 = 1.0
*       Set cloud parameters in scaled units.
      RB2 = RB2/RBAR1
*       Form rms velocity of cluster members in km/sec.
      A1 = 0.047*SQRT(ZMASS*ZMBAR/RBAR1)
*       Define velocity unit and scaled cloud velocity & dispersion.
      A2 = A1/SQRT(0.5D0*ZMASS)
      VCL = VCL/A2
      SIGMA = SIGMA/A2
*
*       Specify conservative cloud integration step from crossing time.
      STEPCL = 0.002*TCR*RB2/VCL
*
*       Quantize the cloud integration step.
      DT = STEPCL
      CALL STEPK(DT,DTN)
      STEPCL = DTN
*
*       Scale radii & masses to model units.
      DO 10 J = 1,NCL
          RCL2(J) = RCL2(J)/RBAR1
          CLM(J) = CLM(J)/ZMBAR
   10 CONTINUE
*
*       Set time-scale for 'sun-rise' 0.05 of the cloud crossing time.
      CLDOT = 0.1*RB2/VCL
      CLDOT = 1.0/CLDOT
*
      WRITE (6,15)  RB2, VCL, SIGMA, STEPCL, (CLM(J),RCL2(J),J=1,NCL)
   15 FORMAT (/,12X,'CLOUDS   ',3F7.1,F9.6,9(F7.1,F5.1))
*
*       Define the square of cloud half-mass radii & growth times.
      DO 20 J = 1,NCL
          RCL2(J) = RCL2(J)**2
          CLMDOT(J) = CLM(J)*CLDOT
   20 CONTINUE
*
*       Include treatment of GMC (specify TBIG and reduce membership).
      IF (KZ(12).EQ.-2) THEN
          WRITE (6,22)  DTBIG, CLM(NCL), SQRT(RCL2(NCL))*RBAR
   22     FORMAT (/,12X,'GMC:    DT =',F7.1,'  MB =',1P,E8.1,
     &                                      '  RB =',0P,F6.1)
          TBIG = DTBIG
          NCL = NCL - 1
      END IF
*
*       Set square boundary radius, impact parameter & inverse cube.
      RB2 = RB2**2
      PCL2 = RB2
      RB3 = 1.0/(RB2*SQRT(RB2))
*
*       Define density centre for routine CLOUD.
      DO 25 K = 1,3
          RDENS(K) = 0.0
   25 CONTINUE
*
*       Initialize new clouds on the boundary.
      DO 30 ICL = 1,NCL
          CALL CLOUD(ICL)
   30 CONTINUE
*
      RETURN
*
      END
