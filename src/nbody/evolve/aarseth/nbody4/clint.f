      SUBROUTINE CLINT
*
*
*       Cloud integration.
*       ------------------
*
      INCLUDE 'common4.h'
      COMMON/CLOUDS/  XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &                CLM(MCL),CLMDOT(MCL),VCL,SIGMA,RB2,RB3,PCL2,
     &                TCL,STEPCL,TBIG,DTBIG,NCL,NEWCL
      REAL*8  FC(3)
*
*
*       Set time-step and update the cloud reference time.
      DT = TIME - TCL
      DT1 = 0.5D0*DT
      TCL = TIME
      IBIG = 0
*
*       Integrate cloud orbits in rotating coordinates with tidal effects.
      DO 20 ICL = 1,NCL
          FC(1) = XCL(1,ICL)*TIDAL(1) + XDOTCL(2,ICL)*TIDAL(4)
          FC(2) = -XDOTCL(1,ICL)*TIDAL(4)
          FC(3) = XCL(3,ICL)*TIDAL(3)
          XCL2 = 0.0
*
          DO 10 K = 1,3
              XCL(K,ICL) = (FC(K)*DT1 + XDOTCL(K,ICL))*DT + XCL(K,ICL)
              XDOTCL(K,ICL) = FC(K)*DT + XDOTCL(K,ICL)
              XCL2 = XCL2 + (XCL(K,ICL) - RDENS(K))**2
   10     CONTINUE
*
*       Determine the minimum impact parameter (squared).
          IF (KZ(12).NE.-2) THEN
              PCL2 = MIN(XCL2,PCL2)
          ELSE IF (CLM(ICL).GT.CLM(1)) THEN
              PCL2 = MIN(XCL2,PCL2)
          END IF
*
*       Check for modification of the cloud mass (increase or reduction).
          IF (XCL2.LT.RB2) THEN
              IF (BODYCL(ICL).LT.CLM(ICL)) THEN
                  BODYCL(ICL) = BODYCL(ICL) + CLMDOT(ICL)*DT
              END IF
          ELSE
*       Reduce cloud mass gradually by 'sun-set' procedure.
              BODYCL(ICL) = BODYCL(ICL) - CLMDOT(ICL)*DT
          END IF
*
*       Initialize a new cloud when current mass becomes negative.
          IF (BODYCL(ICL).LE.0.0) THEN
              IF (CLM(ICL).LE.CLM(1)) THEN
                  CALL CLOUD(ICL)
              ELSE
                  IBIG = 1
              END IF
          END IF
   20 CONTINUE
*
      IF (IBIG.GT.0) NCL = NCL - 1
      IF ((TIME+TOFF)*TSTAR.GT.TBIG) THEN
          NCL = NCL + 1
          TBIG = TBIG + DTBIG
          CALL CLOUD(NCL)
      WRITE (6,25)  TTOT, TPHYS, SQRT(PCL2)*RBAR, CLM(NCL)*SMU
   25 FORMAT (' ADD BIG:   T TP PCM BODY ',2F9.2,F7.2,1P,E9.1)
      PCL2 = RB2
      CALL FLUSH(3)
      END IF
*
      RETURN
*
      END
