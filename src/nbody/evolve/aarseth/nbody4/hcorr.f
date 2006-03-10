      SUBROUTINE HCORR(I,DM,RNEW)
*
*
*       Mass loss correction of KS orbit.
*       ---------------------------------
*
      INCLUDE 'common4.h'
*
*
*       Copy pair index and determine secondary component (note I1 = I here).
      IPAIR = KSPAIR
      I1 = I
      I2 = 2*IPAIR - 1
      IF (I2.EQ.I1) THEN
          I2 = I2 + 1
      END IF
*
*       Define c.m. index and save old semi-major axis & reduced mass.
      J = N + IPAIR
      SEMI0 = -0.5D0*BODY(J)/H(IPAIR)
      ZMU0 = BODY(I1)*BODY(I2)/BODY(J)
*
*       Obtain energy change from M*A = const and H = -M/(2*A) (by DCH 8/96).
      DH = DM/SEMI0*(1.0 - 0.5*DM/BODY(J))
*
*       Reduce mass of c.m. body and subtract energy correction terms.
      BODY(J) = BODY(J) - DM
      ZMU1 = (BODY(I1) - DM)*BODY(I2)/BODY(J)
      EMD0 = EMDOT
      EMDOT = EMDOT - ZMU1*DH - (ZMU1 - ZMU0)*H(IPAIR)
*
*       Include diagnostic warning for large relative energy change.
      IF (DH.GT.0.2*ABS(H(IPAIR))) THEN
          WRITE (6,10) IPAIR, NAME(I), KSTAR(I), KSTAR(J), DH, H(IPAIR),
     &                 R(IPAIR)/SEMI0, DM*ZMBAR, BODY(I)*ZMBAR, RNEW*SU
   10     FORMAT (' WARNING!    LARGE CORRECTION',
     &                     '    KS NM K* DH H R/A DMS M R* ',
     &                          I4,I6,2I4,2F8.2,3F7.2,F8.1)
      END IF
*
*       Update binding energy due to mass loss DM and set new radius.
      H(IPAIR) = H(IPAIR) + DH
      IF (RADIUS(I).GT.0.0) THEN
          RADIUS(I) = RNEW
      END IF
*
*       Modify KS variables at constant eccentricity.
      CALL EXPAND(IPAIR,SEMI0)
*
      RETURN
*
      END
