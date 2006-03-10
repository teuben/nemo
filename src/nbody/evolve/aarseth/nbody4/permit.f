      SUBROUTINE PERMIT(PERIM,IGO)
*
*
*       Check on existing multiple regularization.
*       ------------------------------------------
*
      INCLUDE 'common4.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
*
*
*       Search any existing subsystem.
      ISUB = 0
      DO 10 L = 1,NSUB
*       Distinguish between triple & quad case.
          IF (JCOMP.LE.N.AND.NAMES(4,L).EQ.0) THEN
              ISUB = L
      ELSE IF (JCOMP.GT.N.AND.NAMES(4,L).GT.0) THEN
              ISUB = L
          END IF
   10 CONTINUE
*
*       Do not allow a second regularization of the same type.
      IF (ISUB.GT.0) THEN
          IGO = 1
*       Enforce termination at next extension if new system < RSUM/2.
          IF (PERIM.LT.0.5*RSUM) THEN
              STEPS(ISUB) = 0.0
          END IF
      END IF
*
      RETURN
*
      END
