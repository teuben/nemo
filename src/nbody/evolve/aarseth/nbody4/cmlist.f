      SUBROUTINE CMLIST
*
*
*       Neighbour list for chain c.m.
*       -----------------------------
*
      INCLUDE 'common4.h'
      REAL*8  M,MASS,MC,MIJ,MKK
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX),LISTCM(LMAX)
      REAL*8  W1,W2,W3,RSEP2,RS2
*
*
*       Define characteristic neighbour distance.
      RS = RSCALE/FLOAT(N)**0.3333
*
*       Adopt generous neighbour distance from maximum of 2*RS & 100*RSUM.
      RS2 = MAX(4.0*RS**2,CMSEP2*RSUM**2)
*
*       Make neighbour list for chain c.m (skip ghosts).
    1 NNB1 = 1
      DO 10 J = IFIRST,NTOT
          IF (J.EQ.ICH.OR.BODY(J).LE.0.0D0) GO TO 10
          W1 = X(1,J) - X(1,ICH)
          W2 = X(2,J) - X(2,ICH)
          W3 = X(3,J) - X(3,ICH)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
          IF (RSEP2.LT.RS2) THEN
              NNB1 = NNB1 + 1
              LISTCM(NNB1) = J
              IF (NNB1.EQ.LMAX) THEN
                  RS2 = 0.5*RS2
                  GO TO 1
              END IF
          END IF
   10 CONTINUE
*
*       Increase distance on zero membership and try again.
      IF (NNB1.EQ.1) THEN
          RS2 = 2.0*RS2
          GO TO 1
      END IF
*
*       Save neighbour number.
      LISTCM(1) = NNB1 - 1
*
      RETURN
*
      END
