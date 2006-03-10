      SUBROUTINE REINIT(ISUB)
*
*
*       Re-initialization of chain system.
*       ----------------------------------
*
      INCLUDE 'common4.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,ANG(3),FIRR(3),FD(3)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX),LISTCM(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/ECHAIN/  ECH
*
*
*       Set phase indicator for step reduction (routine STEPS).
      IPHASE = 8
*
*       Obtain current coordinates & velocities of all particles (for FPOLY).
      CALL XVPRED(IFIRST,NTOT)
*
*       Initialize force polynomial and add differential force corrections.
      CALL FPOLY1(ICH,ICH,0)
      DO 5 K = 1,3
          FIRR(K) = 0.0
          FD(K) = 0.0
    5 CONTINUE
      CALL CHF(ICH,X(1,ICH),XDOT(1,ICH),FIRR,FD)
      DO 10 K = 1,3
          F(K,ICH) = F(K,ICH) + 0.5*FIRR(K)
          FDOT(K,ICH) = FDOT(K,ICH) + ONE6*FD(K)
   10 CONTINUE
*
*       See whether body #ICH should be added to NLIST.
*     IF (T0(ICH) + STEP(ICH).LT.TLIST) THEN
*         CALL NLMOD(ICH,1)
*     END IF
*
*       Obtain maximum unperturbed separation based on dominant neighbour.
      CALL EXTEND(ISUB)
*
*       Initialize perturber list for chain (copy from EXTEND).
      NNB = ILIST(1)
      DO 20 L = 2,NNB+1
          LISTCM(L) = ILIST(L)
   20 CONTINUE
      LISTCM(1) = NNB
      CALL CHLIST(ICH)
*
*       Update decision-making variables for chain regularization.
      TS(ISUB) = TIME
      STEPS(ISUB) = 0.01*STEP(ICH)
*
*       Re-calculate new energy of chain system (just in case).
      CALL CONST(XCH,VCH,M,NN,ECH,ANG,ALAG)
*
*       Set phase indicator < 0 to ensure new NLIST in routine INTGRT.
      IPHASE = -1
*
*       Specify non-zero indicator for sending new data to GRAPE in INTGRT.
      ISEND = -1
*
      RETURN
*
      END
