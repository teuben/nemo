      SUBROUTINE DIFSY1(N,EPS,H,X,Y)
*
*
*       Bulirsch-Stoer integrator.
*       --------------------------
*
*       Works if Gamma = (H - E)/L. For other time transformations EPS
*       must be scaled appropriately such that the test is esentially
*       of the form (H - E)/L < EPS.
*       Convergence test: ABS(Q'*DP) < EPS*TFAC & ABS(P'*DQ) < EPS*TFAC.
*       Reference: Mikkola (1987). In 'The Few Body Problem' p. 311.
*       This works only if eqs. are canonical and we have P's & Q's.
*       One additional eq. is allowed (e.g. for t'=...??) but not checked.
*
*
      PARAMETER  (NMX=80)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  Y(N),YA(NMX),YL(NMX),YM(NMX),DY(NMX),DZ(NMX),DT(NMX,7),
     &        D(7),X,XN,H,G,B,B1,U,V,C,TA,W,DABS
      COMMON/INTFAC/  LX,LE,LP,LV,LT,J10,NHALF2
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
*
      LOGICAL  KONV,BO,KL,GR,FYBAD,stopB
      COMMON/SLOW2/  stepl,stopB
*     REAL*8  EP(4)
*     DATA  EP/0.4D-1,0.16D-2,0.64D-4,0.256D-5/
      SAVE
*
*
      JTI=0
      FY=1.
      DO I=1,N
      YA(I)=Y(I)
      END DO
      CALL DERQP(Y(1),Y(LX),Y(LE),Y(LP),Y(LV),Y(LT),
     &           DZ(1),DZ(LX),DZ(LE),DZ(LP),DZ(LV),DZ(LT))
      IF(JC.GT.0)H=DSC
   10 XN=X+H
      BO=.FALSE.
      M=1
      JR=2
      JS=3
*     DO J=1,10
      DO J=1,J10
      IF(BO)THEN
      D(2)=16.D0/9.D0
      D(4)=64.D0/9.D0
      D(6)=256.D0/9.D0
      ELSE
      D(2)=2.25D0
      D(4)=9.D0
      D(6)=3.6D1
      END IF
      IF(J.GT.7)THEN
      L=7
      D(7)=6.4D1
      ELSE
      L=J
      D(L)=M*M
      END IF
      KONV=L.GT.3
      M=M+M
      G=H/(M)
      B=G+G
      M=M-1
      DO I=1,N
      YL(I)=YA(I)
      YM(I)=YA(I)+G*DZ(I)
      END DO
      DO K=1,M
      CALL DERQP(YM(1),YM(LX),YM(LE),YM(LP),YM(LV),YM(LT),
     &           DY(1),DY(LX),DY(LE),DY(LP),DY(LV),DY(LT))
*
*       Include procedure for terminating slow-down treatment.
      if (l.gt.3.and.stopB) then
      stepl=(k-1)*B+G
      return
      end if
*
      DO I=1,N
      U=YL(I)+B*DY(I)
      YL(I)=YM(I)
      YM(I)=U
      END DO
      END DO
      CALL DERQP(YM(1),YM(LX),YM(LE),YM(LP),YM(LV),YM(LT),
     &           DY(1),DY(LX),DY(LE),DY(LP),DY(LV),DY(LT))
      KL=L.LT.2
      GR=L.GT.5
      FS=0.
      DO I=1,N
      V=DT(I,1)
      C=(YM(I)+YL(I)+G*DY(I))*0.5D0
      DT(I,1)=C
      TA=C
      IF(.NOT.KL)THEN
      DO K=2,L
      B1=D(K)*V
      B=B1-C
      W=C-V
      U=V
      IF(B.NE.0.D0)THEN
      B=W/B
      U=C*B
      C=B1*B
      END IF
      V=DT(I,K)
      DT(I,K)=U
      TA=U+TA
      END DO
      IS=I+N/2
      IF(IS.GT.NHALF2)IS=I-(N/2)
      DYIS=DABS(DY(IS))
*     DYIS=DABS(DY(IS))*TFACIN
      IF(I.GT.NHALF2)DYIS=0.0D0
      IF(KONV)THEN
      TEST=DABS((Y(I)-TA)*DYIS)
      IF(TEST.GT.EPS) KONV=.FALSE.
      END IF
      IF(.NOT.GR)THEN
      FV=DABS(W)*DYIS
      IF(FS.LT.FV) FS=FV
      END IF
      END IF
      Y(I)=TA
      END DO
      IF(FS.NE.0.D0)THEN
      FA=FY
      K=L-1
*       Note factor 0.5 suggested by Seppo (22/3/99).
      FY=0.5D0*(EP(K)/FS)**(1./(L+K))
      FA7=0.7*FA
      IF(L.EQ.2)FA7=0.0D0
      FYBAD=.NOT.((FA7.GT.FY).OR.(FY.GT.0.7))
      IF(FYBAD)THEN
      H=H*FY
      JTI=JTI+1
      IF(JTI.GT.5)THEN
      H=0.0D0
      DO I=1,N
      Y(I)=YA(I)
      END DO
      RETURN
      END IF
      GOTO 10
      END IF
      END IF
      IF(KONV)THEN
      X=XN
      H=H*FY
      RETURN
      END IF
      D(3)=4.D0
      D(5)=1.6D1
      BO=.NOT.BO
      M=JR
      JR=JS
      JS=M+M
      END DO
      H=0.5*H
      GOTO 10
      END
