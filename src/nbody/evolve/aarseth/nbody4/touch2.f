      SUBROUTINE TOUCH2(K1,K2,RCR)
*
*
*       Collision detector for chain pairs.
*       -----------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      PARAMETER  (NMX=10,NMX4=4*NMX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/EBSAVE/  EBS
*
*
      WRITE (6,1)  NAMEC(K1), NAMEC(K2), ISTAR(K1), ISTAR(K2),
     &             SIZE(K1), SIZE(K2), RCR, QPERI, EBS
    1 FORMAT (' CHAIN COLL    NAM K* R* RC R EB ',2I6,2I4,1P,5E10.2)
*
*       Set zero radii for binary components to avoid repeated events.
      IF (EBS.LT.0.0) THEN
          SIZE(K1) = 0.0D0
          SIZE(K2) = 0.0D0
      END IF
*
      RETURN
*
      END
