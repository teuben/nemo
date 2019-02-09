      SUBROUTINE XTRNLF(I,FIRR,FD)
*
*
*       External force & first derivative.
*       ----------------------------------
*
      INCLUDE 'common4.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
      REAL*8  XI(3),XIDOT(3),FIRR(3),FD(3),
     &        XG(3),XGDOT(3),FM(3),FMD(3),FS(3),FSD(3)
*
*
*       See whether to include a linearized galactic tidal force (two cases).
      IF (KZ(14).LE.2) THEN
          FIRR(1) = FIRR(1) + TIDAL(1)*X(1,I) + TIDAL(4)*XDOT(2,I)
          FIRR(2) = FIRR(2) - TIDAL(4)*XDOT(1,I)
          FIRR(3) = FIRR(3) + TIDAL(3)*X(3,I)
          FD(1) = FD(1) + TIDAL(1)*XDOT(1,I) + TIDAL(4)*FIRR(2)
          FD(2) = FD(2) - TIDAL(4)*FIRR(1)
          FD(3) = FD(3) + TIDAL(3)*XDOT(3,I)
      END IF
*
*       Consider point-mass, disk and/or logarithmic halo model.
      IF (KZ(14).EQ.3) THEN
          DO 2 K = 1,3
              XI(K) = X(K,I)
              XIDOT(K) = XDOT(K,I)
              XG(K) = RG(K) + XI(K)
              XGDOT(K) = VG(K) + XIDOT(K)
    2     CONTINUE
*       Employ differential instead of linearized forms for better accuracy.
          IF (GMG.GT.0.0D0) THEN
              CALL FNUC(RG,VG,FS,FSD)
              CALL FNUC(XG,XGDOT,FM,FMD)
              DO 10 K = 1,3
                  FIRR(K) = FIRR(K) + (FM(K) - FS(K))
                  FD(K) = FD(K) + (FMD(K) - FSD(K))
   10         CONTINUE
          END IF
*
*       Include Miyamoto disk for positive disk mass.
          IF (DISK.GT.0.0D0) THEN
              CALL FDISK(RG,VG,FS,FSD)
              CALL FDISK(XG,XGDOT,FM,FMD)
              DO 20 K = 1,3
                  FIRR(K) = FIRR(K) + (FM(K) - FS(K))
                  FD(K) = FD(K) + (FMD(K) - FSD(K))
   20         CONTINUE
          END IF
*
*       Check addition of logarithmic halo potential to regular force.
          IF (V02.GT.0.0D0) THEN
              CALL FHALO(RG,VG,FS,FSD)
              CALL FHALO(XG,XGDOT,FM,FMD)
              DO 30 K = 1,3
                  FIRR(K) = FIRR(K) + (FM(K) - FS(K))
                  FD(K) = FD(K) + (FMD(K) - FSD(K))
   30         CONTINUE
          END IF
      END IF
*
      RETURN
*
      END
