      PROGRAM EPIBIN
********************************************************************************
C     This program reads the file EPIRV with positions and velocities created 
C     by EPIC5 and bins the line-of-sight velocities for given position 
C     angles. EPIRV is supposed to contain less than 2000 points and its name 
C     should not contain more than 10 characters.
********************************************************************************
      IMPLICIT NONE

      INTEGER I,J,NPOINT,BETAI,NBOX,K
      REAL PI,XT(2000),YT(2000),VXT(2000),VYT(2000),BETA,
     &     RADVT(2000),VMIN,VMAX,VMEAN,DV,BOX(2000),VBOX(2000)
      REAL BOXMAX,XMIN,XMAX
      CHARACTER EPIRV*10,Y*1

      PI=3.14159265

      WRITE(*,*) '********************************************'
      WRITE(*,*) '** This program reads the file EPIRV and  **'
      WRITE(*,*) '** bins the line-of-sight velocities for  **'
      WRITE(*,*) '** choosen position angles.               **'
      WRITE(*,*) '** Version 2021-04-20                     **'
      WRITE(*,*) '********************************************'

C Read input file with positions and velocities

  200 WRITE(*,*) 'Input file: '
      READ(*,201) EPIRV
  201 FORMAT(A)
      OPEN(23,FILE=EPIRV,STATUS='OLD',ERR=994)
      I=0
 202  READ(23,101,END=10) XT(I+1),YT(I+1),VXT(I+1),VYT(I+1)
      WRITE(*,*) I+1,XT(I+1),YT(I+1),VXT(I+1),VYT(I+1)      
  101 FORMAT(1x,2F8.3,2F12.1)
      I=I+1
      GOTO 202
  994 WRITE(*,*)'CANNOT OPEN FILE', EPIRV
      GOTO 999
 10   NPOINT=I
      WRITE(*,*) 'found ',NPOINT,' entries in ',EPIRV
      IF (NPOINT.EQ.0) GOTO 999

C Compute distribution of line-of-sight velocities 

  213 WRITE(*,*)'GIVE POSITION ANGLE OF LINE OF NODES (BAR=90,TO QUIT
     &WRITE 360 OR MORE (integer)):'
      READ(*,*)BETAI
      IF (BETAI.GE.360) GOTO 999
      WRITE(*,*)'LINE-OF-SIGHT VELOCITIES PA=',BETAI
      BETA=BETAI*PI/180
      DO 3 I=1,NPOINT
         RADVT(I)= VXT(I)*COS(BETA)+VYT(I)*SIN(BETA)
    3 CONTINUE

C Find max and min values

      VMAX=RADVT(1)
      VMIN=RADVT(1)
      VMEAN=RADVT(1)
      DO 204 I=2,NPOINT
        IF (RADVT(I).GT.VMAX) VMAX=RADVT(I)
        IF (RADVT(I).LT.VMIN) VMIN=RADVT(I)
        VMEAN=VMEAN+RADVT(I)
  204 CONTINUE
      VMEAN=VMEAN/NPOINT
      DV=(VMAX-VMIN)/100

C Sort velocities

  212 WRITE(*,*)'VMIN= ',VMIN
      WRITE(*,*)'VMEAN= ',VMEAN
      WRITE(*,*)'VMAX= ',VMAX
      WRITE(*,*)'DV= ',DV
  211 FORMAT(F4.1)
      WRITE(*,*)'Write "y" or choose new vmin,vmax,dv :'
      READ(*,*,ERR=210) VMIN,VMAX,DV
      NBOX=INT((VMAX-VMIN)/DV)
      BOXMAX=0
  210 DO 208 J=1,NBOX
        BOX(J)=0.
        VBOX(J)=VMIN+J*DV
  208 CONTINUE
      DO 207 I=1,NPOINT
        J=INT((RADVT(I)-VMIN)/DV)
        BOX(J)=BOX(J)+1
        IF (BOX(J).GT.BOXMAX) BOXMAX=BOX(J)
 207  CONTINUE
      WRITE(*,*) 'BOXMAX=',BOXMAX

C Plot histogram, for PA=BETAI 

      DO 53 K=1,2
	IF (K.EQ.1) THEN
            CALL PGBEGIN(0,'/XWINDOW',1,1)
            ELSE
           CALL PGBEGIN(0,'histogram.ps/CPS',1,1)
        ENDIF
      CALL PGSCF(2)
      CALL PGENV(XMIN,XMAX,0.,BOXMAX,0,0)
      CALL PGLABEL('km/s','N','Velocity Distribution PA=')
      CALL PGIDEN
      CALL PGBIN(NBOX,VBOX,BOX,.FALSE.)
   53 CONTINUE
      CALL PGEND

C New histogram

      WRITE(*,*)'New histogram,same PA? (y/n) :'
      READ(*,*) Y
      IF (Y.EQ.'Y'.OR.Y.EQ.'y') GOTO 212
      GOTO 213
  999 STOP
      END
