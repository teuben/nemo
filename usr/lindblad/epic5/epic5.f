       MODULE PARAMS
       IMPLICIT NONE
       INTEGER, PARAMETER :: NRAD=1000
       END MODULE PARAMS

       PROGRAM EPIC5_fin
******************************************************************************
C     This program computes positions, velocities and densities along closed
C     orbits of interstellar matter, including frictional forces, in a galaxy 
C     with an arbitrary perturbing potential. Radial velocities are given for 
C     chosen lines of sight.
C     For theory and description see Per A.B. Lindblad, thesis Stockholm 
C     Observatory 1996, Chapt. 2.2 and Paper V (Lindblad & Lindblad 1994, 
C     ASPC 66, p.29), and the "EPIC5 manual".
C
C     Input parameters must be given in a parameter file with the following
C     content (file names with max 10 characters):
C
C         Name of file with rotation curve                    ('VFILE')
C         Name of file with perturbing potential              ('BFILE')
C         Pattern velocity in km/s/kpc                        (OMEGAP)
C         Coefficient of friction at 0 kpc                    (BLAMBD0)
C         Coefficient of friction at RMAX kpc                 (BLAMBDf)
C         Soft corotation coefficient in km/s/kpc             (EOMEG)
C         Minimum value for R in the computations in kpc      (R0)
C         Interval of R in the computations in kpc            (DR)
C         Maximum value for R in the computations in kpc      (RMAX)
C         Number of THETA values (<60) for each R             (NTHETA)
C
C     The rotation curve and the perturbing potential are read by
C     the subroutine ALLPLOT.
C
C     'VFILE' must have the following format:
C
C           NFRAD DRFRAD
C           FRAD(I), (I=1,NFRAD)
C
C         where NFRAD (<NRAD) is the number of velocities and DRFRAD is the 
C         increment between velocities in the rotation curve, and FRAD the 
C         velocities. I=1 corresponds to R = 0.
C
C     'BFILE' must have the following format:
C
C           NPER DRP MMAX NTH NC
C           AMPC(M), (M=1,10) 
C           AMPS(M), (M=1,10)
C           For M=1,MMAX:
C           CP(I,M) (I=1,NPER)
C           SP(I,M) (I=1,NPER)
C
C         where NPER (<NRAD) is the number of values and DRP the increment 
C         between values for the perturbing potential. CP(I,M) and SP(I,M)  
C         are amplitudes of the COS(M*THETA) resp. SIN(M*THETA) components of 
C         of the Fourier development the perturbing potential. AMPC(M) and 
C         AMPS(M) are factors with which these amplitudes will be multiplied. 
C         MMAX (<10) is the maximum M value in the development of the potential. 
C         NTH (<100) and NC are the number of THETA values and contours for 
C         potential plotting.
C
C     Positions and densities will be listed in file EPIDEN and plotted in  
C     file density.ps
C     Positions, x and y velocities will be listed in file EPIRV and plotted 
C     in file velocity.ps
C     Test values are written in file EPITEST
C     For printing the last input plot and orbit plot do (in a separate terminal 
C     window):
C         gv pgplot.ps and print the desired plots
C
C     EPIC first version written by P O Lindblad & P A B Lindblad 1993.08.11
C     This code has later been modified by P O Lindblad, N. Pi–ol-Ferrer, and K Fathi
C
C     Major modification history includes
C     EPIC5 version 2011.05.05 --> write in EPIRV VX and VY
C     EPIC5 version 2011.05.16 --> plotting densities over perturbing potential
C     EPIC5 version 2011.05.17 --> introduced parameter NRAD
C     EPIC5 version 2011.05.29 --> plotting omega vs kappa and omega/kappa vs radius (to be improved)
C     EPIC5 version 2011.06.06 --> plotting the density wave wavelenght vs radius, overplotting 4*dr, which is the allowed limit for the wavelenght.
C     EPIC5 version 2011.06.16 --> making possible plotting in ps format with colour and without names and dates and time (removing CALL PGIDEN)
C     EPIC5 version 2011.06.29 --> Change the scales on axis in plot omega vs kappa. 
C     EPIC5 version 2011.07.11 --> Change in order to allow to rotate in counterclock and clock - wise. 
C     EPIC5 version 2011.07.11 --> Remove plot of wavelenght of spiral arms if there is no spiral in the potential.
C     EPIC5 version 2011.09.07 --> Introduce of the parameter BLAMBDA=LAMBDA/KAPPA
C     EPIC5 version 2011.09.08 --> Correct LAMBDA definition in order to take into account a clock-wise rotation. 
C     EPIC5 version 2011.09.08 --> Remove the question about the direction in the rotation.
C     EPIC5 version 2011.09.09 --> Introduce the linear interpolation for BLAMBDA.
C     EPIC5 version 2011.09.29 --> Change of units in plot OMEGAT-KAPPAT/2 vs R. Change of xrange and yrange in plot KAPPAT vs
C                                  OMEGAT. Add of points in the orbit plot. Choose to plot the damping coeffiencient LAMBDA vs 
C                                  RADIUS before plotting orbits.
C     EPIC5 version 2011.10.03 --> Change of the dependency of the damping coefficient. Now it is just linear, introducing its values at OILR and OLR.
C     EPIC5 version 2012.02.01 --> Remove plot of density wave wavelenght versus radia from input plots. Final major modification.
C     EPIC5 version 2021.04.20 --> Compiled for gfortran v9
       
C
C***************************************************************************
      USE PARAMS
      IMPLICIT NONE


      INTEGER I,J,K,SYMBOL,BETAI,NC,IDIM,JDIM,I1,I2,J1,J2,II,M
      REAL RMAX,BETA,XLL,XUR,YLL,YUR,C(NRAD),CI,RADVT(NRAD,61)
      CHARACTER VFILE*10,BFILE*10,RADVEL*24,Y*1,EPIPAR*12
      EXTERNAL PLOTDEN
      EXTERNAL PLOTRV

      REAL XIA,ETAA,VXIA,VETAA,R1,R2,TH1,TH2

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1 
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP
      COMMON/TEST/R1,R2,TH1,TH2

      REAL CR(5),IILR(5),OILR(5),OLR(5)
      COMMON/RESON/CR,IILR,OILR,OLR

      REAL CPOT(50),APOT(NRAD,100)
      REAL DC,AMAX,DTH
      INTEGER NC2,NPERP
      COMMON/POTENTPLOT/APOT
      COMMON/POTP/DTH
      integer c1, c2, kk
      REAL BAXIS,AAXIS
      REAL LAMBDAV(NRAD)
      COMMON/DAMP/ LAMBDAV

      EXTERNAL PLOTPOT

      PI=3.14159265

      WRITE(*,*) '*****************************************************'
      WRITE(*,*) '** This program computes positions, velocities and **'
      WRITE(*,*) '** densities along closed orbits of interstellar   **'
      WRITE(*,*) '** matter,including frictional forces, in an       **'
      WRITE(*,*) '** arbitrary perturbing potential. Radial velo-    **'
      WRITE(*,*) '** cities are given for chosen lines of sight.     **'
      WRITE(*,*) '** The direction of rotation is                    **'
      WRITE(*,*) '** defined according to the sign of the pattern    **'
      WRITE(*,*) '** speed (positive: counter clock wise rotation,   **'
      WRITE(*,*) '** negative: clock wise rotation).                 **'
      WRITE(*,*) '** EPIC5 Version 2012.02.01                        **'
      WRITE(*,*) '*****************************************************'

C Read parameter file

 200  WRITE(*,*) 'Parameter file: '
      READ(*,201) EPIPAR
 201  FORMAT(A)

  61  OPEN(21,FILE=EPIPAR,STATUS='OLD',ERR=990)
      READ(21,FMT=*,ERR=991)VFILE,BFILE,OMEGAP,BLAMBD0,BLAMBDF,EOMEG,R0,
     &                      DR,RMAX,NTHETA
      GOTO 1
  990 WRITE(*,*)'CANNOT OPEN FILE', EPIPAR,'TRY AGAIN (y/n): '
      READ(*,20)Y
      IF (Y.EQ.'Y'.OR.Y.EQ.'y') GOTO 200
      GOTO 999
  991 WRITE(*,*)'ERROR WHEN READING FILE ',EPIPAR
      CLOSE(21)
      GOTO 999
   1  WRITE(*,*) 'File with rotation curve:       ',VFILE
      WRITE(*,*) 'File with perturbing potential: ',BFILE
      WRITE(*,21) ' Pattern velocity:     ',OMEGAP, ' km/s/kpc'
      WRITE(*,'(A,F7.2,F7.2,A)') ' Friction coefficients at OILR and OLR 
     &: ',BLAMBD0,BLAMBDF, ' km/s/kpc'
      WRITE(*,21) ' Corotation softening: ',EOMEG, ' km/s/kpc'
      WRITE(*,22) ' Min, delta, max R:',R0 ,DR, RMAX,' kpc '
      WRITE(*,*) 'Number of THETA values: ',NTHETA
   20 FORMAT(A1)
   21 FORMAT(A,F7.2,A)
   22 FORMAT(A,3F6.2,A)
      WRITE(*,*) 'ARE PARAMETERS CORRECT ? (y/n):'
      READ(*,20)Y
      IF (Y.EQ.'Y'.OR.Y.EQ.'y') GOTO 60
      WRITE(*,*) 'New parameter file ? (y/n): '
      READ(*,20)Y
      IF (Y.NE.'Y'.AND.Y.NE.'y') GOTO 999
      CLOSE(21)
      GOTO 200

 60   CLOSE(21)

C Open testfile

      OPEN(25,FILE='EPITEST',STATUS='OLD',ERR=989)
      R1 = 0.9
      R2 = 1.0
      TH1 = PI/4
      TH2 = TH1+PI/16
      GOTO 202
  989 WRITE(*,*)'CANNOT OPEN FILE EPITEST'
      GOTO 999

C Read rotation curve and potential. Optional plotting:

 202  RPMAX = RMAX

      CALL ALLPLOT(VFILE,BFILE)

      IF(ERROR.EQ.'y') GOTO 999
      IF(RPMAX.EQ.0.) THEN
        WRITE(*,94) 'Change RMAX: ',RMAX,' ? (y/n)'
 94     FORMAT(A,F6.1,A)
        READ(*,20) Y
        IF(Y.NE.'y'.AND.Y.NE.'Y') GOTO 999
        WRITE(*,*) 'RMAX: '
        READ(*,*) RMAX
        GOTO 202
      ENDIF
      
C Should orbits be computed or stop here?

      WRITE(*,*)'COMPUTE ORBITS ? (y/n):'
      READ(*,20)Y
      IF (Y.NE.'Y'.AND.Y.NE.'y') GOTO 999

  15  CONTINUE

C Compute output tables:

      IF(RMAX.LE.R0) THEN
        WRITE(*,*) 'RMAX less than R0'
        GOTO 999
      ENDIF
      NR=INT((RMAX+0.0001-R0)/DR)+1

      CALL COTABS

C Call plot routines here to plot orbits:

      XUR=1.1*RMAX
      XLL=-XUR
      YUR=XUR
      YLL=XLL
      K=1
!      WRITE(*,*)'SEE DAMPING COEFFICIENT vs RADIUS? (y/n)'
!      READ(*,20) Y
!      IF (Y.EQ.'Y'.OR.Y.EQ.'y') THEN 
!         CALL PGBEGIN(0,'/XWINDOW',1,1)
!         CALL PGSCF(2)  
!         CALL PGWNAD(0,XUR,MINVAL(LAMBDAV)-1,MAXVAL(LAMBDAV)+1)
!         CALL PGBOX('BCMNSTV',0.,0,'BCMNSTV',0.,0)
!         CALL PGLABEL('kpc','km s\u-1\d kpc\u-1\d','DAMPING COEFF.')
c        CALL PGIDEN
!         CALL PGSLS (1)
!         SYMBOL=17
!         CALL PGPT(1,0,LAMBDAV(1),SYMBOL)
!         DO I=2,NPER
!            IF (LAMBDAV(I).GT.0) THEN
!               CALL PGDRAW((I-1)*DRP,LAMBDAV(I))
!               CALL PGPT(1,(I-1)*DRP,LAMBDAV(I),SYMBOL)
!            END IF
!         END DO
!         CALL PGEND
!      END IF
      CALL PGBEGIN(0,'/XWINDOW',1,1)
 12   CALL PGSCF(2)
      CALL PGWNAD(XLL,XUR,YLL,YUR)
      CALL PGBOX('BCMNSTV',0.,0,'BCMNSTV',0.,0)
      CALL PGLABEL('kpc','kpc','ORBITS')
c     CALL PGIDEN
      CALL PGSFS (2)
      SYMBOL=17
      DO I=1,NR
         CALL PGMOVE(XT(I,1),YT(I,1))
C         CALL PGPOINT(1,XT(I,1),YT(I,1),SYMBOL)
         DO J=2,NTHETA
            CALL PGDRAW(XT(I,J),YT(I,J))
C            CALL PGPOINT(1,XT(I,J),YT(I,J),SYMBOL)
         END DO
         CALL PGDRAW(XT(I,1),YT(I,1))
      END DO
      call pgsci(2)
      CALL PGCIRC(0.,0.,IILR(2))
      CALL PGCIRC(0.,0.,IILR(1))
      CALL PGCIRC(0.,0.,OILR(1))
      CALL PGCIRC(0.,0.,OILR(2))
      CALL PGCIRC(0.,0.,CR(1))
      CALL PGCIRC(0.,0.,OLR(1))
      call pgsci(1)
      CALL PGEND
      
      IF (K.EQ.2) GOTO 11

C Print the orbit?

      WRITE(*,*)'PRINT THIS PLOT ? (y/n):'
      READ(*,13)Y
   13 FORMAT(A1)
      IF (Y.NE.'Y'.AND.Y.NE.'y') GOTO 11
      K=2
      CALL PGBEGIN(0,'orbits.ps/CPS',1,1)
      GOTO 12

   11 CONTINUE

C Change RMAX?

      WRITE(*,105) 'CHANGE RMAX (',RMAX,')? (y/n): '
 105  FORMAT(A,F6.1,A)
      READ(*,13) Y
      IF (Y.NE.'Y'.AND.Y.NE.'y') GOTO 14
      WRITE(*,*) 'RMAX:'
      READ(*,*) RMAX
      GOTO 15

C Enlarge arrays in theta

  14  DO 31 I=1,NR
        XT(I,NTHETA+1)=XT(I,1)
        YT(I,NTHETA+1)=YT(I,1)
        SIGMAT(I,NTHETA+1)=SIGMAT(I,1)
   31 CONTINUE
      

C Print output tables for positions and surface densities

      OPEN(22,FILE='EPIDEN',STATUS='OLD',ERR=992)
      WRITE(22,100)((XT(I,J),YT(I,J),SIGMAT(I,J),J=1,NTHETA),
     &I=1,NR-2)
  100 FORMAT(' ',2F8.3,F12.4)
      GOTO 9
  992 WRITE(*,*)'CANNOT OPEN FILE EPIDEN'
      GOTO 999
    9 CONTINUE
      CLOSE(22)

C Plot surface densities

   35 WRITE(*,*)'PLOT DENSITIES ? (y/n):'
      READ(*,34)Y
   34 FORMAT(A1)
      IF (Y.NE.'Y'.AND.Y.NE.'y') GOTO 30
      WRITE(*,*)'NUMBER OF CONTOURS?:'
      READ(*,FMT=*,ERR=302)NC
      WRITE(*,*)'CONTOURS?:'
      READ(*,FMT=*,ERR=302)(C(I),I=1,NC)
      GOTO 301
 302  WRITE(*,*) 'Wrong input'
      GOTO 35
      
 301  BAXIS=MINVAL(sqrt(XT(NR,1:NTHETA)**2+YT(NR,1:NTHETA)**2))
      AAXIS=MAXVAL(sqrt(XT(NR,1:NTHETA)**2+YT(NR,1:NTHETA)**2))
 
      XUR=1.1*AAXIS/SQRT(2.)
      YUR=1.1*BAXIS/SQRT(2.)
      XLL=-XUR
      YLL=-YUR
 
      DO 33 K=1,2
	IF (K.EQ.1) THEN
            CALL PGBEGIN(0,'/XWINDOW',1,1)
            ELSE
           CALL PGBEGIN(0,'density.ps/CPS',1,1)
        ENDIF
        CALL PGSCF(2)
        CALL PGWNAD(XLL,XUR,YLL,YUR)
        CALL PGBOX('BCMNSTV',0.,0,'BCMNSTV',0.,0)
        CALL PGLABEL('kpc','kpc','DENSITIES')
c        CALL PGIDEN
      CALL PGSFS (2)
      IDIM=NRAD
      JDIM=41
      I1=1
      I2=NR-2
      J1=1
      J2=NTHETA+1
      CALL PGCONX(SIGMAT,IDIM,JDIM,I1,I2,J1,J2,C,NC,PLOTDEN)
      call pgsci(2)
      CALL PGCIRC(0.,0.,IILR(2))
      CALL PGCIRC(0.,0.,IILR(1))
      CALL PGCIRC(0.,0.,OILR(1))
      CALL PGCIRC(0.,0.,OILR(2))
      CALL PGCIRC(0.,0.,CR(1))
      CALL PGCIRC(0.,0.,OLR(1))
      call pgsci(1)
      CALL PGEND

   33 CONTINUE


      WRITE(*,*)'PLOT DENSITIES WITH PERTURBING POTENTIAL? (y/n):'
      READ(*,34)Y
      IF (Y.NE.'Y'.AND.Y.NE.'y') GOTO 35
      CALL PGBEGIN(0,'/XWINDOW',1,1)
      AMAX=MAXVAL(APOT)
      DC = AMAX/10
      DO KK=1,10
        CPOT(KK) = KK*DC
        CPOT(10+KK) = -KK*DC
      ENDDO
      NC2 = 2*10
c      CALL PGIDEN
      CALL PGSCF(2)
      CALL PGSFS (2)
      CALL PGWNAD(XLL,XUR,YLL,YUR)
      CALL PGBOX('BCMNSTV',0.,0,'BCMNSTV',0.,0)
      CALL PGLABEL('kpc','kpc','PERTURBING POTENTIAL')
      IDIM = NRAD
      JDIM = 100
      I1 = 1
      I2 = INT(RPMAX/DRP)+1!NPERP
      J1 = 1
      J2 = NTHETA*2+1
      CALL PGCONX(APOT,IDIM,JDIM,I1,I2,J1,J2,CPOT,NC2,PLOTPOT)
      IDIM=NRAD
      JDIM=41
      I1=1
      I2=NR-2
      J1=1
      J2=NTHETA+1
      call pgsci(4)
      CALL PGCONX(SIGMAT,IDIM,JDIM,I1,I2,J1,J2,C,NC,PLOTDEN)
      call pgsci(2)
      CALL PGCIRC(0.,0.,IILR(2))
      CALL PGCIRC(0.,0.,IILR(1))
      CALL PGCIRC(0.,0.,OILR(1))
      CALL PGCIRC(0.,0.,OILR(2))
      CALL PGCIRC(0.,0.,CR(1))
      CALL PGCIRC(0.,0.,OLR(1))
      call pgsci(1)
      CALL PGEND



      GOTO 35
   30 CONTINUE
      
C Compute radial velocities:

    2 WRITE(*,*)'GIVE POSITION ANGLE OF LINE OF NODES (BAR=90,TO QUIT
     &WRITE 360 OR MORE):'
      READ(*,*)BETAI
      IF (BETAI.GE.360) GOTO 999
      WRITE(RADVEL(1:21),'(21A)')'RADIAL VELOCITIES PA='
      WRITE(RADVEL(22:24),'(I3.3)')BETAI
      BETA=BETAI*PI/180
      DO 3 I=1,NR
           DO 4 J=1,NTHETA
                RADVT(I,J)= VXT(I,J)*COS(BETA)+VYT(I,J)*SIN(BETA)
    4      CONTINUE
    3 CONTINUE

C Print output file for radial velocities:

      OPEN(23,FILE='EPIRV',STATUS='OLD',ERR=994)
      WRITE(23,101)((XT(I,J),YT(I,J),VXT(I,J),VYT(I,J),J=1,NTHETA),
     &I=1,NR)
  101 FORMAT(' ',2F8.3,2F12.1)
      GOTO 10
  994 WRITE(*,*)'CANNOT OPEN FILE EPIRV'
      GOTO 999
   10 CONTINUE
      CLOSE(23)


C Enlarge arrays in theta

      DO 51 I=1,NR
        RADVT(I,NTHETA+1)=RADVT(I,1)
   51 CONTINUE

C Plot radial velocities

      WRITE(*,*)'CONTOUR INTERVAL (KM/S)?:'
      READ(*,*)CI
      NC=2*INT(400/CI)+1
      K=(NC-1)/2+1
      C(K)=0.
      DO 55 J=1,NC/2
        C(K-J)=-J*CI
        C(K+J)=J*CI
   55 CONTINUE

      DO 53 K=1,2
	IF (K.EQ.1) THEN
            CALL PGBEGIN(0,'/XWINDOW',1,1)
            ELSE
           CALL PGBEGIN(0,'velocity.ps/CPS',1,1)
        ENDIF
        CALL PGSCF(2)
C	CALL PGSCI(5)
        CALL PGWNAD(XLL,XUR,YLL,YUR)
        CALL PGBOX('BCMNSTV',0.,0,'BCMNSTV',0.,0)
C	CALL PGSCI(8)
        CALL PGLABEL('kpc','kpc','LINE OF SIGHT VELOCITIES')
c        CALL PGIDEN
C	CALL PGSCI(13)
      CALL PGSFS (2)

      IDIM=NRAD
      JDIM=41
      I1=1
      I2=NR-2
      J1=1
      J2=NTHETA+1
      CALL PGCONX(RADVT,IDIM,JDIM,I1,I2,J1,J2,C,NC,PLOTRV)
      IF (k.EQ.1.) THEN 
         call pgsci(4)
         call pgconx(SIGMAT,IDIM,JDIM,I1,I2,J1,J2,[1.4,1.9],2,PLOTDEN)
         call pgsci(1)
      ENDIF
      call pgsci(2)
      CALL PGCIRC(0.,0.,IILR(2))
      CALL PGCIRC(0.,0.,IILR(1))
      CALL PGCIRC(0.,0.,OILR(1))
      CALL PGCIRC(0.,0.,OILR(2))
      CALL PGCIRC(0.,0.,CR(1))
      CALL PGCIRC(0.,0.,OLR(1))
      call pgsci(1)
      CALL PGEND
   53 CONTINUE
   50 CONTINUE

      GOTO 2
  999 CLOSE(25)
      STOP
      END



      SUBROUTINE PREPAR
C***********************************************************************
C     Prepares some input data tables
C***********************************************************************
      USE PARAMS
      IMPLICIT NONE


      INTEGER M,I
      REAL R,DOM

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBD0,BLAMBDF,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,EOMEG
      CHARACTER ERROR*1,SYM*5,PP*1 
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP
      CHARACTER COUNT*1
      COMMON/ROTATION/COUNT
      



C Create tables OMEGAT,KAPPAT,AT

      OMEGAT(1)=(-0.5*FRAD(3)+2*FRAD(2)-1.5*FRAD(1))/DRFRAD
      DO 2 I=2,NFRAD
          R=(I-1)*DRFRAD
          OMEGAT(I)=FRAD(I)/R
    2 CONTINUE

      DO 3 I=1,NFRAD-2
          R=(I-1)*DRFRAD
          DOM=(OMEGAT(I+1)-OMEGAT(I)-(OMEGAT(I+2)-2*OMEGAT(I+1)+
     &                                           OMEGAT(I))/2)/DRFRAD
          KAPPAT(I)=SQRT(2*OMEGAT(I)*(2*OMEGAT(I)+R*DOM))
          AT(I)=-R*DOM/2
    3 CONTINUE

      IF (COUNT.EQ.'Y'.OR.COUNT.EQ.'y') THEN 
         DO I=1,NFRAD
            KAPPAT(I)=-KAPPAT(I)
         ENDDO
      END IF
 

C Create tables DT,CT,ET:

      DO M=1,MMAX
          DT(1,M)=0.
      ENDDO
      DO M=1,MMAX
        DO I=2,NPER-2
          R=(I-1)*DRP
          DT(I,M) = M*PSI(I,M)/R
        ENDDO
      ENDDO

      DO M=1,MMAX
        DO I=1,NPER-2
          CT(I,M)=(PSI(I+1,M)-PSI(I,M)-(PSI(I+2,M)-2*PSI(I+1,M)+
     &                                             PSI(I,M))/2)/DRP
        ENDDO
      ENDDO
      
  5   DO M=1, MMAX
        DO I=1,NPER-2
          ET(I,M) = PSI(I,M)*(THM(I+1,M)-THM(I,M))/DRP
        ENDDO
      ENDDO

  6   RETURN
      END


      SUBROUTINE COTABS
C***********************************************************************
C     Computes output tables
C***********************************************************************
      USE PARAMS 
      IMPLICIT NONE



      INTEGER I,J,K,M,NTH
      REAL R,V,DXI,DETA,THETA,DTHETA,OMEGA,KAPPA,A,R1,R2,TH1,TH2
 
      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP
      COMMON/ROT/OMEGA,KAPPA,A
      COMMON/TEST/R1,R2,TH1,TH2
      REAL LAMBDAV(NRAD)
      COMMON/DAMP/ LAMBDAV


C Compute tables for xi, eta, vxi, veta, x, y, vx, vy

      R=R0
      DTHETA=2*PI/NTHETA
      DO 1 I=1,NR
          CALL INTERP(R,FRAD,V)
          CALL INTERP(R,OMEGAT,OMEGA)
          CALL INTERP(R,KAPPAT,KAPPA)
          CALL INTERP(R,AT,A)



C Test
      IF(R.GE.R1.AND.R.LE.R2) THEN
         WRITE(25,300) 'R=',R,' V=',V,' OMEGA=',OMEGA,' KAPPA=',KAPPA,
     &               ' A=',A
  300    FORMAT(A,F10.7,A,F8.4,A,F8.4,A,F8.4,A,F8.4)
         WRITE(25,*) ' '
      ENDIF

  6     THETA=0
        DO 2 J=1,NTHETA

          CALL COORD(R,V,THETA,XIT(I,J),ETAT(I,J),VXIT(I,J),
     &             VETAT(I,J),XT(I,J),YT(I,J),VXT(I,J),VYT(I,J))

          THETA=THETA+DTHETA
    2   CONTINUE

C Add two extra terms in theta

        DO 5 J=1,2
          K=NTHETA+J
          XIT(I,K)=XIT(I,J)
          ETAT(I,K)=ETAT(I,J)
    5   CONTINUE
        R=R+DR
    1 CONTINUE

C Compute table for sigma

      R=R0
      K=NR-2
      DO 3 I=1,K
          DO 4 J=1,NTHETA
             THETA = (J-1)*DTHETA
             DXI=(XIT(I+1,J)-XIT(I,J)-(XIT(I+2,J)-2*XIT(I+1,J)+
     &            XIT(I,J))/2)/DR
             DETA=(ETAT(I,J+1)-ETAT(I,J)-(ETAT(I,J+2)-2*ETAT(I,J+1)+
     &            ETAT(I,J))/2)/DTHETA
             IF (R.EQ.0) THEN
                   SIGMAT(I,J)=0
                   ELSE
                   SIGMAT(I,J)=1-XIT(I,J)/R-DETA/R-DXI

C  Test
      IF(R.GE.R1.AND.R.LE.R2.AND.THETA.GE.TH1.AND.THETA.LE.TH2) THEN
        WRITE(25,305) 'R=',R,' THETA=',THETA,' I=',I,' J=',J,
     &  ' SIGMAT(I,J)=',SIGMAT(I,J)
  305   FORMAT(A,F10.7,A,F8.6,A,I3,A,I3,A,F8.5)
        WRITE(25,*) ' '
      ENDIF

             ENDIF
    4     CONTINUE
          R=R+DR
    3 CONTINUE
      RETURN
      END      


      SUBROUTINE INTERP(R,TAB,VALUE)
C***********************************************************************
C     Interpolates within a VFILE to a given R
C***********************************************************************      
      USE PARAMS
      IMPLICIT NONE


      INTEGER I,NFRAD,NPER,MMAX
      REAL R,TAB(NRAD),VALUE,X,OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)

      I=INT(R/DRFRAD)+1
      X=R/DRFRAD-INT(R/DRFRAD)
      VALUE=TAB(I)+X*(TAB(I+1)-TAB(I))+X*(X-1)*(TAB(I+2)-2.*TAB(I+1)+
     &      TAB(I))/2.
      RETURN
      END


      SUBROUTINE INTPOT(R,M,TAB,VALUE)
C***********************************************************************
C     Interpolates within two-dimensional potential tables to a given R
C***********************************************************************

      USE PARAMS
      IMPLICIT NONE

      INTEGER I,M,NFRAD,NPER,MMAX
      REAL R,TAB(NRAD,3),VALUE,X,OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
    
      I=INT(R/DRP)+1
      X=R/DRP-INT(R/DRP)
      VALUE=TAB(I,M)+X*(TAB(I+1,M)-TAB(I,M))+X*(X-1)*(TAB(I+2,M)-
     &      2.*TAB(I+1,M)+TAB(I,M))/2.
      RETURN
      END


      SUBROUTINE AMPL(R,M)
C******************************************************************************
C     Computes the amplitudes of the deviations from circular motion
C******************************************************************************
      USE PARAMS
      IMPLICIT NONE


      INTEGER M,NTH,II
      REAL R,OMEGA,KAPPA,KNY,KNYM,F,E,G,H,P,Q,C,D,I,J,K,L,T,EKNYM,A
     &,OM,OMS,DELTAR

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1
      REAL LAMBDA
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP
      COMMON/ROT/OMEGA,KAPPA,A
      COMMON/POT/C,D,E
      REAL LAMBDAV(NRAD)
      COMMON/DAMP/ LAMBDAV

      REAL CR(5),IILR(5),OILR(5),OLR(5)
      COMMON/RESON/CR,IILR,OILR,OLR


      CALL INTPOT(R,M,CT,C)
      CALL INTPOT(R,M,DT,D)
      II=INT(R/DRP)+1
      DELTAR=R/DRP-INT(R/DRP)
      BLAMBDA=BLAMBD0+(BLAMBDF-BLAMBD0)/
     &(MAXVAL(OLR)-MAXVAL(OILR))*((II-1)*DRP-MAXVAL(OILR))
      LAMBDA=BLAMBDA!*abs(OMEGA)!KAPPA
      LAMBDAV(II)=LAMBDA
      E=ET(II,M)+DELTAR*(ET(II+1,M)-ET(II,M))

      OM = M*(OMEGA - OMEGAP)
      OMS = OM/(OM**2 + (M*EOMEG)**2)
      
      H = KAPPA**2 - OM**2
      G = OM*C + 2*OMEGA*D
      K = KAPPA**2 + OM**2
      L = OM*C - 2*OMEGA*D
      P = OMEGA - 2*A
      Q = 3*KAPPA**2 + OM**2 - 8*OMEGA**2
      T = H**2 + 8*LAMBDA**2*K + 16*LAMBDA**4

      DV = (H*G-2*LAMBDA*K*E-4*LAMBDA**2*L-8*LAMBDA**3*E)/T
      GV = ((-H*(2*OMEGA*G-H*D)+4*LAMBDA*(H*A+2*OM**2*OMEGA)*E
     &           -4*LAMBDA**2*(2*P*G-Q*D)+16*LAMBDA**3*A*E)/T)*OMS
      EV = ((OM**2*H*E+2*LAMBDA*(K*G-2*H*OMEGA*D)-4*LAMBDA**2*OM**2*E
     &           +8*LAMBDA**3*OM*C)/T)*OMS     
      FV = ((2*OM*H*OMEGA*E+2*LAMBDA*(4*OM*OMEGA*G+H*(2*A*C-OM*D))
     &      +8*LAMBDA**2*OM*P*E+8*LAMBDA**3*(2*A*C+OM*D))/T)*OMS
      DE = DV*OMS
      GE = GV*OMS
      EE = EV*OMS
      FE = FV*OMS

      RETURN
      END


      SUBROUTINE COORD(R,V,THETA,XI,ETA,VXI,VETA,X,Y,VX,VY)
C***************************************************************************
C     Computes coordinates xi, eta, x, y and corresponding velocities
C***************************************************************************

      USE PARAMS
      IMPLICIT NONE

      INTEGER M,NTH,I
      REAL V,R,THETA,X,Y,VX,VY,RR,THETA1,XI,ETA,VXI,VETA,TH,TH1,TH2
      REAL C,D,E,R1,R2,DELTAR

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1 
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP
      COMMON/POT/C,D,E
      COMMON/TEST/R1,R2,TH1,TH2
      REAL LAMBDAV(NRAD)
      COMMON/DAMP/ LAMBDAV


C Compute xi, eta, vxi, veta:

      XI=0
      ETA=0
      VXI=0
      VETA=0
      DO M=1,MMAX
        CALL AMPL(R,M)

C Test
        IF(R.GE.R1.AND.R.LE.R2.AND.THETA.GE.TH1.AND.THETA.LE.TH2) THEN
           WRITE(25,301) 'R=',R,' M=',M
 301       FORMAT(A,F10.7,A,I2)
           WRITE(25,*) 'C=',C,' D=',D,' E=',E
           WRITE(25,*) 'DE=',DE,' GE=',GE,' EE=',EE,' FE=',FE
           WRITE(25,*) ' '
        ENDIF

        I=INT(R/DRP)+1
        DELTAR=R/DRP-INT(R/DRP)
        TH=THM(I,M)+DELTAR*(THM(I+1,M)-THM(I,M))
        XI=XI+DE*COS(M*THETA-TH)+EE*SIN(M*THETA-TH)
        ETA=ETA+GE*SIN(M*THETA-TH)+FE*COS(M*THETA-TH)
        VXI=VXI-DV*SIN(M*THETA-TH)+EV*COS(M*THETA-TH)
        VETA=VETA+GV*COS(M*THETA-TH)-FV*SIN(M*THETA-TH)
      ENDDO

C Test
      IF(R.GE.R1.AND.R.LE.R2.AND.THETA.GE.TH1.AND.THETA.LE.TH2) THEN
        WRITE(25,303) 'R=',R,' THETA=',THETA
  303   FORMAT(A,F10.7,A,F8.6)
        WRITE(25,*) 'XI=',XI,' ETA=',ETA,' VXI=',VXI,' VETA=',VETA
        WRITE(25,*) ' '
      ENDIF

C Compute x, y, vx, vy:

      RR=R+XI
      THETA1=THETA+ETA/R
      Y=RR*SIN(THETA1)
      X=RR*COS(THETA1)
      VY=VXI*SIN(THETA1)+(V+VETA)*COS(THETA1)*RR/R
      VX=VXI*COS(THETA1)-(V+VETA)*SIN(THETA1)*RR/R

C Test
      IF(R.GE.R1.AND.R.LE.R2.AND.THETA.GE.TH1.AND.THETA.LE.TH2) THEN
        WRITE(25,304) 'R=',R,' THETA=',THETA
  304   FORMAT(A,F10.7,A,F8.6)
        WRITE(25,*) 'X=',X,' Y=',Y,' VX=',VX,' VY=',VY
        WRITE(25,*) ' '
      ENDIF

      RETURN
      END


      SUBROUTINE PLOTDEN(VISBLE,X,Y,Z)
C******************************************************************************
C     This is a subroutine used by PGCONX to convert from a rectangular grid
C     to a grid made up by the computed points along the closed orbits.
C******************************************************************************
      USE PARAMS
      IMPLICIT NONE


      INTEGER VISBLE,I,J,NTH
      REAL X,Y,ID,JD,X1,X2,Y1,Y2,Z,XWORLD,YWORLD,DTHETA

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1  
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP

      I=INT(X)
      ID=X-INT(X)
      J=INT(Y)
      JD=Y-INT(Y)
      X1=XT(I,J)+ID*(XT(I+1,J)-XT(I,J))
      X2=XT(I,J+1)+ID*(XT(I+1,J+1)-XT(I,J+1))
      XWORLD=X1+JD*(X2-X1)
      Y1=YT(I,J)+ID*(YT(I+1,J)-YT(I,J))
      Y2=YT(I,J+1)+ID*(YT(I+1,J+1)-YT(I,J+1))
      YWORLD=Y1+JD*(Y2-Y1)

      IF (Z.EQ.1.) CALL PGSLS(4)
      IF (Z.GT.1.) CALL PGSLS(1)
      IF (Z.LT.1.) CALL PGSLS(2)

      IF (VISBLE.EQ.0) THEN
	CALL PGMOVE(XWORLD,YWORLD)
        ELSE
	CALL PGDRAW(XWORLD,YWORLD)
      END IF
      END

      SUBROUTINE PLOTRV(VISBLE,X,Y,Z)
C******************************************************************************
C     This is a subroutine used by PGCONX to convert from a rectangular grid
C     to a grid made up by the computed points along the closed orbits.
C******************************************************************************
      USE PARAMS
      IMPLICIT NONE


      INTEGER VISBLE,I,J,NTH
      REAL X,Y,ID,JD,X1,X2,Y1,Y2,Z,XWORLD,YWORLD,DTHETA

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1  
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP

      I=INT(X)
      ID=X-INT(X)
      J=INT(Y)
      JD=Y-INT(Y)
      X1=XT(I,J)+ID*(XT(I+1,J)-XT(I,J))
      X2=XT(I,J+1)+ID*(XT(I+1,J+1)-XT(I,J+1))
      XWORLD=X1+JD*(X2-X1)
      Y1=YT(I,J)+ID*(YT(I+1,J)-YT(I,J))
      Y2=YT(I,J+1)+ID*(YT(I+1,J+1)-YT(I,J+1))
      YWORLD=Y1+JD*(Y2-Y1)

      IF (Z.EQ.0.) CALL PGSLS(4)
      IF (Z.GT.0.) CALL PGSLS(1)
      IF (Z.LT.0.) CALL PGSLS(1)

      IF (VISBLE.EQ.0) THEN
	CALL PGMOVE(XWORLD,YWORLD)
        ELSE
	CALL PGDRAW(XWORLD,YWORLD)
      END IF
      END


      SUBROUTINE ALLPLOT(VFILE,BFILE)
C********************************************************************
C*** This subroutine reads given input files, and plots the 
C*** potentials, rotation curve, resonance positions,radial and 
C*** tangential forces.                                       
C*** First version of this subroutine was written by P.A.B. Lindblad in 1992-05-03                  
C********************************************************************
      USE PARAMS
      IMPLICIT NONE


      INTEGER NBTYPE(4),N
      REAL OMPA,AMPC(10),AMPS(10),REDFAC(4),CP(NRAD,10),SP(NRAD,10)
      CHARACTER VFILE*10,BFILE*10,Y*1

      COMMON/POTP/DTH

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1  
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP

      INTEGER NTH,NC,K,I,J,I1,I2,IDIM,J1,J2,JDIM,NCT,
     &        NPERP,M,NC2
      REAL APOT(NRAD,100),DTH,AMAX,CPOT(50),XUR,XLL,YUR,YLL
      REAL DC,R,PHI(NRAD),V,TPOT(NRAD,100),TMAX,RM,RMP

      REAL CR(5),IILR(5),OILR(5),OLR(5)
      COMMON/RESON/CR,IILR,OILR,OLR

      COMMON/POTENTPLOT/APOT
      CHARACTER COUNT*1
      COMMON/ROTATION/COUNT

      EXTERNAL PLOTPOT

C Read the number of points, NFRAD, and the increment between
C points, DRFRAD, for the rotation curve 

      ERROR = 'n'
      OPEN(27,FILE=VFILE,STATUS='OLD',ERR=801)
      GOTO 2
 801  WRITE(*,*) 'Cannot open file ',VFILE
      ERROR = 'y'
      GOTO 997 
 2    READ(27, *) NFRAD, DRFRAD

C Check if desired R-range is permitted

      RM = (NFRAD-3)*DRFRAD
      IF(RM.LT.RPMAX) THEN
       WRITE(*,*) 'RMAX TOO LARGE FOR ROTATION CURVE TABLE'
       RPMAX = 0.
       GOTO 997
      ENDIF

C Read the rotation curve

      READ(27, *) (FRAD(I) , I = 1, NFRAD)

C      WRITE(*,*) 'Clock-wise rotation?(y/n)'
C      READ(*,*) COUNT
      IF (OMEGAP.LT.0) COUNT='y'
      IF (COUNT.EQ.'Y'.OR.COUNT.EQ.'y') THEN 
         DO I=1,NFRAD
            FRAD(I)=-FRAD(I)
         ENDDO
      END IF
      
 

C If a perturbing potential,read its parameters 

      DO M=1,10
        DO I=1,NRAD
          CP(I,M) = 0. 
          SP(I,M) = 0.
          APOT(I,M) = 0.
        ENDDO
        AMPC(M) = 0.
        AMPS(M) = 0.
      ENDDO
      AMAX = 0.
      WRITE(*,*) 'Perturbing potential ? (y/n) :'
      READ(*,91)PP
      IF(PP.NE.'y'.AND.PP.NE.'Y') THEN
        NPER=NFRAD
        DRP=DRFRAD
        MMAX = 10
        NTH = 64
        DTH = 2*PI/NTH
        NC = 5
        GOTO 8
      ENDIF

 700  OPEN(26,FILE=BFILE,STATUS='OLD',ERR=802)
      GOTO 3
 802  WRITE(*,*) 'Cannot open file ',BFILE
      ERROR = 'y'
      GOTO 996 
 3    READ(26, *) NPER, DRP, MMAX, NTH, NC

C Check if desired R-range is permitted

      RMP = (NPER-3)*DRP
      IF(RMP.LT.RPMAX) THEN
        WRITE(*,*) 'RMAX TOO LARGE FOR POTENTIAL TABLE'
        RPMAX = 0.
        GOTO 996
      ENDIF

C Read the perturbing potential 

      READ(26,*) (AMPC(M),M=1,10)
      READ(26,*) (AMPS(M),M=1,10)
      DO M=1,MMAX
        READ(26,*) (CP(I,M),I=1,NPER)
        READ(26,*) (SP(I,M),I=1,NPER)
      ENDDO

C Scale the amplitude of the disturbing potential with the scaling
C factors AMP(M).

  5   DO M=1,MMAX
        DO I=1,NPER
          CP(I,M) = CP(I,M) * AMPC(M)
          SP(I,M) = SP(I,M) * AMPS(M)
        ENDDO
      ENDDO

C Compute data tables

  10  DO M=1,MMAX
        DO I=1,NPER
          PSI(I,M) = SQRT(CP(I,M)**2+SP(I,M)**2)
        ENDDO
        N=0
        DO I=2,NPER
          IF(CP(I,M).EQ.0..AND.SP(I,M).EQ.0.) THEN
            THM(I,M) = 0.
            ELSE 
            THM(I,M) = ATAN2(SP(I,M),CP(I,M))
          ENDIF
           THM(I,M) = THM(I,M) + N*2*PI
          IF(SP(I-1,M).GT.0..AND.SP(I,M).LT.0.) THEN
            IF(CP(I,M).GT.0.) GOTO 100
            THM(I,M) = THM(I,M) + 2*PI
            N=N+1
  100       CONTINUE
          ENDIF
          IF(SP(I-1,M).LT.0..AND.SP(I,M).GT.0.) THEN
            IF(CP(I,M).GT.0.) GOTO 101
            THM(I,M) = THM(I,M) - 2*PI
            N=N-1
  101       CONTINUE
          ENDIF
        ENDDO
        THM(1,M)=THM(2,M)
      ENDDO      

C Create a file with the perturbing potential

      DTH = 2*PI/NTH
      DO I=1,NPER
        DO J=1,(NTH+1)
          APOT(I,J) = 0.
        ENDDO
      ENDDO
      IF(PP.NE.'y'.AND.PP.NE.'Y') GOTO 8
      DO I=1,NPER
        DO J=1,(NTH+1)
          DO M=1,MMAX 
            APOT(I,J) = APOT(I,J) + PSI(I,M)*COS(M*J*DTH-THM(I,M))
          ENDDO
          IF (APOT(I,J).GT.AMAX) THEN
            AMAX = APOT(I,J)
          END IF
        ENDDO
      ENDDO

C Create a file with the total potential

  8   TMAX = 0.
      NCT = 4*NC
      NPERP = INT(RPMAX/DRP)
      R = NPERP*DRP
      CALL INTERP(R,FRAD,V)
      PHI(NPERP) = V**2*DRP/R
      DO I=1,(NPERP-1)
        R = (NPERP-I)*DRP
        CALL INTERP(R,FRAD,V)
        PHI(NPERP-I) = PHI(NPERP+1-I) + V**2*DRP/R
      ENDDO

      DO I=1,NPERP
        DO J=1,(NTH+1)
          TPOT(I,J) = PHI(I) + APOT(I,J)
          IF (TPOT(I,J).GT.TMAX) THEN
            TMAX = TPOT(I,J)
          ENDIF
        ENDDO
      ENDDO

      CALL PREPAR
      CALL RESONANCE

C Should the input be plotted?

 4    WRITE(*,*)'SEE THE INPUT PLOTS ? (y/n):'
      READ(*,91)Y
 91   FORMAT(A1)
      IF (Y.NE.'Y'.AND.Y.NE.'y') GOTO 92
      WRITE(*,94) 'RMAX = ',RPMAX,' ? (y/n): '
 94   FORMAT(A,F6.1,A)
      READ(*,91)Y
      IF (Y.EQ.'Y'.OR.Y.EQ.'y') GOTO 93
      WRITE(*,*) 'RMAX: '
      READ (*,*) RPMAX

C Start the plotting session

  93  IF(AMAX.EQ.0.) GOTO 6
      IF(PP.NE.'y'.AND.PP.NE.'Y') GOTO 6

C Plot the perturbing potential

      WRITE(*,*) 'PERTURBING POTENTIAL:'
      CALL PGBEGIN(0,'?',1,1)
      DC = AMAX/NC
      DO K=1,NC
        CPOT(K) = K*DC
        CPOT(NC+K) = -K*DC
      ENDDO
      NC2 = 2*NC
c      CALL PGIDEN
      CALL PGSCF(2)
      CALL PGSFS (2)
      XUR = 1.1*RPMAX
      XLL = -XUR
      YUR = XUR
      YLL = XLL
      CALL PGWNAD(XLL,XUR,YLL,YUR)
      CALL PGBOX('BCMNSTV',0.,0,'BCMNSTV',0.,0)
      CALL PGLABEL('kpc','kpc','PERTURBING POTENTIAL')
      IDIM = NRAD
      JDIM = 100
      I1 = 1
      I2 = NPERP
      J1 = 1
      J2 = NTH+1
      CALL PGCONX(APOT,IDIM,JDIM,I1,I2,J1,J2,CPOT,NC2,PLOTPOT)
      call pgsci(2)
      CALL PGCIRC(0.,0.,IILR(2))
      CALL PGCIRC(0.,0.,OILR(1))
      CALL PGCIRC(0.,0.,OILR(2))
      CALL PGCIRC(0.,0.,CR(1))
      CALL PGCIRC(0.,0.,OLR(1))
      call pgsci(1)
      CALL PGEND

C Plot the total potential

  6   DC = TMAX/NCT
      DO K=1,NCT
        CPOT(K) = K*DC
      ENDDO
      WRITE(*,*) 'TOTAL POTENTIAL:'
      CALL PGBEGIN(0,'?',1,1)
c      CALL PGIDEN
      CALL PGSCF(2)
      CALL PGSFS (2)
      XUR = 1.1*RPMAX
      XLL = -XUR
      YUR = XUR
      YLL = XLL
      CALL PGWNAD(XLL,XUR,YLL,YUR)
      CALL PGBOX('BCMNSTV',0.,0,'BCMNSTV',0.,0)
      CALL PGLABEL('kpc','kpc','TOTAL POTENTIAL')
      IDIM = NRAD
      JDIM = 100
      I1 = 1
      I2 = NPERP
      J1 = 1
      J2 = NTH+1
      CALL PGCONX(TPOT,IDIM,JDIM,I1,I2,J1,J2,CPOT,NCT,PLOTPOT)
      call pgsci(2)
      CALL PGCIRC(0.,0.,IILR(2))
      CALL PGCIRC(0.,0.,OILR(1))
      CALL PGCIRC(0.,0.,OILR(2))
      CALL PGCIRC(0.,0.,CR(1))
      CALL PGCIRC(0.,0.,OLR(1))
      call pgsci(1)
      CALL PGEND

      WRITE(*,*) 'REMAINING INPUT PLOTS ? (y/n):'
      READ(*,91) Y
      IF(Y.NE.'y'.AND.Y.NE.'Y') GOTO 996
      CALL PGBEGIN(0,'?',1,1)

C Plot the resonances

      CALL RESPLOT

C Plot the bar potential components

      CALL BARPOTC
      
C Plot the radial forces

      CALL RADFOR

C Plot the tangential forces

      CALL TANGFOR

C End the plotting session

      CALL PGEND
      GOTO 4

   92 CONTINUE
 996  CLOSE (26)
 997  CLOSE (27)

      RETURN
      END

      SUBROUTINE PLOTPOT(VISBLE,X,Y,Z)
C***********************************************************************
C** This routine plots the potential in 2D, used by PGPLOT
C***********************************************************************

      USE PARAMS
      IMPLICIT NONE


      INTEGER VISBLE
      INTEGER NPER,MMAX
      REAL X,Y,Z,XWORLD,YWORLD,DRP,DTH
      REAL FRAD,PSI,THM
     

      COMMON/POTP/DTH
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)

      XWORLD = X*DRP*COS(Y*DTH)
      YWORLD = X*DRP*SIN(Y*DTH)

C Test
 !      WRITE(*,*) X,Y,Z,DRP,DTH,NPER,VISBLE
 !      WRITE(*,*) XWORLD,YWORLD

      IF (VISBLE.EQ.0) THEN
        CALL PGMOVE(XWORLD,YWORLD)
      ELSE
        CALL PGDRAW(XWORLD,YWORLD)
      ENDIF
      END


      SUBROUTINE RESONANCE
C***********************************************************************
C** This routine computes the positions of the resonances.           ***
C** First version of this subroutine was written by P.A.B. Lindblad in 1992-03-11                 ***
C***********************************************************************

      USE PARAMS
      IMPLICIT NONE


      INTEGER I,NR1,NR2,NR3,NR4,M
      REAL R(NRAD),OMMKH(NRAD),OMPKH(NRAD),THMM(NRAD,10),THMN(NRAD,10)
      REAL XMAX,YMAX,YMIN,DOMMKH(NRAD),DOMPKH(NRAD)
      REAL OMMKF(NRAD),OMPKF(NRAD),DOMMKF(NRAD),DOMPKF(NRAD),P,X

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1  
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP
      
      REAL CR(5),IILR(5),OILR(5),OLR(5)
      COMMON/RESON/CR,IILR,OILR,OLR
      COMMON/KAPPA/OMMKH,OMPKH,OMMKF,OMPKF,NR1,NR2,NR3,NR4
      CHARACTER COUNT*1
      COMMON/ROTATION/COUNT

      DO I = 1,NFRAD
            R(I) = (I-1)*DRFRAD
      ENDDO

C Compute OMEGA-K/m and OMEGA+K/m (m=2,4) and their derivatives

      DO  I=1,NFRAD-2
            OMMKH(I) = OMEGAT(I) - KAPPAT(I)/2.
            OMPKH(I) = OMEGAT(I) + KAPPAT(I)/2.
            OMMKF(I) = OMEGAT(I) - KAPPAT(I)/4.
            OMPKF(I) = OMEGAT(I) + KAPPAT(I)/4.
      ENDDO


C Search for resonances

      DO I=1,5
        CR(I)=0.
        IILR(I)=0.
        OILR(I)=0.
        OLR(I)=0.
      ENDDO

      NR1 = 0
      NR2 = 0
      NR3 = 0
      NR4 = 0

      IF (COUNT.EQ.'Y'.OR.COUNT.EQ.'y') THEN 
         DO I=1,NFRAD-1
            IF(OMEGAT(I).GE.OMEGAP.AND.OMEGAT(I+1).LT.OMEGAP) THEN
                  NR1 = NR1+1
                  CR(NR1) = R(I)+DRFRAD*(OMEGAP-OMEGAT(I))/(OMEGAT(I+1)
     &                                                      -OMEGAT(I))
            ENDIF

            IF(OMEGAT(I).LE.OMEGAP.AND.OMEGAT(I+1).GT.OMEGAP) THEN
                  NR1 = NR1+1
                  CR(NR1) = R(I)+DRFRAD*(OMEGAP-OMEGAT(I))/(OMEGAT(I+1)
     &                                                      -OMEGAT(I))
            ENDIF
            IF(OMMKH(I).GE.OMEGAP.AND.OMMKH(I+1).LT.OMEGAP) THEN
               NR2 = NR2+1
               P = ABS(OMEGAP-OMMKH(I))
               X = ABS(OMMKH(I)-OMMKH(I+1))
               IILR(NR2) = R(I) + DRFRAD*P/X
            ENDIF
            
            IF(OMMKH(I).LE.OMEGAP.AND.OMMKH(I+1).GT.OMEGAP) THEN
               NR3 = NR3+1
               P = ABS(OMEGAP-OMMKH(I))
               X = ABS(OMMKH(I)-OMMKH(I+1))
               OILR(NR3) = R(I) + DRFRAD*P/X
            ENDIF

            IF(OMPKH(I).LE.OMEGAP.AND.OMPKH(I+1).GT.OMEGAP) THEN
               NR4 = NR4+1
               P = ABS(OMEGAP-OMPKH(I))
               X = ABS(OMPKH(I)-OMPKH(I+1))
               OLR(NR4) = R(I) + DRFRAD*P/X
            ENDIF
         ENDDO
      ELSE    
 
         DO I =1,NFRAD-1

            IF(OMEGAT(I).GE.OMEGAP.AND.OMEGAT(I+1).LT.OMEGAP) THEN
                  NR1 = NR1+1
                  CR(NR1) = R(I)+DRFRAD*(OMEGAP-OMEGAT(I))/(OMEGAT(I+1)
     &                                                      -OMEGAT(I))
            ENDIF

            IF(OMEGAT(I).LE.OMEGAP.AND.OMEGAT(I+1).GT.OMEGAP) THEN
                  NR1 = NR1+1
                  CR(NR1) = R(I)+DRFRAD*(OMEGAP-OMEGAT(I))/(OMEGAT(I+1)
     &                                                      -OMEGAT(I))
            ENDIF

            IF(OMMKH(I).LE.OMEGAP.AND.OMMKH(I+1).GT.OMEGAP) THEN
                  NR2 = NR2+1
                  P = ABS(OMEGAP-OMMKH(I))
                  X = ABS(OMMKH(I)-OMMKH(I+1))
                  IILR(NR2) = R(I) + DRFRAD*P/X
            ENDIF

            IF(OMMKH(I).GE.OMEGAP.AND.OMMKH(I+1).LT.OMEGAP) THEN
                  NR3 = NR3+1
                  P = ABS(OMEGAP-OMMKH(I))
                  X = ABS(OMMKH(I)-OMMKH(I+1))
                  OILR(NR3) = R(I) + DRFRAD*P/X
            ENDIF

            IF(OMPKH(I).GE.OMEGAP.AND.OMPKH(I+1).LT.OMEGAP) THEN
                  NR4 = NR4+1
                  P = ABS(OMEGAP-OMPKH(I))
                  X = ABS(OMPKH(I)-OMPKH(I+1))
                  OLR(NR4) = R(I) + DRFRAD*P/X
            ENDIF

         ENDDO
      END IF


      IF(NR2.EQ.0) GOTO 42
      DO I=1,NR2
        WRITE(*,41) 'IILR :',IILR(I),' kpc'
      ENDDO
  42  IF(NR3.EQ.0) GOTO 43
      DO I=1,NR3
        WRITE(*,41) 'OILR :',OILR(I),' kpc'
      ENDDO
  43  IF(NR1.EQ.0) GOTO 40
      DO I=1,NR1
        WRITE(*,41) 'CR :  ',CR(I),' kpc'
  41    FORMAT(A,F6.2,A)
      ENDDO
  40  IF(NR4.EQ.0) GOTO 44
      DO I=1,NR4
        WRITE(*,41) 'OLR : ',OLR(I),' kpc'
      ENDDO
  44  CONTINUE
      RETURN
      END

      SUBROUTINE RESPLOT
C***********************************************************************
C** This routine plots the positions of the resonances.              ***
C** First version of this subroutine was written by P.A.B. Lindblad in 1992-03-11     ***
C***********************************************************************

      USE PARAMS
      IMPLICIT NONE


      INTEGER I,NR1,NR2,NR3,NR4,M
      REAL R(NRAD),OMMKH(NRAD),OMPKH(NRAD),THMM(NRAD,10),THMN(NRAD,10)
      REAL XMAX,YMAX,YMIN,DOMMKH(NRAD),DOMPKH(NRAD)
      REAL OMMKF(NRAD),OMPKF(NRAD),DOMMKF(NRAD),DOMPKF(NRAD),P,X

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1  
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP
      
      REAL CR(5),IILR(5),OILR(5),OLR(5)
      COMMON/RESON/CR,IILR,OILR,OLR
      COMMON/KAPPA/OMMKH,OMPKH,OMMKF,OMPKF,NR1,NR2,NR3,NR4
      REAL WAVEL(NRAD)


      DO I = 1,NFRAD
            R(I) = (I-1)*DRFRAD
      ENDDO

C Plot the curves. Start with setting the scales

      XMAX = RPMAX
      YMAX = abs(OMEGAP) * 20.

      CALL PGENV(0.,XMAX,0.,YMAX,0,0)
c      CALL PGIDEN
      CALL PGSCF(2)
 11   CALL PGLABEL('Radius (kpc)','Omega (km/s/kpc)','RESONANCES')

C Plot OMEGA(R)

 14   CALL PGSLS(1)
      CALL PGMOVE(R(1),abs(OMEGAT(1)))
      DO I = 2,NFRAD-2
            CALL PGDRAW(R(I),abs(OMEGAT(I)))
      ENDDO

C Plot OMMKH(R)

      CALL PGSLS(2)
      CALL PGMOVE(R(1),abs(OMMKH(1)))
      DO I = 2,NFRAD-2
            CALL PGDRAW(R(I),abs(OMMKH(I)))
      ENDDO

C Plot OMPKH(R)

      CALL PGSLS(3)
      CALL PGMOVE(R(1),abs(OMPKH(1)))
      DO I = 2,NFRAD-2
            CALL PGDRAW(R(I),abs(OMPKH(I)))
      ENDDO

C Plot OMPKF(R)

      CALL PGSLS(4)
      CALL PGMOVE(R(1),abs(OMPKF(1)))
      DO I = 2,NFRAD-2
            CALL PGDRAW(R(I),abs(OMPKF(I)))
      ENDDO

C Plot OMMKF(R)

      CALL PGSLS(5)
      CALL PGMOVE(R(1),abs(OMMKF(1)))
      DO I = 2,NFRAD-2
            CALL PGDRAW(R(I),abs(OMMKF(I)))
      ENDDO

C Plot the constant value of OMEGAP

      CALL PGSLS(1)
      CALL PGMOVE(0.,abs(OMEGAP))
      CALL PGDRAW(R(NFRAD-2),abs(OMEGAP))

C Mark the position of the resonances

      IF(NR1.EQ.0) GOTO 110
      DO I = 1,NR1
            CALL PGMOVE(CR(I),abs(OMEGAP))
            CALL PGDRAW(CR(I),0.)
      ENDDO

 110  IF(NR2.EQ.0) GOTO 120
      DO I = 1,NR2
            CALL PGMOVE(IILR(I),abs(OMEGAP))
            CALL PGDRAW(IILR(I),0.)
      ENDDO

 120  IF(NR3.EQ.0) GOTO 130
      DO I = 1,NR3
            CALL PGMOVE(OILR(I),abs(OMEGAP))
            CALL PGDRAW(OILR(I),0.)
      ENDDO

 130  IF(NR4.EQ.0) GOTO 140
      DO I = 1,NR4
            CALL PGMOVE(OLR(I),abs(OMEGAP))
            CALL PGDRAW(OLR(I),0.)
      ENDDO

 140  CONTINUE

C Plot the rotation curve

      XMAX = RPMAX
      YMAX = maxval(abs(frad))*1.1
      CALL PGENV(0.,XMAX,0.,YMAX,0,0)
 21   CALL PGLABEL('Radius (kpc)',
     &   'Velocity (km/s)','ROTATION CURVE')
c 24   !CALL PGIDEN 
      CALL PGMOVE(R(1),abs(FRAD(1)))
      DO I = 2,NFRAD
            CALL PGDRAW(R(I),abs(FRAD(I)))
      ENDDO

c     PLOTTING OMEGA VS KAPPA 
      YMAX=MAXVAL(abs(OMEGAT-KAPPAT/2.))
      CALL PGENV(0.,XMAX,0.,YMAX,0,0)
      CALL PGLABEL('Radius (kpc)',
     &   '\(0550)-\(0636)/2 (km s\u-1\d kpc\u-1\d)','')
c      CALL PGIDEN 
      CALL PGMOVE(R(1),abs(OMEGAT(1)-(KAPPAT(1)/2.)))
      DO I = 2,NFRAD-2
            CALL PGDRAW(R(I),abs((OMEGAT(I)-(KAPPAT(I)/2.))))
      ENDDO

c     PLOTTING OMEGA VS KAPPA 
      XMAX=MAXVAL(abs(OMEGAT))
      YMAX=MAXVAL(abs(KAPPAT))
      YMIN=MINVAL(ABS(KAPPAT))   
      CALL PGENV(0,XMAX,YMIN,YMAX,0,0)
      CALL PGLABEL('\(0550) (km s\u-1\d kpc\u-1\d)',
     &   '\(0636) (km s\u-1\d kpc\u-1\d) ','')
c      CALL PGIDEN 
      CALL PGMOVE(abs(OMEGAT(1)),abs(KAPPAT(1)))
      DO I = 2,NFRAD-2
            CALL PGDRAW(abs(OMEGAT(I)),abs(KAPPAT(I)))
      ENDDO
      call pgsls(2)
      CALL PGLINE(NFRAD,abs(OMEGAT),2.*abs((OMEGAT-OMEGAP)))
      call pgmtxt('TV',-4.,0.1,
     &0.,'--- Constant: 2*(\(0550)-\(0550)\dp\u)')
      call pgsls(1)


C Plot the THM
      DO M=1,MMAX
         THMN(:,M)=THM(:,M)/M
      ENDDO
      XMAX = RPMAX
      YMAX = MAXVAL(THMN(:,1:MMAX))+0.1
      YMIN = MINVAL(THMN(1:INT(RPMAX/DRP)+1,1:MMAX))-.1
      CALL PGENV(0.,XMAX,YMIN,YMAX,0,0)
 31   CALL PGLABEL('Radius (kpc)',
     &   'PHASE SHIFT (rad)','PHASE SHIFT')
      R(1)=0.
      DO M=1, MMAX     
c 34   !CALL PGIDEN 
      CALL PGMOVE(R(1),THMN(1,M))
        DO I = 2,NPER
           R(I)=(I-1)*DRP
           CALL PGDRAW((I-1)*DRP,THMN(I,M))
        ENDDO
      ENDDO

      CALL PGPT(nper,R,THMN(:,2),-3)
      
C Plotting the wavelenght of the density wave     
      DO I=1,NPER-1
         WAVEL(I)=PI*DRP/(THM(I+1,2)-THM(I,2))
         IF (ABS(THM(I+1,2)-THM(I,2)).LT.0.0001) WAVEL(I)=0.
      ENDDO
      YMAX=maxval(wavel)
      YMIN=minval(wavel)
      IF (YMAX-YMIN.LT.0.5) GOTO 22
C      CALL PGENV(0.,XMAX,YMIN,YMAX,0,0)
C      CALL PGLABEL('Radius (kpc)',
C     &   '\(0637) (kpc)','Wavelenght')
      R(1)=0.
c      CALL PGIDEN 
C      CALL PGMOVE(R(1),WAVEL(1))
C      DO I = 1,NPER
 
C         CALL PGDRAW((I-1)*DRP,WAVEL(I))
C      ENDDO
      
C      CALL PGMOVE(R(1),4.*DRP)
C      DO I = 1,NPER
 
C         CALL PGDRAW((I-1)*DRP,4.*DRP)
C      ENDDO
      
C      CALL PGMOVE(R(1),-4.*DRP)
C      DO I = 1,NPER
 
c         CALL PGDRAW((I-1)*DRP,-4.*DRP)
c      ENDDO
      


 22   RETURN
      END


      SUBROUTINE BARPOTC
C****************************************************************
C*** This routine plots the perturbing potential components.  ***
C****************************************************************

      USE PARAMS
      IMPLICIT NONE


      INTEGER I,M
      REAL RMAX,C2P,RADI(NRAD),PSIMAX

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1  
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP

C Compute the maximum radius, and the radius array

       RMAX = RPMAX
       DO I = 1,NPER
	  RADI(I) = (I-1)*DRP
       ENDDO
      PSIMAX=PSI(1,1)
      DO M=1,MMAX
        DO I=2,NPER
          IF(PSIMAX.GE.PSI(I,M)) GOTO 30
          PSIMAX=PSI(I,M)
  30      CONTINUE
        ENDDO
      ENDDO
      PSIMAX=1.1*PSIMAX

C Plot the components

      CALL PGENV(0.,RMAX,0.,PSIMAX,0,0)
c      CALL PGIDEN
      CALL PGSCF(2)
 31   CALL PGLABEL('Radius (kpc)',
     &   'Potential','PSI amplitudes successive M values')

 34   CALL PGSLS(1)
      DO M=1,MMAX
        CALL PGMOVE(RADI(1),PSI(1,M))
        DO I = 2,NPER
	  CALL PGDRAW(RADI(I),PSI(I,M))
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE RADFOR
C******************************************************************
C** This subroutine plots the axisymmetric radial force from the***
C** rotation curve. RADFO=FRAD**2/R, and the radial forces due  ***
C** to the perturbing fourier components.                       ***
C** First version of this subroutine was written by P.A.B. Lindblad in 1992-05-01              ***
C******************************************************************

      USE PARAMS
      IMPLICIT NONE


      INTEGER I,J,M
      REAL R(NRAD),RAD(NRAD),RADFO(NRAD,10),CP,RFO(NRAD),FACT(5)
      REAL XMAX,YMAX,YMIN

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1  
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP

      DO I = 1,NFRAD
            R(I) = (I-1)*DRFRAD
      ENDDO
      DO I = 1,NPER
        RAD(I) = (I-1)*DRP
      ENDDO

C Compute the radial axisym.force from the rotation curve

      XMAX = RPMAX
      YMIN = 0.
      YMAX = 0.
      DO I = 1,NFRAD-1
	RFO(I) = -FRAD(I+1)*FRAD(I+1)/R(I+1)
      ENDDO

C Compute the radial forces

      DO M=1,MMAX
        DO I=1,NPER-2
          RADFO(I,M) = CT(I,M)*COS(THM(I,M)) 
     &               - ET(I,M)*SIN(THM(I,M))
        ENDDO
      ENDDO

C Plot the forces

      YMAX=0.
      DO M=1,MMAX
        DO I=2,NPER-2
          IF(YMIN.LT.RADFO(I,M)) GOTO 45
          YMIN=RADFO(I,M)
  45      CONTINUE
        ENDDO
      ENDDO
      YMAX=RADFO(I,1)
      DO M=1,MMAX
        DO I=2,NPER-2
          IF(YMAX.GT.RADFO(I,M)) GOTO 46
          YMAX=RADFO(I,M)
  46      CONTINUE
        ENDDO
      ENDDO
C      if (maxval(rfo).gt.ymax) ymax=maxval(rfo)
C      if (minval(rfo).lt.ymin) ymin=minval(rfo)

      YMAX = 1.1*YMAX
      YMIN = 1.1*YMIN
      CALL PGSLS(1)
      CALL PGSCF(2)
      CALL PGENV(0.,XMAX,YMIN,YMAX,0,0)
 41   CALL PGLABEL('Radius (kpc)',
     &             'Radial force ','RADIAL FORCES ALONG PA -90 DEGR')
c 44   !CALL PGIDEN 

C Start with the axisym. force

      CALL PGMOVE(R(2),RFO(1))

      DO I = 2,NFRAD-1
            CALL PGDRAW(R(I+1),RFO(I))
      ENDDO

C Continue with the perturbing components
      DO M=1,MMAX      
        CALL PGMOVE(RAD(1),RADFO(1,M))
        DO I = 2,NPER-2
	  CALL PGDRAW(RAD(I),RADFO(I,M))
        ENDDO
      ENDDO

C Also draw a zero-line

       CALL PGMOVE(0.,0.)
       CALL PGDRAW(XMAX,0.)
       
      RETURN
      END


      SUBROUTINE DERIV (A,A1,H,M,N)
C*************************************************************************
C*** Takes the derivatives of all points from M to N in the 
C*** array A, and put the derivatives in array A1. 
C*** H is table increment.
C*************************************************************************

      DIMENSION A(1),A1(1)


      H60=H*60.

C Take the derivatives in the first M+2 points in the 
C array A, and put the derivatives in array A1. 

      A1(M)=(60.*A(M+7)/7.-70.*A(M+6)+252.*A(M+5)-525.*A(M+4)+700.*A(M+3
     &)-630.*A(M+2)+420.*A(M+1)-155.5714*A(M))/H60

      A1(M+1)=(2.*A(M+6)-15.*A(M+5)+50.*A(M+4)-100.*A(M+3)+150.*A(M+2)-7
     &7.*A(M+1)-10.*A(M))/H60 

      A1(M+2)=(2.*A(M)-24.*A(M+1)-35.*A(M+2)+80.*A(M+3)-30.*A(M+4)+8.*A(
     &M+5)-A(M+6))/H60

      MM=M+3
      NN=N-3

C Take the derivatives at all points from M to N

      DO J=MM,NN
      A1(J)=(-A(J-3)+9.*(A(J-2)-A(J+2))-45.*(A(J-1)-A(J+1 ))+A(J+3))/H60
      ENDDO

C Take the derivatives of the last N - N+2 points in the array

      A1(N-2)=(-2.*A(N)+24.*A(N-1)+35.*A(N-2)-80.*A(N-3)+30.*A(N-4)-8.*A
     &(N-5)+A(N-6))/H60
 
      A1(N-1)=(-2.*A(N-6)+15.*A(N-5)-50.*A(N-4)+100.*A(N-3)-150.*A(N-2)+
     &77.*A(N-1)+10.*A(N))/H60

      A1(N)=(147.*A(N)-360.*A(N-1)+450.*A(N-2)-400.*A(N-3)+225.*A(N-4)-7
     &2.*A(N-5)+10.*A(N-6))/H60
 
      RETURN
      END 



      SUBROUTINE TANGFOR
C******************************************************************
C** This subroutine plots the tangential force from the tables  ***
C** containing the perturbing potential.                        ***
C** First version of this subroutine was written by P.A.B. Lindblad in 1992-05-01              ***
C******************************************************************

      USE PARAMS
      IMPLICIT NONE


      INTEGER I,J,M
      REAL RADI(NRAD),TANFO(NRAD,10)
      REAL CP,FACT(5)
      REAL XMAX,YMAX,YMIN

      INTEGER NFRAD,NPER,MMAX,NR,NTHETA 
      REAL OMEGAP,DRFRAD,RPMAX,FRAD,DRP,PSI,THM,PI,R0,DR,BLAMBDA,EOMEG
     &,OMEGAT,KAPPAT,AT,CT,DT,ET,XIT,ETAT,VXIT,VETAT,XT,YT,VXT,VYT
     &,SIGMAT,DE,GE,EE,FE,DV,GV,EV,FV,BLAMBD0,BLAMBDF
      CHARACTER ERROR*1,SYM*5,PP*1  
      COMMON/TOTAL/ OMEGAP,NFRAD,DRFRAD,RPMAX
      COMMON/POTAL/ FRAD(NRAD),NPER,DRP,MMAX,PSI(NRAD,10),THM(NRAD,10)
      COMMON/GEMEN/PI,R0,DR,NR,NTHETA,BLAMBD0,BLAMBDF,EOMEG,OMEGAT(NRAD)
     &,KAPPAT(NRAD),AT(NRAD),CT(NRAD,10),DT(NRAD,10),ET(NRAD,10)
     & ,XIT(NRAD,62),ETAT(NRAD,62)
     & ,VXIT(NRAD,60),VETAT(NRAD,60),XT(NRAD,61),YT(NRAD,61)
     &,VXT(NRAD,60),VYT(NRAD,60),SIGMAT(NRAD,61),DE,GE
     &,EE,FE,DV,GV,EV,FV
      COMMON/CHAR/ERROR,SYM,PP

C Compute the maximum tangential force 

      DO M=1,MMAX
        DO I = 1,NPER-2
	  RADI(I) = (I-1)*DRP
	  TANFO(I,M) = -DT(I,M)*SIN(M*PI/4-THM(I,M))
        ENDDO
      ENDDO

C Plot the forces

      XMAX = RPMAX
      YMIN = TANFO(1,1)
      DO M=1,MMAX
        DO I=2,NPER-2
          IF(YMIN.LT.TANFO(I,M)) GOTO 45
          YMIN=TANFO(I,M)
  45      CONTINUE
        ENDDO
      ENDDO
      YMAX=TANFO(1,1)
      DO M=1,MMAX
        DO I=2,NPER-2
          IF(YMAX.GT.TANFO(I,M)) GOTO 46
          YMAX=TANFO(I,M)
  46      CONTINUE
        ENDDO
      ENDDO

      YMAX = 1.1*YMAX
      YMIN = 1.1*YMIN
      CALL PGSLS(1)
      CALL PGSCF(2)
      CALL PGENV(0.,XMAX,YMIN,YMAX,0,0)
 51   CALL PGLABEL('Radius (kpc)',
     &  'Tangential force ','TANGENTIAL FORCES ALONG PA -45 DEGR')
c 54   !CALL PGIDEN 

      DO M=1,MMAX
        CALL PGMOVE(RADI(1),TANFO(1,M))
        DO I = 2,NPER-2
          CALL PGDRAW(RADI(I),TANFO(I,M))
        ENDDO
      ENDDO

C Also draw a zero-line

      CALL PGMOVE(0.,0.)
      CALL PGDRAW(XMAX,0.)

      RETURN
      END








