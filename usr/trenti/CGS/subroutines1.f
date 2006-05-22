CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CGS - Collisionless Galactic Simulator - Source Files  C  
C                                                         C
C  Created by M. Trenti & T.van Albada in Fortran77 2003  C
C                                                         C
C  Version 1.0 Alpha                                      C
C                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  PLEASE READ THE READ.ME FILE WHERE YOU CAN FIND 
C  TECNICAL COMMENTS ABOUT THE CODE AND INSTRUCTIONS TO
C  USE IT



C____________________________________________________
C CMASS: COMPUTES THE CENTER OF MASS OF THE SYSTEM USING 
C ALL PARTICLE INSIDE THE GRID

      SUBROUTINE CMASS()

      IMPLICIT NONE

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

C____________________________________________________
C COMMON FOR GRID QUANTITIES
      INCLUDE 'com_grid.blk'

C_____________________________________________________
C COMMON FOR CENTER OF MASS DISPLACEMENTS
      INCLUDE 'com_cmass.blk'

C____________________________________________________
C LOCAL VARIABLES

      REAL XCG,YCG,ZCG        !CARTESIAN COORDINATES OF CENTER OF MASS
      REAL VXCG,VYCG,VZCG     !CARTESIAN VELOCITIES OF CENTER OF MASS
      
      INTEGER NPINGR            !NUMBER OF PARTICLES INSIDE GRID
      INTEGER L                 !INDEX RUNNING OVER PARTICLES

C____________________________________________________
C INIT
      XCG=0.0
      YCG=0.0
      ZCG=0.0
      VXCG=0.0
      VYCG=0.0
      VZCG=0.0
      NPINGR=0


C____________________________________________________
C COMPUTE CM

      DO 130 L=1,NP             !LOOP OVER PARTICLES
         IF (S(L).GT.RR(IL)) GO TO 130
C____________________________________________________
C                  USE ONLY PARTICLES INSIDE GRID
         NPINGR=NPINGR+1
         XCG=XCG+X(L)
         YCG=YCG+Y(L)
         ZCG=ZCG+Z(L)
         VXCG=VXCG+VX(L)
         VYCG=VYCG+VY(L)
         VZCG=VZCG+VZ(L)
 130  CONTINUE

C_____________________________________________________
C NORMALIZATION
      XCG=XCG/FLOAT(NPINGR)
      YCG=YCG/FLOAT(NPINGR)
      ZCG=ZCG/FLOAT(NPINGR)
      VXCG=VXCG/FLOAT(NPINGR)
      VYCG=VYCG/FLOAT(NPINGR)
      VZCG=VZCG/FLOAT(NPINGR)

C______________________________________________________
C THE SYSTEM IS NOW RECENTERED:

      DO L=1,NP                 !LOOP OVER PARTICLES
         X(L)=X(L)-XCG
         Y(L)=Y(L)-YCG
         Z(L)=Z(L)-ZCG
         VX(L)=VX(L)-VXCG
         VY(L)=VY(L)-VYCG
         VZ(L)=VZ(L)-VZCG
         S(L)=SQRT(X(L)*X(L)+Y(L)*Y(L)+Z(L)*Z(L))
         TH(L)=ACOS(Z(L)/S(L))
         PH(L)=ATAN2(Y(L),X(L))
         IF (PH(L).LT.0.0) PH(L)=PH(L)+2.0*PI
      ENDDO

C__________________________________________________________
C the recentering is logged into /CMstore/ common block
      xCMstore(1)=xCMstore(1)+XCG
      xCMstore(2)=xCMstore(2)+YCG
      xCMstore(3)=xCMstore(3)+ZCG
      vCMstore(1)=vCMstore(1)+VXCG
      vCMstore(2)=vCMstore(2)+VYCG
      vCMstore(3)=vCMstore(3)+VZCG


      RETURN
      END    


C_________________________________________________________
C*********************************************************


C____________________________________________________
C CMASS: COMPUTES THE CENTER OF MASS OF THE SYSTEM USING 
C ALL PARTICLE INSIDE THE HALF MASS RADIUS OF THE SYSTEM

      SUBROUTINE CMASS1()

      IMPLICIT NONE

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

C____________________________________________________
C COMMON FOR GRID QUANTITIES
      INCLUDE 'com_grid.blk'

C_____________________________________________________
C COMMON FOR CENTER OF MASS DISPLACEMENTS
      INCLUDE 'com_cmass.blk'

C____________________________________________________
C LOCAL VARIABLES

      REAL XCG,YCG,ZCG        !CARTESIAN COORDINATES OF CENTER OF MASS
      REAL VXCG,VYCG,VZCG     !CARTESIAN VELOCITIES OF CENTER OF MASS
      
      INTEGER NPINGR            !NUMBER OF PARTICLES INSIDE GRID
      INTEGER L                 !INDEX RUNNING OVER PARTICLES

C____________________________________________________
C INIT
      XCG=0.0
      YCG=0.0
      ZCG=0.0
      VXCG=0.0
      VYCG=0.0
      VZCG=0.0
      NPINGR=0


C____________________________________________________
C COMPUTE CM

      DO 130 L=1,NP             !LOOP OVER PARTICLES
         IF (S(L).GT.rhm) GO TO 130
C____________________________________________________
C                  USE ONLY PARTICLES INSIDE GRID
         NPINGR=NPINGR+1
         XCG=XCG+X(L)
         YCG=YCG+Y(L)
         ZCG=ZCG+Z(L)
         VXCG=VXCG+VX(L)
         VYCG=VYCG+VY(L)
         VZCG=VZCG+VZ(L)
 130  CONTINUE

C_____________________________________________________
C NORMALIZATION
      XCG=XCG/FLOAT(NPINGR)
      YCG=YCG/FLOAT(NPINGR)
      ZCG=ZCG/FLOAT(NPINGR)
      VXCG=VXCG/FLOAT(NPINGR)
      VYCG=VYCG/FLOAT(NPINGR)
      VZCG=VZCG/FLOAT(NPINGR)

C______________________________________________________
C THE SYSTEM IS NOW RECENTERED:

      DO L=1,NP                 !LOOP OVER PARTICLES
         X(L)=X(L)-XCG
         Y(L)=Y(L)-YCG
         Z(L)=Z(L)-ZCG
         VX(L)=VX(L)-VXCG
         VY(L)=VY(L)-VYCG
         VZ(L)=VZ(L)-VZCG
         S(L)=SQRT(X(L)*X(L)+Y(L)*Y(L)+Z(L)*Z(L))
         TH(L)=ACOS(Z(L)/S(L))
         PH(L)=ATAN2(Y(L),X(L))
         IF (PH(L).LT.0.0) PH(L)=PH(L)+2.0*PI
      ENDDO


C__________________________________________________________
C the recentering is logged into /CMstore/ common block
      xCMstore(1)=xCMstore(1)+XCG
      xCMstore(2)=xCMstore(2)+YCG
      xCMstore(3)=xCMstore(3)+ZCG
      vCMstore(1)=vCMstore(1)+VXCG
      vCMstore(2)=vCMstore(2)+VYCG
      vCMstore(3)=vCMstore(3)+VZCG


      RETURN
      END    

C_________________________________________________________
C*********************************************************
C____________________________________________________
C CMASS: COMPUTES THE CENTER OF MASS OF THE SYSTEM USING 
C ALL PARTICLE INSIDE THE HALF MASS RADIUS OF THE SYSTEM

      SUBROUTINE CMASSlook()

      IMPLICIT NONE

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

C____________________________________________________
C COMMON FOR GRID QUANTITIES
      INCLUDE 'com_grid.blk'

C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C LOCAL VARIABLES

      REAL XCG,YCG,ZCG        !CARTESIAN COORDINATES OF CENTER OF MASS
      REAL VXCG,VYCG,VZCG     !CARTESIAN VELOCITIES OF CENTER OF MASS
      
      INTEGER NPINGR            !NUMBER OF PARTICLES INSIDE GRID
      INTEGER L                 !INDEX RUNNING OVER PARTICLES

C____________________________________________________
C INIT
      XCG=0.0
      YCG=0.0
      ZCG=0.0
      VXCG=0.0
      VYCG=0.0
      VZCG=0.0
      NPINGR=0


C____________________________________________________
C COMPUTE CM

      DO 130 L=1,NP             !LOOP OVER PARTICLES
         IF (S(L).GT.rhm) GO TO 130
C____________________________________________________
C                  USE ONLY PARTICLES INSIDE GRID
         NPINGR=NPINGR+1
         XCG=XCG+X(L)
         YCG=YCG+Y(L)
         ZCG=ZCG+Z(L)
         VXCG=VXCG+VX(L)
         VYCG=VYCG+VY(L)
         VZCG=VZCG+VZ(L)
 130  CONTINUE

C_____________________________________________________
C NORMALIZATION
      XCG=XCG/FLOAT(NPINGR)
      YCG=YCG/FLOAT(NPINGR)
      ZCG=ZCG/FLOAT(NPINGR)
      VXCG=VXCG/FLOAT(NPINGR)
      VYCG=VYCG/FLOAT(NPINGR)
      VZCG=VZCG/FLOAT(NPINGR)

C______________________________________________________
C write cmass informations
      write(45,456) time, XCG,YCG,ZCG,VXCG,VYCG,VZCG

 456  format(F10.5,1X,F12.6,1X,F12.6,1X,F12.6,1X,F12.6,1X,
     &     F12.6,1X,F12.6,1X)
      RETURN
      END    

C_________________________________________________________
C*********************************************************










C_________________________________________________________

C SPHERICAL HARMONICS POISSON SOLVER: COMPUTES THE COEFFICIENTS OF 
C EXPANSION FOR THE ACCELERATION IN SPHERICAL HARMONICS SOLVING THE 
C SELF CONSISTENT POISSON EQUATION
C__________________________________________________________
      SUBROUTINE POISSON() 

      IMPLICIT NONE
C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

C____________________________________________________
C COMMON FOR GRID QUANTITIES
      INCLUDE 'com_grid.blk'
 
C____________________________________________________
C SPHERICAL HARMONIC EXPANSION COEFFICIENTS
      INCLUDE 'com_har.blk'

C____________________________________________________
C LEGENDRE POLYNOMIALS
      INCLUDE 'com_lp.blk'


C____________________________________________________
C LOCAL VARIABLES


      double precision eps      !paramter of tolerance near cth=+/-1
      parameter(eps=1.d-7)

      real*8  cmphi(LMAX),smphi(LMAX) !COS(M*PHI) & SIN(M*PHI)
   
      real onePmass             !MASS OF SINGLE PARTICLE
      REAL VOLShell(IL)         !VOLUME OF EACH CELL
      
      integer IPart             !INDEX RUNNING OVER PARTICLES
      INTEGER L                 !INDEX RUNNING OVER COEFF. IN SPHE. HAR.
      INTEGER N,M               !INDEXES RUNNING OVER EXPANSION SPHE. HAR. 
C                               !N-1=HARMONIC USED, M-1="THETA INDEX" USED 
      
      INTEGER I                 !RADIAL CELL RUNNING INDEX
      INTEGER II                !RADIAL CELL INDEX OF EACH PARTICLE
      REAL*8 DI                 !NORMALIZED DISTANCE FROM THE CELL 

      INTEGER I0                !NEAREST RADIAL CELL
      INTEGER I1                !NEXT NEAREST RADIAL CELL

      DOUBLE PRECISION CTH      !COS(TH)
      
      REAL*8 AUX                !AUX VARIABLE TO SAVE COMPUTATION TIME
      REAL*8 AN                 !FLOAT(N)

C____________________________________________________
C INIT STUFF:

C______________________________
      onePmass=MTOT/float(NP)   !single particle mass

C____________________________________________________
      DO I=1,IL                 !SPHE. HAR. COEFF.
         L=0
         DO N=1,NHAR+1
            DO M=1,N
               L=L+1
               ANM(L,I)=0.
               BNM(L,I)=0.
            ENDDO
         ENDDO      
      ENDDO
C_____________________________________________________
      do i=1,IL1                !RADIAL SHELL VOLUME
         VOLShell(I)=(RR(I+1)**3-RR(I)**3)/3.
      enddo

C end of init stuff

C______________________________________________________
C                  COEFFICIENTS ANM AND BNM




C______________________________________________________
c loop over particles
      do IPart=1,NP

c check if the particle is out of the grid and eventually skip assignment
         if(IP(IPart).gt.IL) goto 17

C______________________________________________________
c generates legendre polynomials for the particle

         do N=1,lmax
            do M=1,lmax
               plm(N,M)=0.d0
            enddo
         enddo
         cth=Z(IPart)/S(IPart)  !COS THETA


C______________________________________________________
C  VERY EURISTIC BUT SOMEWHAT TESTED...
C  USED TO AVOID CTH=+/-1 

         if(1.-(dabs(cth)).le.eps) then
            if(cth.gt.0.) then
               cth=1.-eps
            else
               cth=-1.+eps
            endif
         else
         endif
         
C______________________________________________________
C LEGENDRE POLYNOMIALS FOR EACH PARTICLE
         call  plgndrMT1(NHAR,CTH) 


C______________________________________________________
C generate the array containing cos(m phi) and sin(m phi)
         do i=1,nhar+1
            cmphi(i)=cos((float(i-1)*PH(IPart)))
            smphi(i)=sin((float(i-1)*PH(IPart)))
         enddo
         
 
C______________________________________________________
C RADIAL CELL ASSIGNMENT (IPart == PARTICLE INDEX)
         II=IP(IPart)  
         DI=(S(IPart)-RR(II))/DR(II) !NORMALIZED DISTANCE FROM CELL


C_______________________________________________________
C CHOOSE TWO NEAREST CELLS TO ASSIGN MASS, I0 & I1
C VAN ALBADA RECIPE FOR ASSIGNMENT         
         IF (DI.GT.0.5) GO TO 11
         I0=II                  
         I1=II-1
c     check if the particle is at the center of the system
         if(I1.eq.0) I1=1

         GO TO 12
 11      DI=1-DI
         I0=II
         I1=II+1
c check if the particle is out of the grid and eventually skip assignment
         if(I1.gt.IL1) goto 17

 12      CONTINUE

         
C______________________________________________________________
C PARTICLE CONTRIBUTION ASSIGNMENT

C______________________________________________________________
c loop over spherical harmonics 
         L=0                    !RUNNING INDEX
         DO N=1,NHAR+1          
            DO M=1,N
               L=L+1

               aux=plm(N,M)*cmphi(M)
               ANM(L,I0)=ANM(L,I0)+aux*(.5+DI)/volshell(I0)
               ANM(L,I1)=ANM(L,I1)+aux*(.5-DI)/volshell(I1)
               
               aux=plm(N,M)*smphi(M)
               BNM(L,I0)=BNM(L,I0)+aux*(.5+DI)/volshell(I0)
               BNM(L,I1)=BNM(L,I1)+aux*(.5-DI)/volshell(I1)

            enddo
         enddo                  !end loop over spher.har.


c     end loop over particles
 17      continue               !skip assignment for out of grid particles
      enddo

C_________________________________________________________
C NORMALIZATION LOOP
      DO I=1,IL                 !LOOP OVER SPHERICAL GRID
      L=0
      DO N=1,NHAR+1             !LOOP OVER SPHE HAR
         DO M=1,N
            L=L+1

            ANM(L,I)=C(N,M)*ANM(L,I)*onePmass
            BNM(L,I)=C(N,M)*BNM(L,I)*onePmass

         enddo
      enddo
      enddo


C NOW I HAVE THE AMN & BMN and all goes on as usual from van Albada code

C_____________________________________________________
C                   COEFFICIENTS CNM AND DNM
C SEE TECNICAL PAPER FOR DETAILS: 

      L=0
      DO N=1,NHAR+1
         DO M=1,N
            L=L+1
            FFA(L,IL)=0.0
            FFB(L,IL)=0.0
            GGA(L,1)=0.0
            GGB(L,1)=0.0
            DO  I=1,IL1
               FFA(L,IL-I)=(FQ(N,IL-I)*DRR(IL-I))*ANM(L,IL-I)+
     1              QQ(N,IL-I)*FFA(L,IL+1-I)
               FFB(L,IL-I)=(FQ(N,IL-I)*DRR(IL-I))*BNM(L,IL-I)+
     1              QQ(N,IL-I)*FFB(L,IL+1-I)
      GGA(L,I+1)=(GQ(N,I)*DRR(I))*ANM(L,I)+(QQ(N,I)*QQ(2,I))*GGA(L,I)
      GGB(L,I+1)=(GQ(N,I)*DRR(I))*BNM(L,I)+(QQ(N,I)*QQ(2,I))*GGB(L,I)
            enddo
         enddo
      enddo

C
      DO I=1,IL
         L=0
         DO N=1,NHAR+1
            AN=FLOAT(N-1)
            DO M=1,N
               L=L+1
               CNM(L,I)=-(PI4G/(2.0*AN+1.0))*(FFA(L,I)+GGA(L,I))
               DNM(L,I)=-(PI4G/(2.0*AN+1.0))*(FFB(L,I)+GGB(L,I))
            enddo
         enddo
      enddo

C________________________________________________________________
C     END OF POISSON SOLVER

      RETURN
      END


C_________________________________________________________
C*********************************************************
C_________________________________________________________

C_________________________________________________________
C     ACCELERATION ASSIGNMENT


      SUBROUTINE accel()

      IMPLICIT NONE
C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

C____________________________________________________
C COMMON FOR GRID QUANTITIES
      INCLUDE 'com_grid.blk'
 
C____________________________________________________
C SPHERICAL HARMONIC EXPANSION COEFFICIENTS
      INCLUDE 'com_har.blk'

C____________________________________________________
C LEGENDRE POLYNOMIALS
      INCLUDE 'com_lp.blk'

C____________________________________________________
C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C LOCAL VARIABLES


      double precision eps      !paramter of tolerance near cth=+/-1
      parameter(eps=1.d-7)
      double precision eps1     
      parameter(eps1=5.d-5)     !parameter of tolerance near S(IPART)=0   

      real*8  cmphi(LMAX),smphi(LMAX) !COS(M*PHI) & SIN(M*PHI)

C____________________________________________________
C MATRIX VERSION OF VECTORIZED COEFFICIENTS OF SPHE HAR EXPANSION
      real*8 AUXFGA(lmax,lmax,NRGMAX),AUXFGB(lmax,lmax,NRGMAX)
     1   ,mtC(lmax,lmax,NRGMAX), mtD(lmax,lmax,NRGMAX)

      INTEGER L                 !INDEX RUNNING OVER COEFF. IN SPHE. HAR

      INTEGER N,M               !INDEXES RUNNING OVER EXPANSION SPHE. HAR. 
C                               !N-1=HARMONIC USED, M-1="THETA INDEX" USED 
      

C_______________________________________________________
C PARTICLES
      INTEGER I                 !RADIAL CELL RUNNING INDEX
      INTEGER IPart             !PARTICLE INDEX
      INTEGER II                !RADIAL CELL INDEX OF EACH PARTICLE
      REAL*8 DI                 !NORMALIZED DISTANCE FROM CELL
      INTEGER I0                !NEAREST CELL INDEX
      INTEGER I1                !NEXT NEAREST CELL INDEX

      double precision cth,sth  !cos(th) & sin(th)

      
C______________________________________________________
C ACCELERATION COMPUTATION
      REAL*8 F1                 !RADIAL FORCE
      REAL*8 F2                 !THETA FORCE
      REAL*8 F3                 !PHI FORCE

      REAL*8 FX,FY,FZ           !CARTESIAN COORDINATES FORCE

      REAL*8 SUMM_A,SUMM_B      !AUX VARIABLES FOR SUMMATION OVER
C                               !SPHERICAL HARM. EXPANSION
      REAL*8 SUMM1,SUMM2        !FIRST AND SECOND CONTRIBUTION 
C                               !FROM INTERPOLATION

      REAL*8 FFR                !RADIAL FORCE FOR PARTICLES OUTSIDE GRID
      REAL*8 GM                 !GRAV*MASS INSIDE GRID


C_______________________________________________________
C SELF FORCE
      INTEGER JJ                !cell index for self force erasure
      REAL*8 SPACE              !spacing in cell for self force 
C                               !(as DI for radial cell)
      

C______________________________________________________
C A CONSTANT:
      INTEGER NTRIG             !STOP VALUE FOR GENERATING COS(m*Phi)...

C_______________________________________________________
C SET N_TRIG TO AVOID PROBLEMS WITH L=0 CALCULATION
      NTRIG=NHAR+1
      IF(NTRIG.le.2) NTRIG=2

C_____________________________________________________
C REWRITES THE COEFF. IN A FRIENDLIER MATRIX BASED WAY

      DO I=1,IL                 !LOOP OVER RADIAL GRID
         L=0
         DO N=1,NHAR+1          !LOOP OVER SPHE HAR.
            DO M=1,N
               L=L+1
               
               AUXFGA(N,M,I)=((N-1)*FFA(L,I)-N*GGA(L,I))/FLOAT(2*N-1)
               AUXFGB(N,M,I)=((N-1)*FFB(L,I)-N*GGB(L,I))/FLOAT(2*N-1)
               mtC(N,M,I)=CNM(L,I)
               mtD(N,M,I)=DNM(L,I)

            ENDDO
         ENDDO
      ENDDO     

C________________________________________________________
C     LOOP OVER PARTICLES

      DO IPart=1,NP

         II=IP(IPart)               !RADIAL CELL ASSIGNMENT
         
         IF (II.GE.IL1) GO TO 20 !OUT OF GRID PARTICLE: ONLY RADIAL ACCEL

C______________________________________________________
c generates legendre polynomials for the particle

         do N=1,lmax
            do M=1,lmax
               plm(N,M)=0.d0
            enddo
         enddo
         CTH=Z(IPART)/S(IPART)  !COS THETA
         STH=DSQRT(DABS(1.D0-CTH*CTH)) !SIN THETA

C______________________________________________________
C  VERY EURISTIC BUT SOMEWHAT TESTED...
C  USED TO AVOID CTH=+/-1 

         if(1.-(DABS(cth)).le.eps) then
            if(cth.gt.0.) then
               cth=1.-eps
               STH=DSQRT(DABS(1.D0-CTH*CTH)) !SIN THETA
            else
               cth=-1.+eps
               STH=DSQRT(DABS(1.D0-CTH*CTH)) !SIN THETA
            endif
         else
         endif
         
C______________________________________________________
C LEGENDRE POLYNOMIALS FOR EACH PARTICLE
         call  plgndrMT(NHAR,CTH) 

C______________________________________________________
C generate the array containing cos(m phi) and sin(m phi)
         do i=1,NTRIG
            cmphi(i)=cos((float(i-1)*PH(IPart)))
            smphi(i)=sin((float(i-1)*PH(IPart)))
         enddo
         
 
C______________________________________________________
C RADIAL CELL ASSIGNMENT (IPart == PARTICLE INDEX)
         II=IP(IPart)  
         DI=(S(IPart)-RR(II))/DR(II) !NORMALIZED DISTANCE FROM CELL

C______________________________________________________
C SELF FORCE: 
         JJ=DI*NDIV
         SPACE=(DI-float(JJ)/float(NDIV))*float(NDIV)
         JJ=JJ+1

C______________________________________________________
C van Albada recipe for assignment to the two nearest cells
         IF (DI.GT.0.5) GO TO 11
         I0=II
         I1=II+1
         GO TO 12
 11      DI=1.0-DI
         I0=II+1
         I1=II
 12      CONTINUE

C_______________________________________________________
C START OF ACCELERATION COMPUTATION:
C ACCELERATION IS CALCULATED WITH LINEAR INTERPOLATION FROM TWO
C NEAREST CELLS, I0 & I1
C FOR DETAILS ABOUT THE PROCEDURE SEE TECNICAL PAPER 

C________________________________________________________
C RADIAL FORCE: FIRST COMPUTE THEN NORMALIZE
C     
         
         SUMM1=0.               !INTERPOLATION FIRST CONTRIBUTE
         
         DO M=1,NHAR+1          !LOOP OVER M (THETA INDEX FOR SPHE HAR)
            
            SUMM_A=0.
            SUMM_B=0.
            
            DO N=M,NHAR+1       !LOOP OVER SPHE HAR INDEX
               
               SUMM_A=SUMM_A+PLM(N,M)*AUXFGA(N,M,I0)
               SUMM_B=SUMM_B+PLM(N,M)*AUXFGB(N,M,I0)
               
               
            ENDDO
            
            SUMM1=SUMM1+SUMM_A*CMPHI(M)+SUMM_B*SMPHI(M)
            
 

         ENDDO
         
         SUMM2=0.               !INTERPOLATION SECOND CONTRIBUTE
         
         DO M=1,NHAR+1
            
            SUMM_A=0.
            SUMM_B=0.
            
            DO N=M,NHAR+1
               
               SUMM_A=SUMM_A+PLM(N,M)*AUXFGA(N,M,I1)
               SUMM_B=SUMM_B+PLM(N,M)*AUXFGB(N,M,I1)
               
            ENDDO
            
            SUMM2=SUMM2+SUMM_A*CMPHI(M)+SUMM_B*SMPHI(M)
            
         ENDDO
         
         
         
C_______________________________________________________                  
         F1=PI4G/S(IPart)*(SUMM1*(1.-DI)+SUMM2*DI) !RADIAL FORCE ASSIGNMENT


C________________________________________________________
C THETA FORCE: FIRST COMPUTE THEN NORMALIZE     
         
         SUMM1=0.               !INTERPOLATION FIRST CONTRIBUTE
         
         DO M=1,NHAR+1
            
            SUMM_A=0.
            SUMM_B=0.
            
            DO N=M,NHAR+1
               
               SUMM_A=SUMM_A+DPLM(N,M)*mtC(N,M,I0)
               SUMM_B=SUMM_B+DPLM(N,M)*mtD(N,M,I0)
               
            ENDDO
            
            SUMM1=SUMM1+SUMM_A*CMPHI(M)+SUMM_B*SMPHI(M)
            
         ENDDO
         
         SUMM2=0.               !INTERPOLATION SECOND CONTRIBUTE
         
         DO M=1,NHAR+1
            
            SUMM_A=0.
            SUMM_B=0.
            
            DO N=M,NHAR+1
               
               SUMM_A=SUMM_A+DPLM(N,M)*mtC(N,M,I1)
               SUMM_B=SUMM_B+DPLM(N,M)*mtD(N,M,I1)
               
            ENDDO
            
            SUMM2=SUMM2+SUMM_A*CMPHI(M)+SUMM_B*SMPHI(M)
            
         ENDDO
         
         
         
         
C____________________________________________________        
         F2=STH*(SUMM1*(1.-DI)+SUMM2*DI)/S(IPart) !THETA FORCE ASSIGNMENT

   

C________________________________________________________
C PHI FORCE: FIRST COMPUTE THEN NORMALIZE   

         
         SUMM1=0.               !INTERPOLATION FIRST CONTRIBUTE
         
         DO M=1,NHAR+1
            
            SUMM_A=0.
            SUMM_B=0.
            
            DO N=M,NHAR+1
               
               SUMM_A=SUMM_A+PLM(N,M)*mtC(N,M,I0)
               SUMM_B=SUMM_B+PLM(N,M)*mtD(N,M,I0)
               
            ENDDO
            
            SUMM1=SUMM1+float(M-1)*(-SUMM_A*SMPHI(M)+SUMM_B*CMPHI(M))
            
         ENDDO
         
         SUMM2=0.               !INTERPOLATION SECOND CONTRIBUTE
         
         DO M=1,NHAR+1
            
            SUMM_A=0.
            SUMM_B=0.
            
            DO N=M,NHAR+1
               
               SUMM_A=SUMM_A+PLM(N,M)*mtC(N,M,I1)
               SUMM_B=SUMM_B+PLM(N,M)*mtD(N,M,I1)
               
            ENDDO

            SUMM2=SUMM2+float(M-1)*(-SUMM_A*SMPHI(M)+SUMM_B*CMPHI(M))
            
         ENDDO

C_______________________________________________________                  
         F3=-(SUMM1*(1.-DI)+SUMM2*DI)/(STH*S(IPart)) !PHI FORCE ASSIGNMENT

C______________________________________________________
C END OF FORCE CALCULATION IN SPHERICAL COORDINATES


C_______________________________________________________
C SELF FORCE ERASURE !mt08/04/03

         F1=F1-(SF(II,JJ)*(1.-SPACE)+SF(II,JJ+1)*SPACE)
                

C________________________________________________________
C FROM SPHERICAL TO CARTESIAN COORDINATES

C________________________________________________________
C FORCE IN CARTESIAN COORD:         
         FX=(F1*STH+F2*CTH)*CMPHI(2) - F3*SMPHI(2)
         FY=(F1*STH+F2*CTH)*SMPHI(2) + F3*CMPHI(2)
         FZ= F1*CTH-F2*STH


C________________________________________________________
C END OF PARTICLES INSIDE GRID NOW PARTICLES OUTSIDE GRID

         GO TO 30

C_________________________________________________________
C  PARTICLES OUTSIDE GRID
C  ONLY RADIAL FORCE FROM TOTAL MASS INSIDE R PARTICLE

C COMPUTE GRAV*MASS_{INSIDE GRID}
   20 GM=GRAVC*MTOT/FLOAT(NP)*FLOAT(NP-NOUT)


      FFR=-GM/(S(IPart)**2)         !RADIAL FORCE

C__________________________________________________________
      FX= FFR*SIN(TH(IPart))*COS(PH(IPart)) !BACK TO CARTESIAN
      FY= FFR*SIN(TH(IPart))*SIN(PH(IPart))
      FZ= FFR*COS(TH(IPart))

C__________________________________________________________
   30 CONTINUE

C__________________________________________________________
C END OF FORCE CALCULATION

C__________________________________________________________
C LEAP FROG INTEGRATION: POSITION & VELOCITY UPDATING


C__________________________________________________________                  
C UPDATE VELOCITIES & POSITIONS
      VX(IPart)= VX(IPart) + FX*DT
      VY(IPart)= VY(IPart) + FY*DT
      VZ(IPart)= VZ(IPart) + FZ*DT
      X(IPart)=  X(IPart) + VX(IPart)*DT
      Y(IPart)=  Y(IPart) + VY(IPart)*DT
      Z(IPart)=  Z(IPart) + VZ(IPart)*DT
      
C__________________________________________________________     
C UPDATE FORCES
      FX0(IPart)=FX
      FY0(IPart)=FY
      FZ0(IPart)=FZ

C__________________________________________________________
C UPDATE POSITIONS IN CARTESIAN SPHERICAL COORDINATES
      S(IPart)= SQRT(X(IPart)*X(IPart)
     &+Y(IPart)*Y(IPart)+Z(IPart)*Z(IPart))


C___________________________________________
C CHECK THAT S(IPART) IS DIFFERENT FROM ZERO         MT 29/03/03 ?IS USEFUL?
      IF(S(IPART).LE.EPS1*EPS1) THEN
         Z(IPART)=0.
         Y(IPART)=0.
         X(IPART)=EPS1
      ELSE
      ENDIF


      TH(IPart)= ACOS(Z(IPart)/S(IPart))
      PH(IPart)= ATAN2(Y(IPart),X(IPart))
      IF (PH(IPart).LT.0.0) PH(IPart)=PH(IPart)+2.0*PI     

         

      ENDDO                     !END LOOP OVER PARTICLES

C____________________________________________________________
C                    UPDATE TIME
      
      TIME=TIME+DT


      RETURN
      END
C____________________________________________________________
C END OF ACCELERATION SUBRUOTINE


C_________________________________________________________
C*********************************************************
C_________________________________________________________


C_________________________________________________________
C ASSIGN PARTICLES TO RADIAL GRID
C TAKEN FROM THE ORIGINAL VAN ALBADA CODE

      SUBROUTINE P2G()

      IMPLICIT NONE
      
C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

C____________________________________________________
C COMMON FOR GRID QUANTITIES
      INCLUDE 'com_grid.blk'

C____________________________________________________
C LOCAL VARIABLES
      INTEGER L                 !INDEX RUNNING OVER PARTICLES 
      INTEGER IR                !INDEX FOR LOOP
      INTEGER II                !CELL NUMBER OF EACH PARTICLE

C____________________________________________________

      NOUT=0                    !NUMBER OF PARTICLES OUT OF GRID INIT

C____________________________________________________
C LOOP OVER PARTICLES
      DO 100 L=1,NP

      IF (S(L).LT.RR(IL)) GO TO 20

C     PARTICLES OUTSIDE GRID !ASSIGN OUT OF GRID INDEX AND GO TO NEXT PART
      NOUT=NOUT+1
      IP(L)=IL+1 
      GO TO 100

C
   20 CONTINUE
C                   PARTICLES INSIDE GRID
C                   REDISTRIBUTE PARTICLES ON GRID
      IF (IP(L).GT.IL1) IP(L)=IL1
      II=IP(L)
      IF (S(L).LT.RR(II)) GO TO 50
      IF (S(L).GE.RR(II+1)) GO TO 60
      GO TO 80
   50 CONTINUE
C                     PARTICLE IN NEW RING
      DO 51 IR=1,IL2
      II=II-1
      IF (S(L).GE.RR(II)) GO TO 80
   51 CONTINUE
C
   60 CONTINUE
      DO 61 IR=1,IL2
      II=II+1
      IF (S(L).LT.RR(II+1)) GO TO 80
   61 CONTINUE
C
   80 CONTINUE


C_________________________________________________________
C INDEX ASSIGNMENT
      IP(L)=II


 100  CONTINUE                  !END LOOP OVER PARTICLES

      RETURN
      END


C_________________________________________________________
C*********************************************************
C_________________________________________________________


C_________________________________________________________
C END OF SIMULATION SUBROUTINE: STORES PARTICLE DATA FOR RERUN
c normally not used anymore
c
      SUBROUTINE ENDSIM()
   
      IMPLICIT NONE
      
C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'


C_____________________________________________________
      INTEGER I                 !INDEX RUNNING OVER PARTICLES

C_____________________________________________________
C OPEN FILES FOR STORING
      OPEN(8,FILE='POS@EOS.DAT')
      OPEN(9,FILE='VEL@EOS.DAT')
      

C_________________________________________________________
c SAVES POSITION & VELOCITIES OF PARTICLES
      do i =1,NP
         write(8,1010) x(i),y(i),z(i)
         write(9,1010) vx(i),vy(i),vz(i)
      enddo

 1010 FORMAT(E16.8,1X,E16.8,1X,E16.8)

C_________________________________________________________
      close(8)                  !close files
      close(9)

      RETURN
      END
