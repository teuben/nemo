CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CGS - Collisionless Galactic Simulator - Source Files  C  
C                                                         C
C  Created by M. Trenti & T.van Albada in Fortran77 2003  C
C                                                         C
C  Version 1.0 Alpha                                      C
C                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C SUBRUOTINE TO INIT VARIUS QUANTITIES: TO BE CALLED ONCE
C___________________________________________________
      SUBROUTINE INITPARAMETER()  
      
      IMPLICIT NONE

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C COMMON FOR GRID QUANTITIES
      INCLUDE 'com_grid.blk'     

C SPHERICAL HARMONIC EXPANSION COEFFICIENTS
C___________________________________________________
      INCLUDE 'com_har.blk'

C_____________________________________________________
C COMMON FOR CENTER OF MASS DISPLACEMENTS
      INCLUDE 'com_cmass.blk'

C___________________________________________________

      PI=3.1415926535898d0      !PI & GRAV CONST DEFINITION                      
      GRAVC=4.4971d0                                                       
      PI4G=4.0*PI*GRAVC

C___________________________________________________
C READ PARAMETER FILE
      
      OPEN(10,FILE='PARAMETER.DAT')

      READ(10,*) IL             !RADIAL GRID NUMBER
      READ(10,*) NP             !NUMBER OF PARTICLES
      READ(10,*) LSTEP          !NUMBER OF STEPS (MAX)
      READ(10,*) INCM           !FREQUENCY OF CMSS CALL 
      READ(10,*) INPR           !FREQUENCY OF DIAGNOSTIC CALL
      READ(10,*) INSNAP         !FREQUENCY OF SNAPSHOT WRITING
      READ(10,*) DT             !TIME STEP
      READ(10,*) TIME           !START TIME OF SIMULATION
      READ(10,*) TMAX           !END TIME OF SIMULATION
      READ(10,*) MTOT           !TOTAL MASS OF GALAXY
      READ(10,*) PLUMMERFLAG    !FLAG FOR PLUMMER INIT
      READ(10,*) MDT            !MAX ALLOWED TIMESTEP
      READ(10,*) MINDT          !MIN ALLOWED TIMESTEP
      IL1=IL-1                                                          
      IL2=IL-2      

      CLOSE(10)

C__________________________________________________________
C CMASS log variables init:
      xCMstore(1)=0.d0
      xCMstore(2)=0.d0
      xCMstore(3)=0.d0
      vCMstore(1)=0.d0
      vCMstore(2)=0.d0
      vCMstore(3)=0.d0


      RETURN
      END


C___________________________________________________
C***************************************************

C THIS SUBRUOTINE PROVIDE THE INIT CONDITION FOR THE CODE:    
C IT FILLS THE ARRAY WITH PARTICLE POSITION AND VELOCITIES, ASSIGNING 
C DUMMY VALUES AT FORCE AND GRID BELONGING.

      SUBROUTINE INITPART()  

      IMPLICIT NONE

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

C____________________________________________________
C LOCAL VARIABLES
      INTEGER L                 !INDEX FOR LOOP OVER PARTICLES
      INTEGER aux               !skip variables
      real aux1,aux2            !skip variables

      REAL*8 ANVEL(3)           !angular velocities for ARBEL INIT
C____________________________________________________
C ANVEL INIT:
      DATA ANVEL/3*1./

C____________________________________________________
C INIT VIRIAL RATIO FOR ARBEL2
      REAL*8 uuu

C____________________________________________________
C BOUNDED FLAG INIT
      DO L=1,NP
         Bflag(L)=1
      ENDDO

C_____________________________________________________
C CHECK FLAG FOR INIT CONDITION TYPE
       
c___________ plummer
      IF(PlummerFLag.eq.1) THEN

         CALL PLUMMER(-195,0.592,20.) !PARTICLES INIT Plummer MODEL
         CALL ROTATE()                !ROTATION OF SYSTEM Jtot==Jz

         RETURN 
      ELSE
      ENDIF

c___________ clumpy init
      IF(PlummerFLag.eq.2) THEN

         write(*,*) "I am reading the virial ratio from init_virial.dat"
         open(98,file='init_virial.dat')
         read(98,*) uuu
         
         write(*,*) 'Virial ratio: ', uuu
         
         CALL ARBEL2(1.D0,1.D0,1.D0,0.00001D0,uuu,-668,8,0.518D0,0)
         CALL ROTATE()          !ROTATION OF SYSTEM Jtot==Jz

         RETURN 
      ELSE
      ENDIF

c__________ load external ICs 
      IF(PlummerFlag.eq.3) THEN
         OPEN(10,file='initPOS.dat')
         DO L=1,NP
            READ(10,*) X(L),Y(L),Z(L)
            S(L)=sqrt(X(L)**2+Y(L)**2+Z(L)**2)

C     DUMMY-VALUES IP AND FORCES                           
         IP(L)=1                                                           
         FX0(L)=0.0                                                        
         FY0(L)=0.0                                                        
         FZ0(L)=0.0 
         ENDDO
         CLOSE(10)

         OPEN(10,file='initVEL.dat')
         DO L=1,NP
            READ(10,*) VX(L),VY(L),VZ(L)
         ENDDO
         CLOSE(10)

         RETURN 
      ELSE
      ENDIF

c__________ 
      IF(PlummerFlag.gt.3) THEN
       WRITE(*,*) "Error In PlummerFlag, Please Check PARAMETER.DAT"
     &        , PlummerFlag
      ELSE
      ENDIF

      RETURN
      END


C___________________________________________________
C***************************************************


C  COMPUTES AUXILIARY QUANTITIES FOR THE RADIAL MESH
C  STORED IN COMMON BLOCK RADMES DECLARED IN THE FILE 
C  'com_grid.blk'
C___________________________________________________
      SUBROUTINE INITGRID()   

      IMPLICIT NONE
      
C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR GRID QUANTITIES
      INCLUDE 'com_grid.blk'

C____________________________________________________
C LOCAL VARIABLES
      
      INTEGER I,N               !INDEX FOR LOOP OVER GRID AND SPHER HAR   
      REAL Q                    !GRID RADIUS RATIO AUX VARIABLE

C____________________________________________________
 
C____________________________________________________
C GRID GENERATION !!MAY BE IMPROVED!!

C     ADOPT INTERVALS FOR RADIAL MESH                     
C     LOGARITHMIC GRID
      RR(1)=0.0                                                         
      RR(2)=0.005                                                        
      DO I=3,IL                                                      
         RR(I)=RR(I-1)*10**0.125
C SET INCREASINGLY GRID SPACING 
         if((RR(I)-RR(I-1)).le.(RR(I-1)-RR(I-2))) 
     &        RR(I-1)=.5*(RR(I)+RR(I-2))

      ENDDO


C__________________________________________________
      WRITE(*,*) 'Max grid @ ',RR(IL) !WRITE MAX GRID POINT

     
      DO  I=1,IL1               !LOOP OVER RADIAL GRID
         R(I)=(RR(I+1)+RR(I))*0.5
         DR(I)=RR(I+1)-RR(I)
         DRR(I)=DR(I)*R(I)
         RIR(I+1)=1.0/RR(I+1)
         Q=RR(I)/RR(I+1)
         QQ(1,I)=1.0
         FQ(1,I)=1.0
         GQ(1,I)=0.5*(1.0+Q)
         DO N=2,LMAX            !LOOP OVER SPHERICAL  HARMONICS 
            QQ(N,I)=Q*QQ(N-1,I)
            FQ(N,I)=2.0*Q/(1.0+Q)*FQ(N-1,I)
            GQ(N,I)=0.5*(1.0+Q)*GQ(N-1,I)
         ENDDO
      ENDDO
C
      RIR(1)=0.0                !FIRST POINT ADJUSTMENT


      RETURN 
      END

C___________________________________________________
C***************************************************


C  COMPUTES AUXILIARY QUANTITIES FOR POISSON SOLVER
C  STORED IN COMMON BLOCK DECLARED IN THE FILE 
C  'com_har.blk'
C  !! WARNING !! ALL IS SET WITH LMAX=7

C___________________________________________________
      SUBROUTINE INITPOISSON()   

      IMPLICIT NONE

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C___________________________________________________
C SPHERICAL HARMONIC EXPANSION COEFFICIENTS
      INCLUDE 'com_har.blk'

C___________________________________________________
C DATA BLOCK FROM VAN ALBADA'S CODE
C !! SET WITH NHAR=7 !!
C COEFF. FOR CALCULATION OF SPHERICAL HARMONICS EXPANSION 

      DATA C(1,1) / .7957747154594719E-01/             
      DATA C(2,1) / .2387324146378402    /                              
      DATA C(2,2) / .2387324146378402    /                              
      DATA C(3,1) / .3978873577297293    /                              
      DATA C(3,2) / .1326291192432496    /                              
      DATA C(3,3) / .3315727981081196E-01/                              
      DATA C(4,1) / .5570423008216316    /                              
      DATA C(4,2) / .9284038347027179E-01/                              
      DATA C(4,3) / .9284038347027179E-02/                              
      DATA C(4,4) / .1547339724504500E-02/                              
      DATA C(5,1) / .7161972439135305    /                              
      DATA C(5,2) / .7161972439135189E-01/                              
      DATA C(5,3) / .3978873577297398E-02/                              
      DATA C(5,4) / .2842052555212405E-03/                              
      DATA C(5,5) / .3552565694015506E-04/                              
      DATA C(6,1) / .8753521870054186    /                              
      DATA C(6,2) / .5835681246702795E-01/                              
      DATA C(6,3) / .2084171873822399E-02/                              
      DATA C(6,4) / .8684049474260097E-04/                              
      DATA C(6,5) / .4824471930144493E-05/                              
      DATA C(6,6) / .4824471930144506E-06/                              
      DATA C(7,1) / 1.034507130097303    /                              
      DATA C(7,2) / .4926224429034898E-01/                              
      DATA C(7,3) / .1231556107258697E-02/                              
      DATA C(7,4) / .3420989186829807E-04/                              
      DATA C(7,5) / .1140329728943299E-05/                              
      DATA C(7,6) / .5183316949742096E-07/                              
      DATA C(7,7) / .4319430791451703E-08/                              


      RETURN
      END


C________________________________________________________
C********************************************************

C  SUBROUTINE TO INIT THE VELOCITIES OF THE PARTICLES BY 
C  SHIFTING BACK THEY BY HALF A TIME STEP AS REQUIRED TO BE
C  CONSISTENT WITH LEAP-FROG METHOD


C________________________________________________________
      SUBROUTINE INIT_VEL()


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


      double precision eps      !parameter of tolerance near cth=+/-1
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


C______________________________________________________

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
         do i=1,nhar+1
            cmphi(i)=cos((float(i-1)*PH(IPart)))
            smphi(i)=sin((float(i-1)*PH(IPart)))
         enddo
         
 
C______________________________________________________
C RADIAL CELL ASSIGNMENT (IPart == PARTICLE INDEX)
         II=IP(IPart)  
         DI=(S(IPart)-RR(II))/DR(II) !NORMALIZED DISTANCE FROM CELL

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
C FOR LEAP FROG INTEGRATION: VELOCITIES UPDATING
C VELOCITIES ARE SHIFTED BACK TO TIME -DT/2
      VX(IPART)=VX(IPART)-0.5*DT*FX
      VY(IPART)=VY(IPART)-0.5*DT*FY
      VZ(IPART)=VZ(IPART)-0.5*DT*FZ

C__________________________________________________________     
C UPDATE FORCES
      FX0(IPart)=FX
      FY0(IPart)=FY
      FZ0(IPart)=FZ
         

      ENDDO                     !END LOOP OVER PARTICLES



      RETURN
      END
C____________________________________________________________



C____________________________________________________________
C     ***************************************************************
      SUBROUTINE ROTATE()
C     ***************************************************************
C FROM VAN ALBADA's CODE: 
C            ROTATE SYSTEM, LET Z-AXIS COINCIDE WITH
C            ANGULAR MOMENTUM VECTOR

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'


C____________________________________________________
C LOCAL VARIABLES     
      real*8 amx,amy,amz,amt,da !ANGULAR MOMENTUM & DA: RATIO XY PLANE/TOT
      INTEGER L                 !INDEX FOR PARTICLE LOOP
      REAL*8 S2                 !SQUARE OF RADIAL DISTANCE FROM CM

      REAL*8 A11,A12,A13,A21,A22,A23,A31,A32,A33 !TRANF MATRIX
      
      REAL*8 XX,YY,ZZ,VXX,VYY,VZZ

C            ROTATE SYSTEM, LET Z-AXIS COINCIDE WITH
C            ANGULAR MOMENTUM VECTOR
      AMX=0.0
      AMY=0.0
      AMZ=0.0

C________________________________________________________________
C COMPUTE ANGULAR MOMENTUM
      DO 100 L=1,NP
      S2=X(L)*X(L)+Y(L)*Y(L)+Z(L)*Z(L)
      AMX=AMX + Y(L)*VZ(L) - Z(L)*VY(L)
      AMY=AMY + Z(L)*VX(L) - X(L)*VZ(L)
      AMZ=AMZ + X(L)*VY(L) - Y(L)*VX(L)
  100 CONTINUE
      AMT=SQRT(AMX*AMX+AMY*AMY+AMZ*AMZ)
      AMX=AMX/AMT
      AMY=AMY/AMT
      AMZ=AMZ/AMT
      DA=SQRT(AMX*AMX+AMY*AMY)

C______________________________________________________________
C            TRANSFORMATION MATRIX (MATHEWS AND WALKER P.404 , GAMMA=0)
C           (11 12 13)    COS TH * COS PH   COS TH * SIN PH   -SIN TH
C       A = (21 22 23)   -SIN PH            COS PH             0
C           (31 32 33)    SIN TH * COS PH   SIN TH * SIN PH    COS TH
      A11= AMZ*AMX/DA
      A12= AMZ*AMY/DA
      A13=-DA
      A21=-AMY/DA
      A22= AMX/DA
      A23= 0
      A31= AMX
      A32= AMY
      A33= AMZ

C________________________________________________________
C APPLY TRANFORMATION
      DO 200 L=1,NP
      XX=X(L)
      YY=Y(L)
      ZZ=Z(L)
      VXX=VX(L)
      VYY=VY(L)
      VZZ=VZ(L)
      X(L)=A11*XX + A12*YY + A13*ZZ
      Y(L)=A21*XX + A22*YY + A23*ZZ
      Z(L)=A31*XX + A32*YY + A33*ZZ
      VX(L)=A11*VXX + A12*VYY + A13*VZZ
      VY(L)=A21*VXX + A22*VYY + A23*VZZ
      VZ(L)=A31*VXX + A32*VYY + A33*VZZ
      S(L)=SQRT(X(L)*X(L)+Y(L)*Y(L)+Z(L)*Z(L))
      TH(L)=ACOS(Z(L)/S(L))
      PH(L)=ATAN2(Y(L),X(L))
      IF (PH(L).LT.0.0) PH(L)=PH(L)+2.0*PI
  200 CONTINUE


C
      RETURN
      END
CCC===================================================



C_____________________________________________________
C*****************************************************

C
C ADJUST GRID ACCORDINGLY TO PARTICLE DISTRIBUTION
C ENSURING THAT EACH CELL CONTAINS THE MOST POSSIBLE 
C UNIFORM NUMBER OF PARTICLES AND THAT THE CENTER IS WELL 
C REPRESENTED 
C
      SUBROUTINE ADJGRID()

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
      REAL*8 MASSFRAD(NRGMAX)   !MASS INSIDE THE FRAC RADIUS
      REAL*8 RRAUX(NRGMAX)      !NEW AUX VALUE OF GRID POINTS
      
      INTEGER J                 !INDEX RUNNING OVER RADIAL CELLS
      INTEGER IS,IG,I           !AUX INDEXES


C____________________________________________________
C FRAC MASS RADII DEFINITION
     
      MASSFRAD(1)=MTOT/FLOAT(NP)*250.
      IF(MASSFRAD(1).GT.(MTOT/FLOAT(IL-7)*.5)) 
     &                    MASSFRAD(1)=.5*MTOT/FLOAT(IL-7)
      MASSFRAD(2)=MTOT/FLOAT(NP)*750.
      IF(MASSFRAD(2).GT.(.8*MTOT/FLOAT(IL-7))) 
     &                    MASSFRAD(2)=.8*MTOT/FLOAT(IL-7)
      
      DO J=3,IL-10
         MASSFRAD(J)=MTOT/FLOAT(IL-7)*FLOAT(J-2)
      ENDDO
C_____________________________________________________
C SEARCH NEW GRID POINTS
      IS=0
      DO J=1,IL-10
         DO I=1,IL1             !GRID ARRAY SEARCHING 
            IF(MASSGRID(I).LE.MASSFRAD(J)) IS=I !POINT FOUND
         ENDDO
         

         IF(IS.EQ.0) THEN       !INSIDE FIRST CELL
            RRAUX(J)=RR(2)*MASSFRAD(J)/MASSGRID(1) !PROBLEM IN MASS 
C                                                   !DEFINITION? 
         ELSE
C_________________________________________________________
C     INTERPOLATION !NEW VERSION: ADDED IS+1 & IG+1 MT020403
            IG=IS+1
            RRAUX(J)=RR(IS+1)+(MASSFRAD(J)-MASSGRID(IS))/
     &           (MASSGRID(IG)-MASSGRID(IS))*(RR(IG+1)-RR(IS+1))
            
         ENDIF

      ENDDO


C_________________________________________________________
C NEW GRID ASSIGNMENT
      RR(1)=0.
      DO J=1,IL-10
         IF(RRAUX(J).LE.RR(J)) THEN
            WRITE(*,*) 'ERROR!',RRAUX(J),RR(J),MASSFRAD(J),MASSFRAD(J-1)
            RR(J+1)=RR(J)*1.1
         ELSE
            RR(J+1)=RRAUX(J)
         ENDIF
C SET INCREASINGLY GRID SPACING 
         if((J.GT.1).and.((RR(J+1)-RR(J)).le.(RR(J)-RR(J-1)))) 
     &        RR(J)=0.5*(RR(J+1)+RR(J-1))

      ENDDO
      
      DO J=IL-8,IL              !END WITH LOG GRID
         RR(J) = RR(J-1)*10.**.15
C      SET INCREASINGLY GRID SPACING 
       if((RR(J)-RR(J-1)).le.(RR(J-1)-RR(J-2))) 
     &        RR(J-1)=0.5*(RR(J)+RR(J-2))

      ENDDO


C__________________________________________________________
C NEW GRID AUX VARIABLES
      DO  I=1,IL1               !LOOP OVER RADIAL GRID
         R(I)=(RR(I+1)+RR(I))*0.5
         DR(I)=RR(I+1)-RR(I)
         DRR(I)=DR(I)*R(I)
         RIR(I+1)=1.0/RR(I+1)
         Q=RR(I)/RR(I+1)
         QQ(1,I)=1.0
         FQ(1,I)=1.0
         GQ(1,I)=0.5*(1.0+Q)
         DO N=2,LMAX            !LOOP OVER SPHERICAL  HARMONICS 
            QQ(N,I)=Q*QQ(N-1,I)
            FQ(N,I)=2.0*Q/(1.0+Q)*FQ(N-1,I)
            GQ(N,I)=0.5*(1.0+Q)*GQ(N-1,I)
         ENDDO

C$$$C writes new grid parameters
C$$$         write(99,*) I,RR(I),R(I),DR(I)



      ENDDO
C
      RIR(1)=0.0                !FIRST POINT ADJUSTMENT


      RETURN
      END


C_____________________________________________________
C*****************************************************

C
C GLOBAL ROUTINE FOR ADAPTIVE GRID
C
      SUBROUTINE ADAPTGRID()

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
      INTEGER I                 !LOOP INDEX


C____________________________________________________
      DO I=1,1

         CALL P2G()             !PARTICEL TO CELL
      
         CALL DENS()            !DENSITY

         CALL ADJGRID()         !ADJUST GRID
      
      ENDDO

      RETURN
      END
