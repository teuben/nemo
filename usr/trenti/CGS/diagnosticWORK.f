CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CGS - Collisionless Galactic Simulator - Source Files  C  
C                                                         C
C  Created by M. Trenti & T.van Albada in Fortran77 2003  C
C                                                         C
C  Version 1.0 Alpha                                      C
C                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C__________________________________________________________
C DIAGNOSTIC SUBROUTINES LIBRARY


C
C MODIFIED TO ALLOW THE COMPUTATION OF THE DISTRIBUTION FUNCTION!
C

      SUBROUTINE DIAGNOSTIC(ESflag)

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
C COMMON FOR ENERGY DISTRIBUTION
      INCLUDE 'com_DF.blk'

C_____________________________________________________
C COMMON FOR CENTER OF MASS DISPLACEMENTS
      INCLUDE 'com_cmass.blk'

C____________________________________________________
C LOCAL VARIABLES

      integer ESflag

      double precision eps      !paramter of tolerance near cth=+/-1
      parameter(eps=1.d-7)

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
      INTEGER IPART             !PARTICLE INDEX
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
      
C_____________________________________________________
C     ENERGY & ANGULAR MOMENTUM: GLOBAL
      
      REAL*8 EPOT               !TOTAL POTENTIAL ENERGY
      REAL*8 EKIN               !TOTAL KINETIC ENERGY
      REAL*8 ENEG               !TOTAL POTENTIAL ENERGY OF BINDED PARTICLES
      REAL*8 ETOT               !TOTAL POTENTIAL ENERGY
      REAL*8 AMX,AMY,AMZ        !TOTAL ANGULAR MOMENTUM COMPONENTS
      REAL*8 ATX,ATY,ATZ        !TOTAL TORQUE VECTOR COMPONENTS
      REAL*8 VIRIAL             !2K/W 2*EKIN/EPOT

C_____________________________________________________
C INERTIA ELLIPSOID
      REAL*8 AI(3,3)            !CARTESIAN INERTIA ELLIPSOID
      REAL*4 AW(3,3)            !AUXILIARY WORK INERTIA ELLIPSOID
C                               !FOR DIAGONALIZATION SUBROUTINE
      REAL*4 DA(3),DW(3)        !EIGENVALUES OF AI AND WORK VERSION 
C                               !OF EIGENVALUES
      INTEGER IERR              !UNUSED FLAG PASSED TO DIAG. SUBROUTINE
      REAL*4 WORK1(50)          !WORKING ARRAY OF THE DIAG. SUBROUTINE ??
      INTEGER J                 !RUNNING INDEX
      
      INTEGER MM(3)             !FLAG TO SORT EIGENVALUES
      REAL*8 P,Q                !AUX VARIABLES TO WRITE AXIS RATIO
      REAL*8 EPS1,EPS2          !AXIS RATIO OF THE SYSTEM (B/A &C/A)
     

C_____________________________________________________
C     ENERGY & ANGULAR MOMENTUM: SINGLE PARTICLE
      
      REAL*8 EP                 !SINGLE PARTICLE POTENTIAL ENERGY
      REAL*8 EK                 !SINGLE PARTICLE KINETIC ENERGY
      REAL*8 EB                 !SINGLE PARTICLE BINDING ENERGY
      REAL*8 AMSQ               !|JTOT|_{SINGLE PARTICLE}
      REAL*8 SHAPE              !ORBITAL SHAPE OF EACH PARTICLE
      REAL*8 MEANSHAPE          !MEAN ORBITAL SHAPE

C_____________________________________________________
C    SINGLE PARTICLE:

      REAL*8 VX0,VY0,VZ0        !VELOCITIES AT T+DT AS POSITIONS      
      REAL*8 AMX_SP,AMY_SP,AMZ_SP !SINGLE PARTICLE ANGULAR MOMENTUM COMPONENTS

C______________________________________________________
C A CONSTANT:
      INTEGER NTRIG             !STOP VALUE FOR GENERATING COS(m*Phi)...


C____________________________________________________
C INERTIA ELLIPSOID 
      REAL*8 COST1              !CONSTANT FOR INERTIA ELLIPSOID TRUNCATION
      PARAMETER(COST1=3.)

C__________________
      REAL*8 fraux
      integer ifr

C____________________________________________________
C CALL VARIOUS OTHER DIAGNOSTIC SUBROUTINES
      
      CALL DENS()
      CALL DENSbound()
      CALL VELDIS()
      CALL VELDISbound()
      CALL FRMRAD()

c______________________
      if(esflag.eq.10) then
         write(90,110) NP,time
      else
      endif
 110  format(I12,1X,E14.7)

C_______________________________________________________
C SET N_TRIG TO AVOID PROBLEMS WITH L=0 CALCULATION
      NTRIG=NHAR+1
      IF(NTRIG.le.2) NTRIG=2
    
C_____________________________________________________
C INIT OF ALL QUANTITIES:
      EPOT=0.0
      EKIN=0.0
      ENEG=0.0
      ETOT=0.0
      AMX=0.0
      AMY=0.0
      AMZ=0.0
      ATX=0.0
      ATY=0.0
      ATZ=0.0
      meanshape=0.

      DO I=1,3
         DO J=1,3
            AI(I,J)=0.0
         ENDDO
      ENDDO

C_____________________________________________________
C BINDING ENERGY ARRAY & LABEL INIT
      DO I=1,60
         IEBIN(I)=0
         TEKST(I)=EMIN+FLOAT(I-1)*(EMAX-EMIN)/60.0
      ENDDO
      TEKST(61)=EMAX


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

      DO IPART=1,NP

         II=IP(IPART)               !RADIAL CELL ASSIGNMENT
         
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
         F1=PI4G/S(IPART)*(SUMM1*(1.-DI)+SUMM2*DI) !RADIAL FORCE ASSIGNMENT

C_______________________________________________________
C SELF FORCE ERASURE !mt08/04/03
         F1=F1-(SF(II,JJ)*(1.-SPACE)+SF(II,JJ+1)*SPACE)


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
         F2=STH*(SUMM1*(1.-DI)+SUMM2*DI)/S(IPART) !THETA FORCE ASSIGNMENT

   

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
         F3=-(SUMM1*(1.-DI)+SUMM2*DI)/(STH*S(IPART)) !PHI FORCE ASSIGNMENT

C______________________________________________________
C END OF FORCE CALCULATION IN SPHERICAL COORDINATES


C________________________________________________________
C FROM SPHERICAL TO CARTESIAN COORDINATES

C________________________________________________________
C FORCE IN CARTESIAN COORD:         
         FX=(F1*STH+F2*CTH)*CMPHI(2) - F3*SMPHI(2)
         FY=(F1*STH+F2*CTH)*SMPHI(2) + F3*CMPHI(2)
         FZ= F1*CTH-F2*STH

C__________________________________________________________
C END OF FORCE CALCULATION


C__________________________________________________________
C  CALCULATE POTENTIAL ENERGY: 
C  SAME SCHEME AS FORCE CALCULATIONL INTERPOLATION IN TWO 
C  RADIAL GRID POINTS 
         
C__________________________________________________________
C INSIDE GRID

         SUMM1=0.               !INTERPOLATION FIRST CONTRIBUTE
         
         DO M=1,NHAR+1
            
            SUMM_A=0.
            SUMM_B=0.
            
            DO N=M,NHAR+1
               
               SUMM_A=SUMM_A+PLM(N,M)*mtC(N,M,I0)
               SUMM_B=SUMM_B+PLM(N,M)*mtD(N,M,I0)
               
            ENDDO
            
            SUMM1=SUMM1+SUMM_A*CMPHI(M)+SUMM_B*SMPHI(M)
            
         ENDDO
         
         SUMM2=0.               !INTERPOLATION SECOND CONTRIBUTE
         
         DO M=1,NHAR+1
            
            SUMM_A=0.
            SUMM_B=0.
            
            DO N=M,NHAR+1
               
               SUMM_A=SUMM_A+PLM(N,M)*mtC(N,M,I1)
               SUMM_B=SUMM_B+PLM(N,M)*mtD(N,M,I1)
               
            ENDDO
            
            SUMM2=SUMM2+SUMM_A*CMPHI(M)+SUMM_B*SMPHI(M)
            
         ENDDO
         
         
         EP=(SUMM1*(1.-DI)+SUMM2*DI) !POTENTIAL ENERGY ASSIGNMENT

C_______________________________________________________
C SELF POTENTIAL ENERGY ERASURE !mt08/04/03
         EP=EP-(SEP(II,JJ)*(1.-SPACE)+SEP(II,JJ+1)*SPACE)


         
C________________________________________________________
C END OF PARTICLES INSIDE GRID NOW PARTICLES OUTSIDE GRID

         GO TO 30

C_________________________________________________________
C  PARTICLES OUTSIDE GRID
C  ONLY RADIAL FORCE FROM TOTAL MASS INSIDE R PARTICLE

C COMPUTE GRAV*MASS_{INSIDE GRID}
 20      GM=GRAVC*MTOT/FLOAT(NP)*FLOAT(NP-NOUT)

         EP=-GM/S(IPART)        !POTENTIAL ENERGY
         FFR=-GM/(S(IPART)**2)  !RADIAL FORCE
C__________________________________________________________
         FX= FFR*SIN(TH(IPART))*COS(PH(IPART)) !BACK TO CARTESIAN
         FY= FFR*SIN(TH(IPART))*SIN(PH(IPART))
         FZ= FFR*COS(TH(IPART))

C__________________________________________________________
 30      CONTINUE
C__________________________________________________________
C END OF FORCE & POTENTIAL CALCULATION


C__________________________________________________________                  
C UPDATE VELOCITIES 
      VX0=VX(IPART)+0.5*DT*FX
      VY0=VY(IPART)+0.5*DT*FY
      VZ0=VZ(IPART)+0.5*DT*FZ


C___________________________________________________________
C CALCULATE KINETIC ENERGY AND BINDING ENERGY
      EK=VX0*VX0+VY0*VY0+VZ0*VZ0
      EB=0.5*EK+EP

C___________________________________________________________
C write down particle snapshot, MT addition for NEMO compatibility
      if(ESFLAG.eq.10) then
         
         write(90,111) IPART,X(IPART)+xCMstore(1),Y(IPART)+xCMstore(2)
     &        ,Z(IPART)+xCMstore(3),VX0+vCMstore(1),VY0+vCMstore(2)
     &        ,VZ0+vCMstore(3),FX,FY,FZ,EP
      else
      endif
 111  format(I11,10(1X,E13.5))

C______________________________________________
C SET A FLAG FOR BOUNDED PARTICLES

      if(EB.gt.0) then 
         Bflag(IPART)=0
      else
         Bflag(IPART)=1
      endif


C___________________________________________________________
C BINDING ENERGY ARRAY CONSTRUCTION:
      IND=INT((EB-EMIN)/(EMAX-EMIN)*60.0) + 1
      IF (IND.LT.1) IND=1
      IF (IND.GT.60) IND=60
      IEBIN(IND)=IEBIN(IND)+1
      
C___________________________________________________________
c ANGULAR MOMENTUM,SINGLE PARTICLE

      amx_sp= Y(IPART)*VZ0 -Z(IPART)*VY0
      amy_sp= Z(IPART)*VX0 -X(IPART)*VZ0
      amz_sp= X(IPART)*VY0 -Y(IPART)*VX0


C___________________________________________________________
c shape parameter single particle:
c shape =e*e-1=2*EB*Jsp**2/(Psi**2*r**2)
      amsq=(amx_sp**2.+amy_sp**2.+amz_sp**2.)
      if(abs(EP)*S(IPART).gt.eps) then   !used to avoid INF
         shape=2.*EB*amsq/(EP**2.*S(IPART)**2.)
      else
      endif

C___________________________________________________________
      shape=sqrt(abs(shape+1.))      !SHAPE SINGLE PARTICLE VALUE

      

C___________________________________________________________ 
      meanshape=meanshape+shape !mean shape parameter

C____________________________________________________________
C GLOBAL QUANTITIES UPDATING 

C___________________________________________________________
C ENERGY
      IF (EB.GT.0.0) GO TO 91
      ENEG=ENEG+EK+EP
   91 EPOT=EPOT+EP
      EKIN=EKIN+EK
      
C___________________________________________________________
C                   ANGULAR MOMENTUM VECTOR
      AMX =AMX + amx_sp
      AMY =AMY + amy_sp
      AMZ =AMZ + amz_sp

C___________________________________________________________
C                   TORQUE VECTOR = D/DT (ANG.MOM)  
      ATX =ATX + Y(IPART)*FZ - Z(IPART)*FY
      ATY =ATY + Z(IPART)*FX - X(IPART)*FZ
      ATZ =ATZ + X(IPART)*FY - Y(IPART)*FX

C___________________________________________________________
C                   MOMENT OF INERTIA TENSOR  FOR R<RHM*COST

      IF(S(IPART).LE.(COST1*rhm)) THEN  
C                             !rhm radius of half mass, even. * const.
         AI(1,1)=AI(1,1) + Y(IPART)**2 + Z(IPART)**2
         AI(2,2)=AI(2,2) + X(IPART)**2 + Z(IPART)**2
         AI(3,3)=AI(3,3) + X(IPART)**2 + Y(IPART)**2
         AI(1,2)=AI(1,2) - X(IPART)*Y(IPART)
         AI(1,3)=AI(1,3) - X(IPART)*Z(IPART)
         AI(2,3)=AI(2,3) - Y(IPART)*Z(IPART) 
      ELSE
      ENDIF

C___________________________________________________________
      ENDDO                     !END OF LOOP OVER PARTICLES
C___________________________________________________________

C___________________________________________________________
C ENERGY RELATED QUANTITIES
      EKIN=0.5*MTOT/FLOAT(NP)*EKIN !KINETIC ENERGY
      EPOT=0.5*MTOT/FLOAT(NP)*EPOT !POTENTIAL ENERGY
      ENEG=0.5*MTOT/FLOAT(NP)*ENEG !TOTAL ENERGY OF BOUNDED PARTICLES
      ETOT=EKIN+EPOT            !TOTAL ENERGY
      VIRIAL=-2.0*EKIN/EPOT     !VIRIAL COEFFICIENT: 2K/W
      MEANSHAPE=MEANSHAPE/FLOAT(NP) !MEAN SHAPE NORMALIZATION

C__________________________________________________________
      AMX=AMX*MTOT/FLOAT(NP)    !angular momentum normalization
      AMY=AMY*MTOT/FLOAT(NP) 
      AMZ=AMZ*MTOT/FLOAT(NP)

      ATX=ATX*MTOT/FLOAT(NP)    !TORQUE normalization
      ATY=ATY*MTOT/FLOAT(NP) 
      ATZ=ATZ*MTOT/FLOAT(NP)

C__________________________________________________________
C   MOMENT OF INERTIA: SYMETRICAL MATRIX
      AI(2,1)=AI(1,2)
      AI(3,1)=AI(1,3)
      AI(3,2)=AI(2,3)

      DO M=1,3                  !MASS NORMALIZATION
         DO N=1,3
            AI(M,N)=MTOT/FLOAT(NP)*AI(M,N)
         ENDDO
      ENDDO


C___________________________________________________________
C MOMENT OF INERTIA: SEARCH OF EIGENVALUES
      DO M=1,3
         DO N=1,3
            AW(M,N)=AI(M,N)     !AW: WORK AUXILIARY ARRAY
         ENDDO
      ENDDO

      IERR=0                    !UNUSED FLAG IN THIS VERSION OF F02AAF

C___________________________________________________________
      CALL F02AAF(AW,3,3,DW,WORK1,IERR) ! EIGENVALUES SEARCH SUBROUTINE

C___________________________________________________________
      DO  M=1,3
         DA(M)=DW(M)            !EIGENVALUES
      enddo

C___________________________________________________________
C            SORT EIGENVALUES MOMENT OF INERTIA             
C            M(1)  IS INDEX SMALLEST MOMENT OF INERTIA        
C            I.E. M(1) CORRESPONDS TO LONGEST AXIS, ETC 
      MM(1)=1                                                          
      MM(2)=2                                                        
      MM(3)=3
      IF (DA(2).LT.DA(1)) MM(1)=2    
      IF (DA(3).LT.DA(MM(1))) MM(1)=3 
      IF (DA(2).GT.DA(3)) MM(3)=2    
      IF (DA(1).GT.DA(MM(3))) MM(3)=1 
      MM(2)=6-MM(1)-MM(3)

C___________________________________________________________              
C            AXIAL RATIOS           
C            EPS1 = B/A       EPS2 = C/A 
      P=DA(MM(1))/DA(MM(2))                
      Q=DA(MM(1))/DA(MM(3))                
      EPS1=(-P*Q-P+Q)/(P*Q-P-Q)          
      EPS2=(-P*Q+P-Q)/(P*Q-P-Q)          
      IF (EPS1.LT.0.0) EPS1=0.0          
      IF (EPS2.LT.0.0) EPS2=0.0          
      EPS1=SQRT(EPS1)           !NORMALIZATION                    
      EPS2=SQRT(EPS2)




C____________________________________________________________
C WRITE A FILE WITH DIAGNOSTIC, STORED IN FORT.*:
      
C MAIN STORING
      WRITE(20,1095) TIME,ETOT,VIRIAL,SQRT(AMX**2+AMY**2+AMZ**2)
      WRITE(19,1095) TIME,AMX,AMY,AMZ,meanshape

C AXIS RATIO & EIGENVALUES OF INERTIA ELLIPSOID 
      write(16,1097) time,DA(MM(1)),DA(MM(2)),DA(MM(3))
      write(17,1096) time,EPS1,EPS2

 1097 FORMAT(F12.6,1X,F12.6,1X,F12.6,1X,F12.6)
 1096 FORMAT(F12.6,1X,F12.6,1X,F12.6)
 1095 FORMAT(F11.6,1X,F13.8,1X,F13.8,1X,F13.8,1X,F13.8)    

      call flush(16)
      call flush(17) 
      CALL FLUSH(20)
      CALL FLUSH(19)
      
c___________________________________________________________
c old Van Albada's fort.2 file
      WRITE(2,2000)
      WRITE(2,2010) TIME,EKIN,EPOT,ETOT,ENEG,VIRIAL
      WRITE(2,2000)
      WRITE(2,2016) AMX, (AI(1,M),M=1,3), ATX
      WRITE(2,2017) AMY, (AI(2,M),M=1,3), ATY
      WRITE(2,2016) AMZ, (AI(3,M),M=1,3), ATZ
      WRITE(2,2011)
      WRITE(2,2012) (TEKST(I),I=1,21)
      WRITE(2,2013) (IEBIN(I),I=1,20)
      WRITE(2,2012) (TEKST(I),I=21,41)
      WRITE(2,2013) (IEBIN(I),I=21,40)
      WRITE(2,2012) (TEKST(I),I=41,61)
      WRITE(2,2013) (IEBIN(I),I=41,60)
      WRITE(2,2000)
      WRITE(2,2015) (R(I),I=1,IL1)
      WRITE(2,2014) (NPAR(I),I=1,IL1)
      WRITE(2,*) 'Nout = ',NOUT

      WRITE(2,2000)
      

 2000 FORMAT(1H0)
 2010 FORMAT(1H0,'TIME=',F6.3,5X,'EKIN=',F10.4,3X,'EPOT=',F10.4,3X,
     1 'ETOT=',F10.4,3X,'ENEG=',F10.4,3X,3X,'2T/W=',F9.4)
 2011 FORMAT(1H0,'DISTRIBUTION OF BINDING ENERGIES')
 2012 FORMAT(1H0,21F6.1)
 2013 FORMAT(1H ,2X,20I6)
 2016 FORMAT(1H ,4X,F13.6,12X,3(F13.6,2X),11X,F13.6)
 2017 FORMAT(1H ,'AM= ',F13.6,3X,'INERTIA= ',3(F13.6,2X),3X,
     1 'TORQUE= ',F13.6)
 2014 FORMAT(1H ,'NPAR',6X,15I8)
 2015 FORMAT(1H ,'R  ',3X,15F8.4)

      RETURN
      END


C*********************************************************
C_________________________________________________________
C*********************************************************

C_________________________________________________________
C TAKEN FROM THE ORIGINAL VAN ALBADA CODE
C COMPUTES DENSITY AND MASS INSIDE RADIAL SHELL

      SUBROUTINE DENS()

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
C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C LOCAL VARIABLES
      INTEGER L                 !INDEX RUNNING OVER PARTICLES 
      INTEGER IR                !INDEX FOR LOOP
      INTEGER II                !CELL NUMBER OF EACH PARTICLE
      INTEGER I                 !INDEX FOR DENSITY CELL

      REAL*8 Q                  !FRACTION OF PARTICLE TO BE 
C                               !ASSIGNED TO MAIN CELL
      DOUBLE PRECISION VOL      !VOLUME OF RADIAL CELL

      REAL*8 RHOzero(NRGMAX)    !density calculated with zero order kernel

C____________________________________________________

      NOUT=0                    !NUMBER OF PART OUTSIDE GRID INIT


C____________________________________________________
C      RHO & NUMBER OF PARTICLE INIT 
      DO I=1,IL
         RHO(I)=0.
         NPAR(I)=0
      ENDDO

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

C MASS ASSIGNMENT
      Q=(S(L)-RR(II))/DR(II)
      
      IF (Q-0.5) 90,91,92
 90   IF (II.EQ.1) GO TO 91
      RHO(II)=RHO(II)+0.5+Q
      RHO(II-1)=RHO(II-1)+0.5-Q
      GO TO 99
   91 RHO(II)=RHO(II)+1.0
      GO TO 99
   92 IF (II.EQ.IL1) GO TO 91

      RHO(II)=RHO(II)+1.5-Q
      RHO(II+1)=RHO(II+1)-0.5+Q
   99 CONTINUE

C_________________________________________________________
C INDEX ASSIGNMENT
      IP(L)=II
      NPAR(II)=NPAR(II)+1


 100  CONTINUE                  !END LOOP OVER PARTICLES

C_________________________________________________________
C MASS GRID LOOP

      MASSGRID(1)=RHO(1)        !PARTICLES INSIDE FIRST CELL
      DO I=2,IL1                !ADD PARTICLES
         MASSGRID(I)=MASSGRID(I-1)+RHO(I)
      ENDDO
      
C_________________________________________________________
C NORMALIZATION
      
      DO I=1,IL1
         MASSGRID(I)=MASSGRID(I)*MTOT/FLOAT(NP) !mass normalization
         VOL= 4.d0/3.d0*PI*(RR(I+1)**3-RR(I)**3) !shell volume
         RHO(I)= MTOT/FLOAT(NP)*RHO(I)/VOL !normalization
         RHOzero(I)=NPAR(I)*MTOT/FLOAT(NP)/VOL
      ENDDO



C__________________________________________________________
C WRITES OUTPUT
      write(18,*) time, ' TIME R RHO M(R) RHOzero NPAR '

      DO I=1,IL1
         write(18,1033) R(I),RHO(I),MASSGRID(I), RHOzero(I),NPAR(I)
      ENDDO

 1033 FORMAT(F11.6,1X,F12.6,1X,F10.6,1X,F12.6,1X,I8)

      RETURN
      END
C_________________________________________________________


C_________________________________________________________
C TAKEN FROM THE ORIGINAL VAN ALBADA CODE
C COMPUTES DENSITY AND MASS INSIDE RADIAL SHELL

      SUBROUTINE DENSbound()

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
C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C LOCAL VARIABLES
      INTEGER L                 !INDEX RUNNING OVER PARTICLES 
      INTEGER IR                !INDEX FOR LOOP
      INTEGER II                !CELL NUMBER OF EACH PARTICLE
      INTEGER I                 !INDEX FOR DENSITY CELL

      REAL*8 Q                  !FRACTION OF PARTICLE TO BE 
C                               !ASSIGNED TO MAIN CELL
      DOUBLE PRECISION VOL      !VOLUME OF RADIAL CELL

      REAL*8 RHOzero(NRGMAX)    !density calculated with zero order kernel

C____________________________________________________

      NOUT=0                    !NUMBER OF PART OUTSIDE GRID INIT


C____________________________________________________
C      RHO & NUMBER OF PARTICLE INIT 
      DO I=1,IL
         RHO(I)=0.
         NPAR(I)=0
      ENDDO

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


C________________________________________
C UNBOUNDED PARTICLES: SKIP ASSIGNMENT
      IF(Bflag(L).ne.1) goto 100


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

C MASS ASSIGNMENT
      Q=(S(L)-RR(II))/DR(II)
      
      IF (Q-0.5) 90,91,92
 90   IF (II.EQ.1) GO TO 91
      RHO(II)=RHO(II)+0.5+Q
      RHO(II-1)=RHO(II-1)+0.5-Q
      GO TO 99
   91 RHO(II)=RHO(II)+1.0
      GO TO 99
   92 IF (II.EQ.IL1) GO TO 91

      RHO(II)=RHO(II)+1.5-Q
      RHO(II+1)=RHO(II+1)-0.5+Q
   99 CONTINUE

C_________________________________________________________
C INDEX ASSIGNMENT
      IP(L)=II
      NPAR(II)=NPAR(II)+1


 100  CONTINUE                  !END LOOP OVER PARTICLES

C_________________________________________________________
C MASS GRID LOOP

      MASSGRID(1)=RHO(1)        !PARTICLES INSIDE FIRST CELL
      DO I=2,IL1                !ADD PARTICLES
         MASSGRID(I)=MASSGRID(I-1)+RHO(I)
      ENDDO
      
C_________________________________________________________
C NORMALIZATION
      
      DO I=1,IL1
         MASSGRID(I)=MASSGRID(I)*MTOT/FLOAT(NP) !mass normalization
         VOL= 4.d0/3.d0*PI*(RR(I+1)**3-RR(I)**3) !shell volume
         RHO(I)= MTOT/FLOAT(NP)*RHO(I)/VOL !normalization
         RHOzero(I)=NPAR(I)*MTOT/FLOAT(NP)/VOL
      ENDDO



C__________________________________________________________
C WRITES OUTPUT
      write(28,*) time, ' TIME R RHO M(R) RHOzero NPAR '

      DO I=1,IL1
         write(28,1032) R(I),RHO(I),MASSGRID(I), RHOzero(I),NPAR(I)
      ENDDO

 1032 FORMAT(F11.6,1X,F12.6,1X,F10.6,1X,F12.6,1X,I8)
      
      RETURN
      END
C_________________________________________________________



C_________________________________________________________
C*********************************************************
C_________________________________________________________

C
C COMPUTE FRACTIONAL MASS RADII 
C

C_________________________________________________________
      SUBROUTINE FRMRAD()

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
C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C LOCAL VARIABLES
      REAL*8 MASSFRAD(14)       !MASS INSIDE THE FRAC RADIUS
      INTEGER I,J,IS,IG         !VARIOUS INDEXES
      
C____________________________________________________
C FRAC MASS RADII DEFINITION
      MASSFRAD(1)=MTOT/FLOAT(NP)*50.
      IF(MASSFRAD(1).GT.1.E-3) MASSFRAD(1)=1.E-3
      MASSFRAD(2)=MTOT/FLOAT(NP)*100.
      IF(MASSFRAD(2).GT.2.E-3) MASSFRAD(2)=2.E-3
      MASSFRAD(3)=5.E-3*MTOT
      MASSFRAD(4)=1.E-2*MTOT
      MASSFRAD(5)=2.E-2*MTOT
      MASSFRAD(6)=5.E-2*MTOT
      MASSFRAD(7)=1.E-1*MTOT
      MASSFRAD(8)=3.E-1*MTOT
      MASSFRAD(9)=5.E-1*MTOT
      MASSFRAD(10)=7.E-1*MTOT
      MASSFRAD(11)=9.E-1*MTOT
      MASSFRAD(12)=9.5E-1*MTOT
      MASSFRAD(13)=9.8E-1*MTOT
      MASSFRAD(14)=9.9E-1*MTOT
      
      
C________________________________________________________
C LOOP OVER FRAC MASS RADII



      IS=0
      DO J=1,14
         DO I=1,IL1              !GRID ARRAY SEARCHING 
            IF(MASSGRID(I).LE.MASSFRAD(J)) IS=I !POINT FOUND
         ENDDO
         
         IF(IS.EQ.0) THEN       !INSIDE FIRST CELL
            RADFM(J)=RR(2)*MASSFRAD(J)/MASSGRID(1) !PROBLEM IN MASS 
C                                                   !DEFINITION? 
         ELSE
C_________________________________________________________
C     INTERPOLATION
            IG=IS+1
            RADFM(J)=RR(IS)+(MASSFRAD(J)-MASSGRID(IS))/
     &           (MASSGRID(IG)-MASSGRID(IS))*(RR(IG)-RR(IS))
            
         ENDIF
         
      ENDDO


C________________________________________________________
      rhm=radfm(9)


      write(*,*) 'In radfm rhm = ', rhm ,radfm(i)

C_________________________________________________________
C WRITE OUTPUT IN FORT.11 & FORT.12
      WRITE(11,2012) TIME,RADFM(1),RADFM(2),RADFM(4),RADFM(6)
      WRITE(12,2012) TIME,RADFM(9),RADFM(10),RADFM(12),RADFM(14)

C_________________________________________________________
C WRITE OUTPUT IN FORT.2
      WRITE(2,2010) (MASSFRAD(I),I=1,14)
      WRITE(2,2011) (RADFM(I),I=1,14) 

   
 2010 FORMAT(1H0,'FRAC MASS',15F8.4)
 2011 FORMAT(1H ,'RADIUS',3X,15F8.4)
 2012 FORMAT(F8.4,1X,F11.6,1X,F11.6,1X,F11.6,1X,F11.6) 
      RETURN
      END



C_________________________________________________________
C*********************************************************
C_________________________________________________________

C
C COMPUTE VELOCITY DISPERSION IN RADIAL SHELL
C ALSO COMPUTE ANISOTROPY PARAMETER AND KR/KT
C 

C_____________________________________________________
      SUBROUTINE VELDIS()

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
C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C LOCAL VARIABLES
      REAL*8 V2,VR,VR2,VT2      !VARIOUS VELOCITIES
      
      REAL*8 SR(NRGMAX)         !MEAN RADIAL VELOCITY DISPERSION
      REAL*8 ST(NRGMAX)         !MEAN TANGENTIAL VELOCITY DISPERSION
      REAL*8 ANIS(NRGMAX)       !ANISOTROPY PARAMETER

      INTEGER L                 !INDEX FOR LOOP OVER PARTICLES
      INTEGER II                 !INDEX FOR LOOP OVER SHELL

      REAL*8 KKR,KKT            !RADIAL & TANGENTIAL KINETIC ENERGY

C____________________________________________________
C DISPERSION INIT
      DO II=1,IL
         SR(II)=0.
         ST(II)=0.
         NPAR(II)=0.
      ENDDO

C____________________________________________________
C LOOP OVER PARTICLE INSIDE GRID
      DO L=1,NP
         II=IP(L)
           
         IF(II.LE.IL) THEN      !INSIDE GRID
            NPAR(II)=NPAR(II)+1
            V2=VX(L)*VX(L)+VY(L)*VY(L)+VZ(L)*VZ(L) !V*V
            VR=(VX(L)*COS(PH(L))+VY(L)*SIN(PH(L)))*SIN(TH(L))+VZ(L) !RADIALVEL 
     &*COS(TH(L)) 
            VR2=VR*VR           
            VT2=V2-VR2          !TANGENTIAL VEL ^2

            SR(II)=SR(II)+VR2
            ST(II)=ST(II)+VT2
         ELSE
         ENDIF
         
      ENDDO                     !END LOOP OVER PART


C____________________________________________________
C COMPUTE KINETIC TANGENTIAL ENERGY  & KINETIC RADIAL ENERGY

      KKR=0.                    !INIT 
      KKT=0.
      DO II=1,IL                !LOOP OVER CELLS

         IF(NPAR(II).GT.0) THEN 
   
            KKR=KKR+SR(II)      !RADIAL & TANGENTIAL KINETIC ENERGY
            KKT=KKT+ST(II)
   
            SR(II)=SR(II)/FLOAT(NPAR(II)) !NORMALIZATION
            ST(II)=ST(II)/FLOAT(NPAR(II)) 

            IF (SR(II).LT.1.0E-10) SR(II)=1.0E-10 !USED TO AVOID INF
            
            ANIS(II)=ST(II)/SR(II) !ANISOTROPY PARAMETER
             
         ELSE
            ANIS(II)=0.
         ENDIF


      ENDDO



C____________________________________________________
C     NORMALIZATION LOOP
         KKR=.5*KKR*MTOT/FLOAT(NP)
         KKT=.5*KKT*MTOT/FLOAT(NP)
         
      
  

C_____________________________________________________
C WRITE OUTPUT
      write(14,1099) time,KKR,KKT,2.*KKR/KKT
      
      write(13,*) time, ' TIME ___R ___VR ___VT ___ANIS ___'
      DO II=1,IL 
         WRITE(13,1099) R(II),SR(II),ST(II),ANIS(II)
      ENDDO

 1099 FORMAT(F12.6,1X,F12.6,1X,F12.6,1X,F12.6)

      RETURN
      END



C_____________________________________________________________
C*************************************************************


C
C COMPUTE VELOCITY DISPERSION IN RADIAL SHELL
C ALSO COMPUTE ANISOTROPY PARAMETER AND KR/KT
C 

C_____________________________________________________
      SUBROUTINE VELDISbound()

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
C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C LOCAL VARIABLES
      REAL*8 V2,VR,VR2,VT2      !VARIOUS VELOCITIES
      
      REAL*8 SR(NRGMAX)         !MEAN RADIAL VELOCITY DISPERSION
      REAL*8 ST(NRGMAX)         !MEAN TANGENTIAL VELOCITY DISPERSION
      REAL*8 ANIS(NRGMAX)       !ANISOTROPY PARAMETER

      INTEGER L                 !INDEX FOR LOOP OVER PARTICLES
      INTEGER II                 !INDEX FOR LOOP OVER SHELL

      REAL*8 KKR,KKT            !RADIAL & TANGENTIAL KINETIC ENERGY

C____________________________________________________
C DISPERSION INIT
      DO II=1,IL
         SR(II)=0.
         ST(II)=0.
         NPAR(II)=0.
      ENDDO

C____________________________________________________
C LOOP OVER PARTICLE INSIDE GRID
      DO L=1,NP

C________________________________________
C UNBOUNDED PARTICLES: SKIP ASSIGNMENT
         IF(Bflag(L).ne.1) goto 99

         II=IP(L)
           
         IF(II.LE.IL) THEN      !INSIDE GRID
            NPAR(II)=NPAR(II)+1
            V2=VX(L)*VX(L)+VY(L)*VY(L)+VZ(L)*VZ(L) !V*V
            VR=(VX(L)*COS(PH(L))+VY(L)*SIN(PH(L)))*SIN(TH(L))+VZ(L) !RADIALVEL 
     &*COS(TH(L)) 
            VR2=VR*VR           
            VT2=V2-VR2          !TANGENTIAL VEL ^2

c_________________________



            SR(II)=SR(II)+VR2
            ST(II)=ST(II)+VT2
         ELSE
         ENDIF 


 99      continue

      ENDDO                     !END LOOP OVER PART


C____________________________________________________
C COMPUTE KINETIC TANGENTIAL ENERGY  & KINETIC RADIAL ENERGY

      KKR=0.                    !INIT 
      KKT=0.
      DO II=1,IL                !LOOP OVER CELLS

         IF(NPAR(II).GT.0) THEN 
   
            KKR=KKR+SR(II)      !RADIAL & TANGENTIAL KINETIC ENERGY
            KKT=KKT+ST(II)
   
            SR(II)=SR(II)/FLOAT(NPAR(II)) !NORMALIZATION
            ST(II)=ST(II)/FLOAT(NPAR(II)) 

            IF (SR(II).LT.1.0E-10) SR(II)=1.0E-10 !USED TO AVOID INF
            
            ANIS(II)=ST(II)/SR(II) !ANISOTROPY PARAMETER
             
         ELSE
            ANIS(II)=0.
         ENDIF


      ENDDO



C____________________________________________________
C     NORMALIZATION LOOP
         KKR=.5*KKR*MTOT/FLOAT(NP)
         KKT=.5*KKT*MTOT/FLOAT(NP)
         
      
  

C_____________________________________________________
C WRITE OUTPUT
      write(34,1098) time,KKR,KKT,2.*KKR/KKT
      
      write(33,*) time, ' TIME ___R ___VR ___VT ___ANIS ___'
      DO II=1,IL 
         WRITE(33,1098) R(II),SR(II),ST(II),ANIS(II)
      ENDDO

 1098 FORMAT(F12.6,1X,F12.6,1X,F12.6,1X,F12.6)
      
      RETURN
      END



