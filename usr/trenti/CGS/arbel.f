CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CGS - Collisionless Galactic Simulator - Source Files  C  
C                                                         C
C  Created by M. Trenti & T.van Albada in Fortran77 2003  C
C                                                         C
C  Version 1.0 Alpha                                      C
C                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c
c modified to spherically average the clumpy distribution if desidered
c

C
C  FROM van Albada code, modified by M.Stiavelli, adapted for CGS by M.Trenti
C
C     GENERATE AN ARBITRARY ELLIPSOID      


C_________________________________________________________
      SUBROUTINE ARBEL (A,B,C,RAND,CIR,ANVEL,seed,MODE)  
      
      IMPLICIT NONE

C__________________________________________________________
C INPUT
      REAL*8 A,B,C              !AXIS SCALE VALUES
      REAL*8 RAND               !RATIO ENERGYRANDOM/WPOT
      REAL*8 CIR                !RATION ENERGYROTATION/WPOT
      REAL*8 ANVEL(3)           !ANGULAR VELOCITY FOR ROTATION   
      INTEGER SEED              !SEED FOR RANDOM NUMBER GENERATOR
      INTEGER MODE              !OPTION FLAG: 0 --> UNIFORM ELLIPSOID
C                               !             1 --> NON OMOGENEUS ELLIPSOID   

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'


C____________________________________________________
C LOCAL VARIABLES

      REAL*8 E1,E2              !AXIS RATIOS SQUARED
      REAL*8 RAD                !RADIUS OF THE ELLIPSOID
      REAL*8 SIG,SIG2           !VELOCITY DISPERSION & SQUARE
      REAL*8 OME2               !ROTATIONAL VELOCITY SQUARE
      REAL*8 RV(3)              !VELOCITY AS SAMPLED FROM A GAUSSIAN
      REAL*8 XX,YY,ZZ,AA        !AUX POSITIONS OF PARTICLE & "ELLIPTICAL" RADIUS 
      REAL*8 A0                 !AUX FOR RAD
      REAL*8 P,U                !RANDOM NUMBERS
      REAL*8 GP                 !SOME DISTRIBUTION CUT ?
      REAL*8 GV,V               !SOME DISTRIBUTION CUT, SAMPLED VELOCITY
      REAL*8 VPHI               !SOLID BODY VELOCITY ROTATION ALONG Z AXIS
      REAL*8 CUT                !CUT POINT IN VELOCITY DISTRIBUTION
      
      INTEGER L                 !INDEX FOR PARTICLES LOOP
      INTEGER I                 !INDEX FOR LOOP 

c random number common 
        integer idum
        common /ised/ idum

        real ran2
        external ran2

C     MODIFIED FOR IBM BY M. STIAVELLI                                  

        write(*,*) A,B,C,RAND,CIR,ANVEL,seed,MODE


c random number init
        idum=seed
  
cc
      E1=B*B/(A*A)                                                      
      E2=C*C/(A*A)                                                      
      RAD=(A*B*C)**(1.0/3.0)                                            
C                   RAND IS RATIO ERAND/EPOT ; CIR IS RATIO ECIR/EPOT.  
      SIG2=0.4*(GRAVC*MTOT/RAD)*RAND                                    
      SIG =SQRT(ABS(SIG2))                                              
      OME2=3.0*CIR*(GRAVC*MTOT/RAD**3)                                  



C______________________________________________________________________
C                    GENERATE PARTICLES IN  ELLIPSOID ,  RHO(A) AS A**-1
      DO L=1,NP

    1 XX=2.0*RAN2()-1.0  
      YY=2.0*RAN2()-1.0  
      ZZ=2.0*RAN2()-1.0  
      AA=SQRT(XX*XX+YY*YY/E1+ZZ*ZZ/E2)                                  
      IF (AA.GT.1.0) GO TO 1
                                            
      IF (MODE.EQ.0) THEN                                               
C             HOMOGENEOUS ELLIPSOID                                     
         A0=RAD                                                         
         AA=1.                                                          
      ELSE                                                              
C             ARBITRARY (E.G. NON-HOMOGENEOUS)                          
C                   CHOOSE VALUE FOR SEMI MAJOR AXIS, A0                
C            RHO(A0) AS IN CLUMPY SPHERE                                
   20 P=RAN2()     
      IF (MODE.EQ.1) GOTO 25                                            
      GP=(EXP(-((P-0.5)**2/0.10246))-0.08716)/0.91284                   
      U=RAN2()     
      IF (U.GT.GP) GO TO 20                                             
   25 A0=P*RAD                                                          
      ENDIF                                                             

C_________________________________________________________________
C                   SCALE COORDINATES XX,YY,ZZ TO VALUE OF A0           
      X(L)=A0/AA*XX                                                     
      Y(L)=A0/AA*YY                                                     
      Z(L)=A0/AA*ZZ                                                     

C                    COORDINATE TRANSFORMATION                          
      S(L)=SQRT(X(L)*X(L)+Y(L)*Y(L)+Z(L)*Z(L))                          
	IF(S(L).EQ.0.) WRITE(2,*) 'S(L) = 0 FOR L=',L
c
      TH(L)=ACOS(Z(L)/S(L))                                             
      PH(L)=ATAN2(Y(L),X(L))                                            
      IF (PH(L).LT.0.0) PH(L)=PH(L)+2.0*PI                              
C                    SOLID BODY ROTATION (Z-AXIS)                       
      VPHI=SQRT(OME2*(X(L)**2+Y(L)**2))                                 
      VX(L)=-VPHI*SIN(PH(L))                                            
      VY(L)= VPHI*COS(PH(L))                                            
      VZ(L)=0.0                                                         
C                     RANDOM VELOCITY COMPONENTS                        
C                     GAUSS DISTRIBUTION                                
C                     G(V)=EXP(-0.5*V**2/SIG2)                          
      CUT=3.0*SIG                                                       
C                                                                       
      DO I=1,3                                                       
 2       V=CUT*(2.0*RAN2()-1.0)      
         U=RAN2()                    
         GV=EXP(-V**2/(2.0*SIG2))                                          
         IF (U.GT.GV) GO TO 2                                              
C     
         RV(I)=V                                                           
   10 ENDDO                                                          
C                     INITIAL VELOCITIES                                
      VX(L)=VX(L)+RV(1) *ANVEL(1)                                       
      VY(L)=VY(L)+RV(2) *ANVEL(2)                                       
      VZ(L)=VZ(L)+RV(3) *ANVEL(3)                                       
C________________________________________________________
C                   DUMMY-VALUE FOR IP               
      IP(L)=1                                                                                                 
     
      ENDDO                     !END LOOP OVER PARTICLES

      RETURN                                                            
      END                                                               



C_____________________________________________________________
C*************************************************************

C
C  FROM van Albada code, modified by M.Stiavelli, adapted for CGS by M.Trenti
C
C  modified 08/04/04 flag_sm feature added
C
C     GENERATE AN ARBITRARY CLUMPSY ELLIPSOID      


C     ***************************************************************
      SUBROUTINE ARBEL2(A,B,C,REKPAR,REKCL,SEED,NCL,RCL,flag_sm)
C     ***************************************************************
C               GENERATE AN ARBITRARY CLUMPSY ELLIPSOID

C            REKPAR IS RATIO ERAND/EPOT FOR PARTICLES
C            REKCL  IS RATIO ERAND/EPOT FOR CLUMPS
C            EPOT IS TOTAL POTENTIAL ENERGY
C            NCL IS NUMBER OF CLUMPS (EVEN) , MAX VALUE 100
C            NP MUST BE INTEGER MULTIPLE OF NCL
C            RCL IS RADIUS CLUMP (SAME VALUE FOR ALL CLUMPS)
C            DETERMINE NUMBER OF STARS PER CLUMP FROM
C            UNIFORM DISTRIBUTION FROM 0 TO 2*MNS
C            MNS IS MEAN NUMBER OF STARS PER CLUMP
C            flag_sm is a flag to spherically average the positions
C              and velocities of the particle for comparison runs    

C__________________________________________________________
C INPUT
      REAL*8 A,B,C              !AXIS SCALE VALUES
      REAL*8 REKPAR             !RATIO ENERGYRANDOM/WPOT
      REAL*8 REKCL              !RATIO ENERGYRANDOM_clumps/WPOT  
      INTEGER SEED              !SEED FOR RANDOM NUMBER GENERATOR
      INTEGER NCL               !NUMBER OF CLUMPS  
      REAL*8 RCL                !RADIUS OF A CLUMP (SAME FOR ALL CLUMPS)
      INTEGER flag_sm           !flag for smoothing (1=true) 
C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

c random number common 
      integer idum
      common /ised/ idum
      
      double precision ran2
      external ran2
      
      DIMENSION RV(3),XCL(1000),YCL(1000),ZCL(1000),
     1     VXCL(1000),VYCL(1000),VZCL(1000),NSCL(1000)
      

c_________________________
c additional variable for smoothing
      double precision u1,u2,u3,ss,v


c     random number init
      idum=seed
      
      WRITE(*,*) A,B,C,REKPAR,REKCL,SEED,NCL,RCL
      
C     AVERAGE PARTICLE MASS
      AM=MTOT/FLOAT(NP)
C     
      E1=B*B/(A*A)
      E2=C*C/(A*A)
      RAD=(A*B*C)**(1.0/3.0)
      MNS=NP/NCL
      NCLH=NCL/2
      DO 10 I=1,NCLH
      ANS=FLOAT(I-1)*FLOAT(MNS)/FLOAT(NCLH) + 1.0E-7
      INS=INT(ANS)
      IF (INS.GT.(MNS-2)) INS=MNS-2
      NSCL(2*I-1)=MNS+INS
      NSCL(2*I)  =MNS-INS
   10 CONTINUE
C            DETERMINE COORDINATES CLUMP CENTERS
C            UNIFORM DISTRIBUTION WITHIN ELLIPSOID
      DO 20 I=1,NCL
   19 XX=2.0*RAN2()-1.0
      YY=2.0*RAN2()-1.0
      ZZ=2.0*RAN2()-1.0
      AA=SQRT(XX*XX+YY*YY/E1+ZZ*ZZ/E2)
      IF (AA.GT.1.0) GO TO 19
      XCL(I)=XX*(A-RCL)
      YCL(I)=YY*(B-RCL)
      ZCL(I)=ZZ*(C-RCL)
   20 CONTINUE
C            DETERMINE VELOCITIES CLUMP CENTERS
C            ISOTROPIC GAUSSIAN
      SIG2=0.4*(GRAVC*MTOT/RAD)*REKCL
      SIG=SQRT(SIG2)
      CUT=3.0*SIG
      DO 30 I=1,NCL
      DO 29 J=1,3
   28 V=CUT*(2.0*RAN2()-1.0)
      U=RAN2()
      GV=EXP(-V**2/(2.0*SIG2))
      IF (U.GT.GV) GO TO 28
      RV(J)=V
   29 CONTINUE
      VXCL(I)=RV(1)
      VYCL(I)=RV(2)
      VZCL(I)=RV(3)
   30 CONTINUE
      WRITE(2,2011) (XCL(I),YCL(I),ZCL(I),VXCL(I),VYCL(I),VZCL(I),      
     1               NSCL(I),I=1,NCL)                                   
 2011 FORMAT(1H ,6F12.4,I8)                                             
C            DETERMINE POSITIONS OF PARTICLES                           
C            L  COUNTS TOTAL NUMBER OF PARTICLES                        
      L=1                                                               
      DO 200 I=1,NCL                                                    
      NN=NSCL(I)                                                        
C            PUT NN PARTICLES IN CLUMP                                  
      DO 190 N=1,NN                                                     
  180 XX=2.0*RAN2()-1.0
      YY=2.0*RAN2()-1.0
      ZZ=2.0*RAN2()-1.0
      AA=SQRT(XX*XX+YY*YY+ZZ*ZZ)                                        
      IF (AA.GT.1.0) GO TO 180     
      X(L)=XX*RCL + XCL(I)                                             
      Y(L)=YY*RCL + YCL(I)                                             
      Z(L)=ZZ*RCL + ZCL(I)                                  
      VX(L)=VXCL(I)                                                    
      VY(L)=VYCL(I)                                                    
      VZ(L)=VZCL(I)                                                    
      L=L+1                                                             
      IF ((L-1).EQ.(NP)) GO TO 210                                 
  190 CONTINUE                                                          
  200 CONTINUE                                                          
C                                                                       
  210 CONTINUE             

   
C            COORDINATE TRANSFORMATION                    
      DO 300 L=1,NP
      S(L)=SQRT(X(L)*X(L)+Y(L)*Y(L)+Z(L)*Z(L))                          
      TH(L)=ACOS(Z(L)/S(L))                                             
      PH(L)=ATAN2(Y(L),X(L))                                            
      IF (PH(L).LT.0.0) PH(L)=PH(L)+2.0*PI                              

C            DUMMY-VALUE IP                           
      IP(L)=1    


C            DETERMINE VELOCITIES OF PARTICLES                           
C            ISOTROPIC GAUSSIAN                                         
      SIG2=0.4*(GRAVC*MTOT/RAD)*REKPAR                                  
      SIG=SQRT(ABS(SIG2))                                               
      CUT=3.0*SIG                                                       
      DO 220 I=1,3                                                      
  219 V=CUT*(2.0*RAN2()-1.0)
      U=RAN2()
      GV=EXP(-V**2/(2.0*SIG2))                                          
      IF (U.GT.GV) GO TO 219                                            
      RV(I)=V                                                           
  220 CONTINUE                                                          
      VX(L)=VX(L)+RV(1)                                                 
      VY(L)=VY(L)+RV(2)                                                 
      VZ(L)=VZ(L)+RV(3)                                                 
  300 CONTINUE                                                          
C  



C____________________________________
C positions generated: look for smoothing flag
C the procedure is at the end to have exactly the same positions and velocities
C as in the clumpy case from the random number generator
C additional random are drawn here at the end

      if(flag_sm.eq.1) then
         do l=1,NP
         
C
 22         U1= 2.0*RAN2()-1.0
            U2= 2.0*RAN2()-1.0
            U3= 2.0*RAN2()-1.0
            SS= SQRT(U1*U1+U2*U2+U3*U3)
            IF ((SS.GT.1.0).OR.(SS.LT.0.01))GO TO 22
C     
C       **DIRECTION COSINES ACCEPTED**
C     
            X(L)= S(L)*U1/SS
            Y(L)= S(L)*U2/SS
            Z(L)= S(L)*U3/SS 
            TH(L)=ACOS(Z(L)/S(L))
            PH(L)=ATAN2(Y(L),X(L))
            IF (PH(L).LT.0.0) PH(L)=PH(L)+2.0*PI
            
C-------------- velocities
                V=sqrt(VX(L)*VX(L)+VY(L)*VY(L)+VZ(L)*VZ(L))
C
C       **DETERMINE VX,VY,VZ**
C
 44             U1= 2.0*RAN2()-1
                U2= 2.0*RAN2()-1
                U3= 2.0*RAN2()-1
                SS= SQRT(U1*U1+U2*U2+U3*U3)
                IF ((SS.GT.1.0).OR.(SS.LT.0.01)) GO TO 44
C
C       **DIRECTION COSINES ACCEPTED**
C
                VX(L)= V*U1/SS
                VY(L)= V*U2/SS
                VZ(L)= V*U3/SS
 

         enddo

      else
      endif
C-------------------------------------
                                         

c___________debug
      do L=1,NP
         write(68,*) L,S(L),sqrt(VX(L)*VX(L)+VY(L)*VY(L)+VZ(L)*VZ(L))
      enddo

                                                                     
      RETURN                                                            
      END
