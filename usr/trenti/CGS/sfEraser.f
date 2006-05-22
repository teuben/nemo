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


C 
C SUBROUTINE SELF FORCES ERASES
C ERASE PARTICLE-ITSELF FORCE INDUCED BY THE GRID

C___________________________________________________
      SUBROUTINE SFER()

      IMPLICIT NONE
C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

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
      REAL VOLShell(NRGMAX)         !VOLUME OF EACH CELL
      
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
C MATRIX VERSION OF VECTORIZED COEFFICIENTS OF SPHE HAR EXPANSION
      real*8 AUXFGA(lmax,lmax,NRGMAX),AUXFGB(lmax,lmax,NRGMAX)
     1   ,mtC(lmax,lmax,NRGMAX), mtD(lmax,lmax,NRGMAX)

C______________________________________________________
C ACCELERATION COMPUTATION
      REAL*8 F1                 !RADIAL FORCE
C
      REAL*8 SUMM_A,SUMM_B      !AUX VARIABLES FOR SUMMATION OVER
C                               !SPHERICAL HARM. EXPANSION
      REAL*8 SUMM1,SUMM2        !FIRST AND SECOND CONTRIBUTION 
C                               !FROM INTERPOLATION
      REAL*8 SS                 !POSITION AT WHICH THE SELF FORCE IS COMPUTED

C______________________________________________________
      INTEGER JJ                !INDEX FOR LOOP OVER EACH CELL SUBDIVISIONS

C____________________________________________________
C INIT STUFF:

C______________________________
      onePmass=MTOT/float(NP)   !single particle mass



C____________________________________________________
      DO II=1,IL1               !LOOP OVER CELLS
         DO  JJ=1,NDIV+1
            DI=FLOAT(JJ-1)/FLOAT(NDIV) !NORMALIZED DISTANCE FROM CELL
            SS=RR(II)+DI*(RR(II+1)-RR(II))

C____________________________________________________
         DO I=1,IL              !SPHE. HAR. COEFF.
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
         do i=1,IL1             !RADIAL SHELL VOLUME
            VOLShell(I)=(RR(I+1)**3-RR(I)**3)/3.
         enddo
         
C     end of init stuff
         
C______________________________________________________
C     COEFFICIENTS ANM AND BNM
         

C______________________________________________________
c     generates legendre polynomials for the particle
         
         do N=1,lmax
            do M=1,lmax
               plm(N,M)=0.d0
            enddo
         enddo
         cth=.1D0               !COS THETA
       
C______________________________________________________
C LEGENDRE POLYNOMIALS FOR EACH PARTICLE
         call  plgndrMT(NHAR,CTH) 
         
         
C______________________________________________________
C generate the array containing cos(m phi) and sin(m phi)
         do i=1,nhar+1
            cmphi(i)=1.D0
            smphi(i)=0.D0
         enddo
         
C______________________________________________________
C     RADIAL CELL ASSIGNMENT
         

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

 

C________________________________________________________________
C ACCELERATION

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
      DI=FLOAT(JJ-1)/FLOAT(NDIV) !NORMALIZED DISTANCE FROM CELL


C______________________________________________________
C van Albada recipe for assignment to the two nearest cells
         IF (DI.GT.0.5) GO TO 31
         I0=II
         I1=II+1
         GO TO 32
 31      DI=1.0-DI
         I0=II+1
         I1=II
 32      CONTINUE

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
         if(SS.le.eps) SS=eps
         SF(II,JJ)=PI4G/SS*(SUMM1*(1.-DI)+SUMM2*DI) !RADIAL FORCE
C                                                         ! ASSIGNMENT

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
         
         
         SEP(II,JJ)=(SUMM1*(1.-DI)+SUMM2*DI) !POTENTIAL ENERGY ASSIGNMENT     

C___________________________________________-
C     test
         SEP(II,JJ) = 0.0
C________________________________________________________
C     TEST
         SF(II,JJ) = 0.0

         



      ENDDO                     !END LOOP OVER DI FOR EACH CELL
      ENDDO                     !END LOOP OVER CELLS


      RETURN
      END
