CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CGS - Collisionless Galactic Simulator - Source Files  C  
C                                                         C
C  Created by M. Trenti & T.van Albada in Fortran77 2003  C
C                                                         C
C  Version 1.0 Alpha                                      C
C                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C__________________________________________________________
C ADAPTIVE TIME STEP SUBROUTINE

      SUBROUTINE TIMESTEP()

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

      REAL*8 VR                 !RADIAL VELOCITY

      REAL*8 MAXDT              !MAXIMUM ALLOWED TIMESTEP
      REAL*8 AUXDT              !AUX STORING TIME TO CROSS SHELL

      INTEGER L                 !INDEX RUNNING OVER PARTICLES
      INTEGER II                !CELL INDEX


C_________________________________________________
      MAXDT=MDT                 !MAX ALLOWED DT!
C_______________________________________________________
C LOOP OVER PARTICLE
      DO L=1,NP
         II=IP(L)
           
         IF(II.LE.IL) THEN      !INSIDE GRID
            VR=(VX(L)*COS(PH(L))+VY(L)*SIN(PH(L)))*SIN(TH(L))+VZ(L) !RADIALVEL 
     &*COS(TH(L))


            AUXDT=DR(II)/ABS(VR)     !TIME USED TO CROSS THE CELL
            IF(AUXDT.LT.MAXDT) MAXDT=AUXDT !NEW MAX. TIME STEP
         ELSE
         ENDIF 
      ENDDO                     !END LOOP OVER PART


C_________________________________________________
C CHECK FOR MINUMUM ALLOWED TIMESTEP
      IF(MAXDT.LE.MINDT) MAXDT=MINDT
      
C_________________________________________________
      DT=MAXDT                  !NEWTIMESTEP ASSIGNMENT


      WRITE(*,*) 'NEW TIME STEP =', MAXDT

      RETURN
      END
