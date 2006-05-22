CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CGS - Collisionless Galactic Simulator - Source Files  C  
C                                                         C
C  Created by M. Trenti & T.van Albada in Fortran77 2003  C
C                                                         C
C  Version 1.1                                            C
C                                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C INSERT LICENCE AGREEMENT

C INSERT PAPER ACKNOWLEDGEMENT


C STRUCTURE OF THE CODE:
C
C   MAIN DRIVER
C       |
C       |---INIT PARAMETERS OF SIMULATION  
C       |---INIT PARTICLES   
C       |---INIT GRID
C       |---INIT POISSON SOLVER/ ACCELERATION ASSIGNMENT
C       |
C      (DO)---+
C       ^     |---PARTICLE TO GRID ASSIGNMENT 
C       |     |---POISSON SOLVER    
C       |     |---UPDATE POSITIONS & VELOCITIES 
C       |     |---CHOOSE NEW TIMESTEP
C       |     |
C       |  (IF REQUIRED)|---GRID RECENTERING @ CM OF THE SYSTEM 
C                       |---DIAGNOSTIC SUBROUTINES
C       |     |
C       |  (CHECK EXIT )---STORE PARTICLE DATA FOR RE-RUN
C       |     |
C       +-----+
C

C UNITS IN THE CODE:
C                                                
C  MASS:     10**11   MSUN                             
C  LENGTH:   10       KPC                              
C  TIME:     10**8    YR                               
C  (VELOCITY: 97.8     KM/S)                            
C


C____________________________________________________
      PROGRAM CGS 

      IMPLICIT NONE 

C____________________________________________________
C GLOBAL VARIABLES:

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

C____________________________________________________
C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C COMMON FOR GRID QUANTITIES
      INCLUDE 'com_grid.blk'
 
C SPHERICAL HARMONIC EXPANSION COEFFICIENTS
C___________________________________________________
      INCLUDE 'com_har.blk'

C____________________________________________________
C LOCAL VARIABLES
      INTEGER I                 !INDEX FOR LOOP OVER TIME STEPS
      

C____________________________________________________
C EXECUTABLE STATEMENTS HERE:

      CALL INITPARAMETER()      !PARAMETER INIT
      
      CALL INITPART()           !PARTICLES INIT
      
      CALL INITGRID()           !GRID INIT ACCORDINGLY TO PARAMETERS GIVEN

      CALL ADAPTGRID()          !ADAPT GRID TO THE GIVEN PARTICLE DISTRIBUTION
      
      CALL INITPOISSON()        !CREATE LEGENDRE POL. & TRIG. FUNC. TABLES
 
      CALL SFER()               !SELF FORCES ERASER COMPUTATION

C_____________________________________________________
C FIRST STEP

      CALL P2G()                !ASSIGN PARTICLES TO RADIAL GRID
 
      CALL CMASS()              !CENTER COORDINATES AT THE CENTER OF SYSTEM

      CALL POISSON()            !SOLVE POISSON EQUATION

      CALL INIT_VEL()           !INIT VELOCITIES TO -DT/2 FOR LEAP-FROG

      CALL DIAGNOSTIC(0)        !PERFORM DIAGNOSTIC (ETOT,JTOT,2K/W...)

      CALL TIMESTEP()           !CHECK TIME STEP EXPERIMENTAL mt020403
      
      CALL WritePart1()          !snapshot writing
C____________________________________________________
C NOW BEGIN THE LOOP OVER TIME STEPS
      DO I=1,LSTEP

C____________________________________________________
C WRITE THE STEP NUMBER
         write(*,*) 'Performing step ',I
         
         CALL P2G()             !ASSIGN PARTICLES TO RADIAL GRID

         CALL POISSON()         !SELFCONSISTENT POISSON SOLVER 
         
         CALL ACCEL()           !PARTICLES ARE ADVANCED IN TIME          
         
C____________________________________________________  
         IF(I.EQ.(INPR*(I/INPR))) THEN 
            CALL DIAGNOSTIC(0)  !CHECK FOR DIAGNOSTIC
         ELSE
         ENDIF

         IF(I.EQ.(INSNAP*(I/INSNAP))) THEN 
            CALL WritePart1()    !snapshot writing
         ELSE
         ENDIF

cpjt        CALL FLUSH()     !! intel doesn't like the no-argument here....

         IF(I.EQ.(INCM*(I/INCM))) THEN

            CALL CMASS1()       !CHECK FOR SYSTEM RECENTERING 
            CALL ADAPTGRID()    !ADAPT GRID TO THE GIVEN PARTICLE DISTRIBUTION
            CALL SFER()         !SELF FORCES ERASER COMPUTATION
            CALL CMASS1()       !CHECK FOR SYSTEM RECENTERING JUST ANOTHER CHECK TO BE SURE
            CALL TIMESTEP()     !CHECK TIME STEP (experimental)

         ELSE
         ENDIF

         IF(TIME.GE.TMAX) GOTO 1 !CHECK FOR TIME OF END SIMULATION

      ENDDO

C____________________________________________________
C END SIMULATION PROCEDURES

 1    CALL P2G()                !ASSIGN PARTICLES TO RADIAL GRID

      CALL POISSON()            !SELFCONSISTENT POISSON SOLVER 
      
      CALL DIAGNOSTIC(1)        !FINAL DIAGNOSTIC IS PERFORMED
      
c     CALL ENDSIM()             !STORES PARTICLES DATA FOR RERUN
      
      open(98,file='OK.sim')
      write(98,*) 'Simulation done!'
      close(98)

      END

