*
*             F I R S T N 
*             ***********
*
*
*       Basic von Hoerner code (Z.f.Astrophys. 50, 184, 1960).
*       ------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge (12/2000).
*       ......................................................
*
      PROGRAM FIRSTN
*
      INCLUDE 'common1.h'
      EXTERNAL SCALE
*
*
*       Initialize the timer.
*     CALL CPUTIM(CPU0) 
*
*       Read input parameters and perform initial setup.
      CALL INPUT
      CALL DATA
      CALL SCALE
*
*       Evaluate total energy and produce output.
    1 CALL OUTPUT
*
*       Advance solutions until next output.
      CALL INTGRT
      GO TO 1
*
      END
