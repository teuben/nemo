c
c     experimental subroutine to write the particle positions,
c     velocities, forces & potential calling diagnostic (this is
c     necessary to evaluate forces and shift particles of half a time-step
c
      subroutine WritePart1()

      IMPLICIT NONE

      call diagnostic(10)

      return
      end



c
c     experimental subroutine to write the particle positions and
c     velocities at diagnostic time
c
      subroutine WritePart()

      IMPLICIT NONE

C____________________________________________________
C COMMON FILE WITH SOME PARAMETERS DECLARATION & CONSTANTS
      INCLUDE 'common.blk'

C____________________________________________________
C COMMON FOR TIME
      INCLUDE 'com_time.blk'

C____________________________________________________
C COMMON FOR PARTICLES DATA
      INCLUDE 'com_part.blk'

C______________________________________
C     local variables
      integer I

      write(90,110) NP,time
      DO I=1,NP
         write(90,111) I,X(I),Y(I),Z(I),VX(I),VY(I),VZ(I)
      ENDDO


 110  format(I8,1X,E14.7)
 111  format(I7,6(1X,E13.5))


      return
      end
