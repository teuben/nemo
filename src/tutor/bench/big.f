      program bigi
c
c  see how big integers can grow
c

      integer i

      i=1

 100  continue
          i = i*2
          write(*,*) i
          if (i.gt.0) goto 100
      end
