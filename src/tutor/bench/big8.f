      program bigi
c
c  see how big integers*8 can grow, or even if the compiler allows it
c

      integer*8 i

      i=1

 100  continue
          i = i*2
          write(*,*) i
          if (i.gt.0) goto 100
      end
