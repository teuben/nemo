      interface to subroutine C_yobbarays[c](onoff)
      integer *2 onoff
      end
      subroutine yobbarays(i)
      call C_yobbarays(i)
      end

      subroutine yobbar(i)
      call C_yobbarays(i)
      end
