      interface to subroutine C_draw[c](x, y, z)
      real *8 x, y, z
      end

      interface to subroutine C_draws[c](x, y, z)
      integer *2 x, y, z
      end

      interface to subroutine C_drawi[c](x, y, z)
      integer x, y, z
      end

      interface to subroutine C_draw2[c](x, y)
      real *8 x, y
      end

      interface to subroutine C_draw2s[c](x, y)
      integer *2 x, y
      end

      interface to subroutine C_draw2i[c](x, y)
      integer x, y
      end

      interface to subroutine C_rdr[c](x, y, z)
      real *8 x, y, z
      end

      interface to subroutine C_rdrs[c](x, y, z)
      integer *2 x, y, z
      end

      interface to subroutine C_rdri[c](x, y, z)
      integer x, y, z
      end

      interface to subroutine C_rdr2[c](x, y)
      real *8 x, y
      end

      interface to subroutine C_rdr2s[c](x, y)
      integer *2 x, y
      end

      interface to subroutine C_rdr2i[c](x, y)
      integer x, y
      end

      subroutine draw(x, y, z)
      call C_draw(x, y, z)
      end

      subroutine draws(x, y, z)
      call C_draws(x, y, z)
      end

      subroutine drawi(x, y, z)
      call C_drawi(x, y, z)
      end

      subroutine draw2(x, y)
      call C_draw2(x, y)
      end

      subroutine draw2s(x, y)
      call C_draw2s(x, y)
      end

      subroutine draw2i(x, y)
      call C_draw2i(x, y)
      end

      subroutine rdr(x, y, z)
      call C_rdr(x, y, z)
      end

      subroutine rdrs(x, y, z)
      call C_rdrs(x, y, z)
      end

      subroutine rdri(x, y, z)
      call C_rdri(x, y, z)
      end

      subroutine rdr2(x, y)
      call C_rdr2(x, y)
      end

      subroutine rdr2s(x, y)
      call C_rdr2s(x, y)
      end

      subroutine rdr2i(x, y)
      call C_rdr2i(x, y)
      end
