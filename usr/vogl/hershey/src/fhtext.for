
      interface to subroutine C_hfont[c](name[reference])
      character *(*) name
      end

      interface to subroutine C_htextsize[c](x, y)
      real *8 x, y
      end

      interface to subroutine C_hboxtext[c](x, y, a, b,
     +  s[reference])
      real*8	x, y, a, b
      character *(*)	s
      end
      
      interface to subroutine C_hboxfit[c](x, y, n)
      real *8 x, y
      integer *2 n
      end

      interface to subroutine C_htextang[c](x)
      real *8 x
      end
      
      interface to subroutine	C_hdrawchar[c](c)
      integer *2 c
      end

      interface to subroutine C_hcharstr[c](str[reference])
      character *(*) str
      end
      
      interface to subroutine C_hgetfontsize[c](x[reference],
     +                                      y[reference])
      real x, y
      end

      interface to function C_hgetfontheight[c]()
      real C_hgetfontheight
      end

      interface to function C_hgetfontwidth[c]()
      real C_hgetfontwidth
      end

      interface to function C_hgetdecender[c]()
      real C_hgetdecender
      end

      interface to function C_hgetascender[c]()
      real C_hgetascender
      end

      interface to subroutine C_hgetcharsize[c](c, x[reference],
     +                                      y[reference])
      character *1 c
      real x, y
      end

      interface to subroutine C_hfixedwidth[c](onoff)
      integer *2 onoff
      end

      interface to subroutine C_hcentertext[c](onoff)
      integer *2 onoff
      end

      interface to subroutine C_hleftjustify[c](onoff)
      integer *2 onoff
      end

      interface to subroutine C_hrightjustify[c](onoff)
      integer *2 onoff
      end

      interface to function C_hnumchars[c]()
      integer *2 C_hnumchars
      end

      interface to function C_hstrlength[c](str[reference])
      real C_hstrlength
      character *(*) str
      end
      
      subroutine hfont(name, l)
      character *(*) name
      character *128 buf

	buf = name
	buf(l + 1 : l + 1) = char(0)
	print*,'buf = ', buf
      call C_hfont(buf)
      return 
      end
      
      subroutine htextsize(x, y)
      call C_htextsize(x, y)
      end

      subroutine htexts(x, y)
      call C_htextsize(x, y)
      end

      subroutine hboxtext(a, b, c, d, name, l)
      character *(*) name
      character *128 buf
	buf = name
	buf(l + 1 : l + 1) = char(0)
      call C_hboxtext(a, b, c, d, buf)
      return 
      end

      subroutine hboxte(a, b, c, d, name, l)
      character *(*) name
      character *128 buf
	buf = name
	buf(l + 1 : l + 1) = char(0)
      call C_hboxtext(a, b, c, d, buf)
      return 
      end
      
      subroutine hboxfit(x, y, n)
      call C_hboxfit(x, y, n)
      end

      subroutine hboxfi(x, y, n)
      call C_hboxfit(x, y, n)
      end
      
      subroutine htextang(x)
      call C_htextang(x)
      end

      subroutine htexta(x)
      call C_htextang(x)
      end

      subroutine hdrawchar(c)
      character *1 c
      call C_hdrawchar(ichar(c))
      end

      subroutine hdrawc(c)
      character *1 c
      call C_hdrawchar(ichar(c))
      end
      
      subroutine hcharstr(name, l)
      character *(*) name
      character *128 buf

	buf = name
	buf(l + 1 : l + 1) = char(0)
      call C_hcharstr(buf)
      end
      
      subroutine hchars(name, l)
      character *(*) name
      character *128 buf

	buf = name
	buf(l + 1 : l + 1) = char(0)
      call C_hcharstr(buf)
      end

      subroutine hgetfontsize(x, y)
      call C_hgetfontsize(x, y)
      end

      subroutine hgetfs(x, y)
      call C_hgetfontsize(x, y)
      end


      subroutine hgetcharsize(c, x, y)
      character *1 c
      call C_hgetcharsize(c, x, y)
      end

      subroutine hgetcs(c, x, y)
      character *1 c
      call C_hgetcharsize(c, x, y)
      end

      real function hgetfontwidth()
	real C_hgetfontwidth
      hgetfontwidth = C_hgetfontwidth()
      end

      real function hgetfw()
	real C_hgetfontwidth
      hgetfw = C_hgetfontwidth()
      end

      real function hgetfontheight()
	real C_hgetfontheight
      hgetfontheight = C_hgetfontheight()
      end

      real function hgetfh()
	real C_hgetfontheight
      hgetfh = C_hgetfontheight()
      end

      real function hgetascender()
	real C_hgetascender
      hgetascender = C_hgetascender()
      end

      real function hgetas()
	real C_hgetascender
      hgetas = C_hgetascender()
      end

      real function hgetdecender()
	real C_hgetdecender
      hgetdecender = C_hgetdecender()
      end

      real function hgetde()
	real C_hgetdecender
      hgetde = C_hgetdecender()
      end

      subroutine hfixedwidth(onoff)
      integer onoff
      call C_hfixedwidth(onoff)
      end

      subroutine hfixed(onoff)
      integer onoff
      call C_hfixedwidth(onoff)
      end

      subroutine hcentertext(onoff)
      integer onoff
      call C_hcentertext(onoff)
      end

      subroutine hcente(onoff)
      integer onoff
      call C_hcentertext(onoff)
      end

      subroutine hleftjustify(onoff)
      integer onoff
      call C_hleftjustify(onoff)
      end

      subroutine hleftj(onoff)
      integer onoff
      call C_hleftjustify(onoff)
      end

      subroutine hrightjustify(onoff)
      integer onoff
      call C_hrightjustify(onoff)
      end

      subroutine hright(onoff)
      integer onoff
      call C_hrightjustify(onoff)
      end

      integer function hnumchars()
      integer*2 C_hnumchars
      hnumchars = C_hnumchars()
      end

      integer function hnumch()
      integer*2 C_hnumchars
      hnumch = C_hnumchars()
      end
      
      real function hstrlength(name, l)
      character *(*) name
      character *128 buf
      real C_hstrlength

	buf = name
	buf(l + 1 : l + 1) = char(0)
      hstrlength = C_hstrlength(buf)
      return
      end
      
      real function hstrle(name, l)
      character *(*) name
      character *128 buf
      real C_hstrlength

	buf = name
	buf(l + 1 : l + 1) = char(0)
      hstrle = C_hstrlength(buf)
      return
      end

