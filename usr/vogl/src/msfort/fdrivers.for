	interface to subroutine C_voutput[c](s[reference])
	character *(*) s
	end

	interface to subroutine C_vinit[c](s[reference])
	character *(*) s
	end

	interface to subroutine C_ginit[c]()
	end

	interface to subroutine C_winopen[c](s[reference])
	character *(*) s
	end

	interface to subroutine C_vnewdev[c](s[reference])
	character *(*) s
	end

	interface to subroutine C_pushdev[c](s[reference])
	character *(*) s
	end

	interface to subroutine C_popdev[c]()
	end

	interface to subroutine C_vgetdev[c](s[reference])
	character *(*) s
	end

	interface to subroutine C_gexit[c]()
	end

	interface to subroutine C_clear[c]()
	end

	interface to subroutine C_color[c](i)
	integer *2 i
	end

	interface to subroutine C_colorf[c](i)
	real *8 i
	end

	interface to subroutine C_mapcolor[c](i, r, g, b)
	integer *2 i, r, g, b
	end

	interface to function C_getplanes[c]()
	integer C_getplanes
	end

	interface to function C_getvaluator[c](n)
	integer C_getvaluator
	integer n
	end

	interface to function C_getbutton[c](n)
	integer*2 C_getbutton
	integer n
	end

	interface to function C_getgdesc[c](n)
	integer C_getgdesc
	integer n
	end

	interface to subroutine C_foreground[c]()
	end

	interface to subroutine C_vsetflush[c](i)
	integer *2 i
	end

	interface to subroutine C_vflush[c]()
	end


	subroutine voutput(s, n)
	character *(*) s
	character *128 buf

	buf = s
	buf(n + 1 : n + 1) = char(0)

	call C_voutput(buf)
	end

	subroutine voutpu(s, n)
	character *(*) s
	character *128 buf

	buf = s
	buf(n + 1 : n + 1) = char(0)

	call C_voutput(buf)
	end

	subroutine vinit(s, n)
	character *(*) s
	integer n

	character *128 buf

	buf = s
	buf(n + 1 : n + 1) = char(0)

	call C_vinit(buf)
	end

	subroutine ginit()
	call C_ginit()
	end

	subroutine winopen(s, n)
	character *(*) s
	character *128 buf

	buf = s
	buf(n + 1 : n + 1) = char(0)

	call C_winopen(buf)
	end

	subroutine winope(s, n)
	character *(*) s
	character *128 buf

	buf = s
	buf(n + 1 : n + 1) = char(0)

	call C_winopen(buf)
	end

	subroutine vnewdev(s, n)
	character *(*) s
	character *128 buf

	buf = s
	buf(n + 1 : n + 1) = char(0)

	call C_vnewdev(buf)
	end

	subroutine vnewde(s, n)
	character *(*) s
	character *128 buf

	buf = s
	buf(n + 1 : n + 1) = char(0)

	call C_vnewdev(buf)
	end

	subroutine pushdev(s)
	character *(*) s
	character *128 t_string

	call C_pushdev(t_string(s))
	end

	subroutine pushde(s)
	character *(*) s
	character *128 t_string

	call C_pushdev(t_string(s))
	end

	subroutine popdev
	call C_popdev
	end


	subroutine vgetdev(s)
	character *(*) s
	call C_vgetdev(s)
c
c	find terminating null and space fill it all
c
	l = len(s)
	do 10 i = 1, l
	  if (ichar(s(i : i)) .eq. 0) then
	      do 5, j = i, l
		  s(j : j) = ' '
5             continue
	  end if
10      continue
	end

	subroutine vgetde(s)
	character *(*) s
	call C_vgetdev(s)
c
c	find terminating null and space fill it all
c
	l = len(s)
	do 10 i = 1, l
	  if (ichar(s(i : i)) .eq. 0) then
	      do 5, j = i, l
		  s(j : j) = ' '
5             continue
	  end if
10      continue
	end

	subroutine gexit
	call C_gexit
	end

	subroutine clear
	call C_clear
	end

	subroutine color(i)
	call C_color(i)
	end

	subroutine colorf(c)
	call C_color(c)
	end

	subroutine mapcolor(i, r, g, b)
	integer i, r, g, b
	call C_mapcolor(i, r, g, b)
	end

	subroutine mapcol(i, r, g, b)
	integer i, r, g, b
	call C_mapcolor(i, r, g, b)
	end

	integer function getplanes()
	integer C_getplanes
	getplanes = C_getplanes()
	end

	integer function getpla()
	integer C_getplanes
	getpla = C_getplanes()
	end

	integer function getvaluator(n)
	integer C_getvaluator
	getvaluator = C_getvaluator(n)
	end

	integer function getval(n)
	integer C_getvaluator
	getval = C_getvaluator(n)
	end

	integer function getbutton(n)
	integer*2 C_getbutton
	getbutton = C_getbutton(n)
	end

	integer function getbut(n)
	integer*2 C_getbutton
	getbut = C_getbutton(n)
	end

	subroutine gconfig
	end

	subroutine gconfi
	end

	subroutine reshapeviewport
	end

	subroutine reshap
	end

	subroutine winconstraints
	end

	subroutine wincon
	end

	subroutine shademodel(m)
	end

	subroutine shadem(m)
	end

	integer function getgdesc(m)
	integer C_getgdesc
	getgdesc = C_getgdesc(m)
	return
	end

	integer function getgde(m)
	integer C_getgdesc
	getgde = C_getgdesc(m)
	return
	end

	subroutine foreground
	end

	subroutine foregr
	end

	subroutine vsetflush(i)
	call C_vsetflush(i)
	end

	subroutine vsetfl(i)
	call C_vsetflush(i)
	end

	subroutine vflush
	call C_vflush
	end
