	interface to subroutine C_font[c](id)
	integer *2 id
	end

	interface to subroutine C_charstr[c](str[reference])
	character *(*) str
	end

	interface to subroutine C_cmov[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_cmovs[c](x, y, z)
	integer *2 x, y, z
	end

	interface to subroutine C_cmovi[c](x, y, z)
	integer x, y, z
	end

	interface to subroutine C_cmov2[c](x, y)
	real *8 x, y
	end

	interface to subroutine C_cmov2s[c](x, y)
	integer *2 x, y
	end

	interface to subroutine C_cmov2i[c](x, y)
	integer x, y
	end

	interface to function C_getwidth[c]()
	integer C_getwidth
	end

	interface to function C_getheight[c]()
	integer C_getheight
	end

	interface to function C_strwidth[c](s[reference])
	integer C_strwidth
	character *(*) s
	end

	interface to subroutine C_getcpos[c](x[reference], 
     +	y[reference])
	integer*2 x, y
	end

	subroutine font(id)
	call C_font(id)
	end

	subroutine charstr(name, l)
	character *(*) name
	integer l
	character *128 buf

	buf = name
	buf(l + 1 : l + 1) = char(0)

	call C_charstr(buf)
	end

	subroutine charst(name, l)
	character *(*) name
	integer l
	character *128 buf

	buf = name
	buf(l + 1 : l + 1) = char(0)

	call C_charstr(buf)
	end
	subroutine cmov(x, y, z)
	call C_cmov(x, y, z)
	end

	subroutine cmovs(x, y, z)
	integer *2 x, y, z
	call C_cmovs(x, y, z)
	end

	subroutine cmovi(x, y, z)
	integer x, y, z
	call C_cmovi(x, y, z)
	end

	subroutine cmov2(x, y)
	call C_cmov2(x, y)
	end

	subroutine cmov2s(x, y)
	integer *2 x, y
	call C_cmov2s(x, y)
	end

	subroutine cmov2i(x, y)
	integer x, y
	call C_cmov2i(x, y)
	end

	integer function getwidth()
	integer C_getwidth
	getwidth = C_getwidth()
	end

	integer function getwid()
	integer C_getwidth
	getwid = C_getwidth()
	end

	integer function getheight()
	integer C_getheight
	getheight = C_getheight()
	end

	integer function gethei()
	integer C_getheight
	gethei = C_getheight()
	end

	integer function strwid(s, l)
	integer C_strwidth
	character *(*) s
	integer l
	character *128 buf

	buf = s
	buf(l + 1 : l + 1) = char(0)

	l = C_strwidth(buf)

	strwid = l
	return
	end

	integer function strwidth(s, l)
	integer C_strwidth
	character *(*) s
	integer l
	character *128 buf

	buf = s
	buf(l + 1 : l + 1) = char(0)

	l = C_strwidth(buf)

	strwidth = l
	return
	end

	subroutine getcpo(ix, iy)
	integer *2 ix, iy

	call C_getcpos(ix, iy)
	end
	
	subroutine getcpos(ix, iy)
	integer *2 ix, iy

	call C_getcpos(ix, iy)
	end

