	interface to subroutine C_pushviewport[c]()
	end

	interface to subroutine C_popviewport[c]()
	end

	interface to subroutine C_viewport[c](a, b, c, d)
	integer *2 a, b, c, d
	end

	interface to subroutine C_getviewport[c](
     +  a[reference], b[reference], c[reference], d[reference])
	integer *2 a, b, c, d
	end

	subroutine pushviewport
	call C_pushviewport
	end

	subroutine pushvi
	call C_pushviewport
	end

	subroutine popviewport
	call C_popviewport
	end

	subroutine popvie
	call C_popviewport
	end

	subroutine viewport(a, b, c, d)
	integer a, b, c, d
	call C_viewport(a, b, c, d)
	end

	subroutine viewpo(a, b, c, d)
	integer a, b, c, d
	call C_viewport(a, b, c, d)
	end

	subroutine getviewport(a, b, c, d)
	integer *2 a, b, c, d
	call C_getviewport(a, b, c, d)
	end

	subroutine getvie(a, b, c, d)
	integer *2 a, b, c, d
	call C_getviewport(a, b, c, d)
	end

