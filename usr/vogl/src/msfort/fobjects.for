	interface to subroutine C_makeobj[c](n)
	integer n
	end

	interface to subroutine C_closeobj[c]()
	end

	interface to subroutine C_delobj[c](n)
	integer n
	end

	interface to function C_genobj[c]()
	integer C_genobj
	end

	interface to function C_getopenobj[c]()
	integer C_getopenobj
	end

	interface to subroutine C_callobj[c](n)
	integer n
	end

	interface to function C_isobj[c](n)
	integer n
	integer *2 C_isobj
	end

	subroutine makeobj(n)
	call C_makeobj(n)
	end

	subroutine makeob(n)
	call C_makeobj(n)
	end

	subroutine closeobj
	call C_closeobj
	end

	subroutine closeo
	call C_closeobj
	end

	subroutine delobj(n)
	call C_delobj(n)
	end

	integer function genobj()
	integer C_genobj
	genobj = C_genobj()
	end

	integer function getopenobj()
	integer C_getopenobj
	getopenobj = C_getopenobj()
	end

	integer function getope()
	integer C_getopenobj
	getope = C_getopenobj()
	end

	subroutine callobj(n)
	call C_callobj(n)
	end

	subroutine callob(n)
	call C_callobj(n)
	end

	integer function isobj(n)
	integer *2  C_isobj
	isobj = C_isobj(n)
	end
