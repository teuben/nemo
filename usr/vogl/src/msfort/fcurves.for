	interface to subroutine C_curvebasis[c](id)
	integer *2 id
	end

	interface to subroutine C_curveprecision[c](nsegs)
	integer *2 nsegs
	end

	interface to subroutine C_rcrv[c](geom[reference])
	real geom(1,1)
	end

	interface to subroutine C_crv[c](geom[reference])
	real geom(1,1)
	end

	interface to subroutine C_crvn[c](n, 
     + geom[reference])
	integer n
	real geom(1,1)
	end

	interface to subroutine C_rcrvn[c](n, 
     + geom[reference])
	integer n
	real geom(1,1)
	end

	interface to subroutine C_curveit[c](n)
	integer *2 n
	end

	subroutine curvebasis(id)
	call C_curvebasis(id)
	end

	subroutine curveprecision(nsegs)
	call C_curveprecision(nsegs)
	end

	subroutine rcrv(geom)
	real geom(4, 4)
	call C_rcrv(geom)
	end

	subroutine crv(geom)
	real geom(4, 3)
	call C_crv(geom)
	end

	subroutine crvn(n, geom)
	real geom(1, 1)
	call C_crvn(n, geom)
	end

	subroutine rcrvn(n, geom)
	real geom(1, 1)
	call C_rcrvn(n, geom)
	end

	subroutine curveit(n)
	call C_curveit(n)
	end

	subroutine curvei(n)
	call C_curveit(n)
	end

	subroutine curveb(id)
	call C_curvebasis(id)
	end

	subroutine curvep(nsegs)
	call C_curveprecision(nsegs)
	end

