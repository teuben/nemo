	interface to subroutine C_defbasis[c](id, b[reference])
	integer *2 id
	real b(1, 1)
	end

	interface to subroutine C_patchbasis[c](b1, b2)
	integer b1, b2
	end

	interface to subroutine C_patchprecision[c](n1, n2)
	integer n1, n2
	end

	interface to subroutine C_patchcurves[c](n1, n2)
	integer n1, n2
	end

	interface to subroutine C_rpatch[c](geomx[reference],
     +                                    geomy[reference],
     +                                    geomz[reference],
     +                                    geomw[reference])
	real geomx(1, 1), geomy(1, 1), geomz(1, 1), geomw(1, 1)
	end

	interface to subroutine C_patch[c](geomx[reference],
     +                                    geomy[reference],
     +                                    geomz[reference])
	real geomx(1, 1), geomy(1, 1), geomz(1, 1)
	end

	subroutine patchbasis(basis1, basis2)
	integer basis1, basis2
	call C_patchbasis(basis1, basis2)
	end

	subroutine patchb(basis1, basis2)
	integer basis1, basis2
	call C_patchbasis(basis1, basis2)
	end

	subroutine patchprecision(n1, n2)
	call C_patchprecision(n1, n2)
	end

	subroutine patchp(n1, n2)
	call C_patchprecision(n1, n2)
	end

	subroutine patchcurves(n1, n2)
	call C_patchcurves(n1, n2)
	end

	subroutine patchc(n1, n2)
	call C_patchcurves(n1, n2)
	end

	subroutine rpatch(geomx, geomy, geomz, geomw)
	real geomx(4, 4),geomy(4, 4),geomz(4, 4),geomw(4, 4)
	call C_rpatch(geomx, geomy, geomz, geomw)
	end

	subroutine patch(geomx, geomy, geomz)
	real geomx(4, 4),geomy(4, 4),geomz(4, 4)
	call C_patch(geomx, geomy, geomz)
	end

	subroutine defbasis(id, b)
	integer id
	real b(4, 4)
	call C_defbasis(id, b)
	end

	subroutine defbas(id, b)
	integer id
	real b(4, 4)
	call C_defbasis(id, b)
	end
