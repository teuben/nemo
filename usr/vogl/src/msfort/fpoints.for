	interface to subroutine C_pnt[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_pnts[c](x, y, z)
	integer *2 x, y, z
	end

	interface to subroutine C_pnti[c](x, y, z)
	integer x, y, z
	end

	interface to subroutine C_pnt2[c](x, y)
	real *8 x, y
	end

	interface to subroutine C_pnt2s[c](x, y)
	integer *2 x, y
	end

	interface to subroutine C_pnt2i[c](x, y)
	integer x, y
	end

	subroutine pnt(x, y, z)
	call C_pnt(x, y, z)
	end

	subroutine pnts(x, y, z)
	integer *2 x, y, z
	call C_pnts(x, y, z)
	end

	subroutine pnti(x, y, z)
	integer x, y, z
	call C_pnti(x, y, z)
	end

	subroutine pnt2(x, y)
	call C_pnt2(x, y)
	end

	subroutine pnt2s(x, y)
	integer *2 x, y
	call C_pnt2s(x, y)
	end

	subroutine pnt2i(x, y)
	integer x, y
	call C_pnt2i(x, y)
	end
