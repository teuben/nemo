	interface to subroutine C_arc[c](x, y, r, s, e)
	real *8 x, y, r	
	integer*2 s, e
	end

	interface to subroutine C_arcs[c](x, y, r, s, e)
	integer *2 x, y, r
	integer*2 s, e
	end

	interface to subroutine C_arci[c](x, y, r, s, e)
	integer x, y, r
	integer*2 s, e
	end

	interface to subroutine C_arcf[c](x, y, r, s, e)
	real *8 x, y, r	
	integer*2 s, e
	end

	interface to subroutine C_arcfs[c](x, y, r, s, e)
	integer *2 x, y, r
	integer*2 s, e
	end

	interface to subroutine C_arcfi[c](x, y, r, s, e)
	integer x, y, r
	integer*2 s, e
	end

	interface to subroutine C_circ[c](x, y, r)
	real *8 x, y, r
	end

	interface to subroutine C_circs[c](x, y, r)
	integer *2 x, y, r
	end

	interface to subroutine C_circi[c](x, y, r)
	integer x, y, r
	end

	interface to subroutine C_circf[c](x, y, r)
	real *8 x, y, r
	end

	interface to subroutine C_circfs[c](x, y, r)
	integer *2 x, y, r
	end

	interface to subroutine C_circfi[c](x, y, r)
	integer x, y, r
	end

	interface to subroutine C_arcprecision[c](n)
	integer *2 n
	end

	interface to subroutine C_circleprecision[c](n)
	integer *2 n
	end

	subroutine arc(x, y, r, s, e)
	real x, y, r
	integer s, e
	call C_arc(x, y, r, s, e)
	end

	subroutine arcs(x, y, r, s, e)
	integer *2  x, y, r
	integer s, e
	call C_arcs(x, y, r, s, e)
	end

	subroutine arci(x, y, r, s, e)
	integer  x, y, r
	integer s, e
	call C_arci(x, y, r, s, e)
	end

	subroutine arcf(x, y, r, s, e)
	real x, y, r
	integer s, e
	call C_arcf(x, y, r, s, e)
	end

	subroutine arcfs(x, y, r, s, e)
	integer *2  x, y, r
	integer s, e
	call C_arcfs(x, y, r, s, e)
	end

	subroutine arcfi(x, y, r, s, e)
	integer  x, y, r
	integer s, e
	call C_arcfi(x, y, r, s, e)
	end

	subroutine circ(x, y, r)
	call C_circ(x, y, r)
	end

	subroutine circs(x, y, r)
	integer *2 x, y, r
	call C_circs(x, y, r)
	end

	subroutine circi(x, y, r)
	integer x, y, r
	call C_circi(x, y, r)
	end

	subroutine circf(x, y, r)
	call C_circf(x, y, r)
	end

	subroutine circfs(x, y, r)
	integer *2 x, y, r
	call C_circfs(x, y, r)
	end

	subroutine circfi(x, y, r)
	integer x, y, r
	call C_circfi(x, y, r)
	end

	subroutine arcprecision(n)
	call C_arcprecision(n)
	end

	subroutine arcpre(n)
	call C_arcprecision(n)
	end

	subroutine circpr(n)
	call C_circleprecision(n)
	end
