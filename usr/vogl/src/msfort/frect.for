	interface to subroutine C_rect[c](a, b, c, d)
	real *8 a, b, c, d
	end

	interface to subroutine C_rects[c](a, b, c, d)
	integer *2 a, b, c, d
	end

	interface to subroutine C_recti[c](a, b, c, d)
	integer a, b, c, d
	end

	interface to subroutine C_rectf[c](a, b, c, d)
	real *8 a, b, c, d
	end

	interface to subroutine C_rectfs[c](a, b, c, d)
	integer *2 a, b, c, d
	end

	interface to subroutine C_rectfi[c](a, b, c, d)
	integer a, b, c, d
	end

	subroutine rect(a, b, c, d)
	call C_rect(a, b, c, d)
	end

	subroutine rects(a, b, c, d)
	integer *2 a, b, c, d
	call C_rects(a, b, c, d)
	end

	subroutine recti(a, b, c, d)
	integer a, b, c, d
	call C_recti(a, b, c, d)
	end

	subroutine rectf(a, b, c, d)
	call C_rectf(a, b, c, d)
	end

	subroutine rectfs(a, b, c, d)
	integer *2 a, b, c, d
	call C_rectfs(a, b, c, d)
	end

	subroutine rectfi(a, b, c, d)
	integer a, b, c, d
	call C_rectfi(a, b, c, d)
	end

