	interface to subroutine C_polarview[c](a, b, c, d)
	real*8 a
	integer *2 b, c, d
	end

	interface to subroutine C_lookat[c](a, b, c, d, e, f, g)
	real*8 a, b, c, d, e, f
	integer *2 g
	end

	interface to subroutine C_perspective[c](a, b, c, d)
	integer *2 a
	real*8 b, c, d
	end

	interface to subroutine C_window[c](a, b, c, d, e, f)
	real*8 a, b, c, d, e, f
	end

	interface to subroutine C_ortho[c](a, b, c, d, e, f)
	real*8 a, b, c, d, e, f
	end

	interface to subroutine C_ortho2[c](a, b, c, d)
	real*8 a, b, c, d
	end

	subroutine polarview(a, b, c, d)
	real a
	integer b, c, d
	call C_polarview(a, b, c, d)
	end

	subroutine polarv(a, b, c, d)
	real a
	integer b, c, d
	call C_polarview(a, b, c, d)
	end

	subroutine lookat(a, b, c, d, e, f, g)
	integer g
	call C_lookat(a, b, c, d, e, f, g)
	end

	subroutine perspective(a, b, c, d)
	integer a
	call C_perspective(a, b, c, d)
	end

	subroutine perspe(a, b, c, d)
	integer a
	call C_perspective(a, b, c, d)
	end

	subroutine window(a, b, c, d, e, f)
	call C_window(a, b, c, d, e, f)
	end

	subroutine ortho(a, b, c, d, e, f)
	call C_ortho(a, b, c, d, e, f)
	end

	subroutine ortho2(a, b, c, d)
	call C_ortho2(a, b, c, d)
	end
