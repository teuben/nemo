	interface to subroutine C_translate[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_scale[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_rotate[c](x, a)
	integer *2 x
	character *1 a
	end

	interface to subroutine C_rot[c](x, a)
	real *8 x
	character *1 a
	end

	subroutine translate(x, y, z)
	call C_translate(x, y, z)
	end

	subroutine transl(x, y, z)
	call C_translate(x, y, z)
	end

	subroutine scale(x, y, z)
	call C_scale(x, y, z)
	end

	subroutine rotate(x, a)
	integer x
	character *1 a
	call C_rotate(x, a)
	end

	subroutine rot(x, a)
	character *1 a
	call C_rot(x, a)
	end

