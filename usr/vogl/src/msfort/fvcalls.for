	interface to subroutine C_v4d[c](v[reference])
	real*8 v(4)
	end

	interface to subroutine C_v3d[c](v[reference])
	real*8 v(3)
	end

	interface to subroutine C_v2d[c](v[reference])
	real*8 v(2)
	end
	interface to subroutine C_v4f[c](v[reference])
	real v(4)
	end

	interface to subroutine C_v3f[c](v[reference])
	real v(3)
	end

	interface to subroutine C_v2f[c](v[reference])
	real v(2)
	end

	interface to subroutine C_v4s[c](v[reference])
	integer *2 v(4)
	end

	interface to subroutine C_v3s[c](v[reference])
	integer *2 v(3)
	end

	interface to subroutine C_v2s[c](v[reference])
	integer *2 v(2)
	end

	interface to subroutine C_v4i[c](v[reference])
	integer v(4)
	end

	interface to subroutine C_v3i[c](v[reference])
	integer v(3)
	end

	interface to subroutine C_v2i[c](v[reference])
	integer v(2)
	end

	subroutine v4d(v)
	real*8 v(4)
	call C_v4d(v)
	end

	subroutine v3d(v)
	real*8 v(3)
	call C_v3d(v)
	end

	subroutine v2d(v)
	real*8 v(2)
	call C_v2d(v)
	end

	subroutine v4f(v)
	real v(4)
	call C_v4f(v)
	end

	subroutine v3f(v)
	real v(3)
	call C_v3f(v)
	end

	subroutine v2f(v)
	real v(2)
	call C_v2f(v)
	end

	subroutine v4i(v)
	integer v(4)
	call C_v4i(v)
	end

	subroutine v3i(v)
	integer v(3)
	call C_v3i(v)
	end

	subroutine v2i(v)
	integer v(2)
	call C_v2i(v)
	end

	subroutine v4s(v)
	integer *2 v(4)
	call C_v4s(v)
	end

	subroutine v3s(v)
	integer *2 v(3)
	call C_v3s(v)
	end

	subroutine v2s(v)
	integer *2 v(2)
	call C_v2s(v)
	end
