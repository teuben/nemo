	interface to subroutine C_prefposition[c](x, y, a, b)
	integer x, y, a, b
	end

	interface to subroutine C_prefsize[c](x, y)
	integer x, y
	end

	subroutine prefposition(x, y, a, b)
	integer x, y, a, b
	call C_prefposition(x, y, a, b)
	end

	subroutine prefpo(x, y, a, b)
	integer x, y, a, b
	call C_prefposition(x, y, a, b)
	end

	subroutine prefsize(x, y)
	integer x, y
	call C_prefsize(x, y)
	end

	subroutine prefsi(x, y)
	integer x, y
	call C_prefsize(x, y)
	end
