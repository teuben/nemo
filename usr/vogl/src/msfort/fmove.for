	interface to subroutine C_move[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_moves[c](x, y, z)
	integer *2 x, y, z
	end

	interface to subroutine C_movei[c](x, y, z)
	integer x, y, z
	end

	interface to subroutine C_move2[c](x, y)
	real *8 x, y
	end

	interface to subroutine C_move2s[c](x, y)
	integer *2 x, y
	end

	interface to subroutine C_move2i[c](x, y)
	integer x, y
	end

	interface to subroutine C_rmv[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_rmvs[c](x, y, z)
	integer *2 x, y, z
	end

	interface to subroutine C_rmvi[c](x, y, z)
	integer x, y, z
	end

	interface to subroutine C_rmv2[c](x, y)
	real *8 x, y
	end

	interface to subroutine C_rmv2s[c](x, y)
	integer *2 x, y
	end

	interface to subroutine C_rmv2i[c](x, y)
	integer x, y
	end

	subroutine move(x, y, z)
	call C_move(x, y, z)
	end

	subroutine moves(x, y, z)
	integer *2 x, y, z
	call C_move(x, y, z)
	end

	subroutine movei(x, y, z)
	integer x, y, z
	call C_move(x, y, z)
	end

	subroutine move2(x, y)
	call C_move2(x, y)
	end

	subroutine move2s(x, y)
	integer *2 x, y
	call C_move2s(x, y)
	end

	subroutine move2i(x, y)
	integer x, y
	call C_move2i(x, y)
	end

	subroutine rmv(x, y, z)
	call C_rmv(x, y, z)
	end

	subroutine rmvs(x, y, z)
	integer *2 x, y, z
	call C_rmvs(x, y, z)
	end

	subroutine rmvi(x, y, z)
	integer x, y, z
	call C_rmvi(x, y, z)
	end

	subroutine rmv2(x, y)
	call C_rmv2(x, y)
	end

	subroutine rmv2s(x, y)
	integer *2 x, y
	call C_rmv2s(x, y)
	end

	subroutine rmv2i(x, y)
	integer x, y
	call C_rmv2i(x, y)
	end
