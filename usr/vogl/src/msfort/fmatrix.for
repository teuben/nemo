	interface to subroutine C_pushmatrix[c]()
	end

	interface to subroutine C_popmatrix[c]()
	end

	interface to subroutine C_multmatrix[c](m[reference])
	real m(1,1)
	end

	interface to subroutine C_getmatrix[c](m[reference])
	real m(1,1)
	end

	interface to subroutine C_loadmatrix[c](m[reference])
	real m(1,1)
	end

	subroutine pushmatrix
	call C_pushmatrix
	end

	subroutine pushma
	call C_pushmatrix
	end

	subroutine popmatrix
	call C_popmatrix
	end

	subroutine popmat
	call C_popmatrix
	end

	subroutine multmatrix(m)
	real m(4, 4)
	call C_multmatrix(m)
	end

	subroutine multma(m)
	real m(4, 4)
	call C_multmatrix(m)
	end

	subroutine getmatrix(m)
	real m(4, 4)
	call C_getmatrix(m)
	end

	subroutine getmat(m)
	real m(4, 4)
	call C_getmatrix(m)
	end

	subroutine loadmatrix(m)
	real m(4, 4)
	call C_loadmatrix(m)
	end

	subroutine loadma(m)
	real m(4, 4)
	call C_loadmatrix(m)
	end
