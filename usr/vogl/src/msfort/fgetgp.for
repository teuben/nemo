	interface to subroutine C_getgp[c](x[reference],
     +                              y[reference],
     +                              z[reference])
	end

	interface to subroutine C_getgpos[c](x[reference],
     +                              y[reference],
     +                              z[reference])
	end

	subroutine getgp(x, y, z)
	call C_getgp(x, y, z)
	end

	subroutine getgpos(x, y, z)
	call C_getgpos(x, y, z)
	end

	subroutine getgpo(x, y, z)
	call C_getgpos(x, y, z)
	end
