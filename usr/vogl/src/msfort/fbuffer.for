	interface to subroutine C_backbuffer[c](i)
	integer *2 i
	end

	interface to subroutine C_frontbuffer[c](i)
	integer *2 i
	end

	interface to subroutine C_swapbuffers[c]()
	end

	interface to subroutine C_doublebuffer[c]()
	end

	interface to subroutine C_singlebuffer[c]()
	end

	subroutine backbuffer(i)
	call C_backbuffer(i)
	end

	subroutine frontbuffer(i)
	call C_frontbuffer(i)
	end

	subroutine swapbuffers
	call C_swapbuffers
	end

	subroutine doublebuffer
	call C_doublebuffer
	end

	subroutine singlebuffer
	call C_singlebuffer
	end

	subroutine backbu(i)
	call C_backbuffer(i)
	end

	subroutine frontb(i)
	call C_frontbuffer(i)
	end

	subroutine swapbu
	call C_swapbuffers
	end

	subroutine double
	call C_doublebuffer
	end

	subroutine single
	call C_singlebuffer
	end
