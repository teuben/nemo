	interface to subroutine C_deflinestyle[c](n, ls)
	integer n, ls
	end

	interface to subroutine C_setlinestyle[c](n)
	integer n
	end

	interface to subroutine C_linewidth[c](n)
	integer n
	end

	subroutine setlinestyle(n)
	integer n
	call C_setlinestyle(n)
	end

	subroutine deflinestyle(n, ls)
	integer n, ls
	call C_deflinestyle(n, ls)
	end

	subroutine setlin(n)
	integer n
	call C_setlinestyle(n)
	end

	subroutine deflin(n, ls)
	integer n, ls
	call C_deflinestyle(n, ls)
	end

	subroutine linewi(n)
	integer n
	call C_linewidth(n)
	end

	subroutine linewidth(n)
	integer n
	call C_linewidth(n)
	end
