	interface to subroutine C_polymode[c](mode)
	integer mode
	end

	interface to subroutine C_poly2[c](n , p[reference])
	integer n
	real p(1, 1)
	end

	interface to subroutine C_poly2s[c](n , p[reference])
	integer n
	integer *2  p(1, 1)
	end

	interface to subroutine C_poly2i[c](n , p[reference])
	integer n
	integer p(1, 1)
	end


	interface to subroutine C_poly[c](n , p[reference])
	integer n
	real p(1, 1)
	end

	interface to subroutine C_polys[c](n , p[reference])
	integer n
	integer *2 p(1, 1)
	end

	interface to subroutine C_polyi[c](n , p[reference])
	integer n
	integer p(1, 1)
	end

	interface to subroutine C_polf2[c](n , p[reference])
	integer n
	real p(1, 1)
	end

	interface to subroutine C_polf2s[c](n , p[reference])
	integer n
	integer *2  p(1, 1)
	end

	interface to subroutine C_polf2i[c](n , p[reference])
	integer n
	integer p(1, 1)
	end


	interface to subroutine C_polf[c](n , p[reference])
	integer n
	real p(1, 1)
	end

	interface to subroutine C_polfs[c](n , p[reference])
	integer n
	integer *2 p(1, 1)
	end

	interface to subroutine C_polfi[c](n , p[reference])
	integer n
	integer p(1, 1)
	end

	interface to subroutine C_pmv[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_pmvs[c](x, y, z)
	integer *2 x, y, z
	end

	interface to subroutine C_pmvi[c](x, y, z)
	integer x, y, z
	end

	interface to subroutine C_pmv2[c](x, y)
	real *8 x, y
	end

	interface to subroutine C_pmv2s[c](x, y)
	integer *2 x, y
	end

	interface to subroutine C_pmv2i[c](x, y)
	integer x, y
	end

	interface to subroutine C_rpmv[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_rpmvs[c](x, y, z)
	integer *2 x, y, z
	end

	interface to subroutine C_rpmvi[c](x, y, z)
	integer x, y, z
	end

	interface to subroutine C_rpmv2[c](x, y)
	real *8 x, y
	end

	interface to subroutine C_rpmv2s[c](x, y)
	integer *2 x, y
	end

	interface to subroutine C_rpmv2i[c](x, y)
	integer x, y
	end

	interface to subroutine C_pdr[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_pdrs[c](x, y, z)
	integer *2 x, y, z
	end

	interface to subroutine C_pdri[c](x, y, z)
	integer x, y, z
	end

	interface to subroutine C_pdr2[c](x, y)
	real *8 x, y
	end

	interface to subroutine C_pdr2s[c](x, y)
	integer *2 x, y
	end

	interface to subroutine C_pdr2i[c](x, y)
	integer x, y
	end

	interface to subroutine C_rpdr[c](x, y, z)
	real *8 x, y, z
	end

	interface to subroutine C_rpdrs[c](x, y, z)
	integer *2 x, y, z
	end

	interface to subroutine C_rpdri[c](x, y, z)
	integer x, y, z
	end

	interface to subroutine C_rpdr2[c](x, y)
	real *8 x, y
	end

	interface to subroutine C_rpdr2s[c](x, y)
	integer *2 x, y
	end

	interface to subroutine C_rpdr2i[c](x, y)
	integer x, y
	end

	interface to subroutine C_backface[c](onoff)
	integer *2 onoff
	end

	interface to subroutine C_frontface[c](onoff)
	integer *2 onoff
	end

	interface to subroutine C_pclos[c]()
	end

	subroutine polymode(mode)
	integer mode
	call C_polymode(mode)
	end

	subroutine polymo(mode)
	integer mode
	call C_polymode(mode)
	end

	subroutine backface(n)
	call C_backface(n)
	end

	subroutine backfa(n)
	call C_backface(n)
	end

	subroutine frontface(n)
	call C_frontface(n)
	end

	subroutine frontf(n)
	call C_frontface(n)
	end

	subroutine poly2(n, p)
	real p(1, 1)
	call C_poly2(n , p)
	end

	subroutine poly2s(n, p)
	integer *2 p(1, 1)
	call C_poly2s(n , p)
	end

	subroutine poly2i(n, p)
	integer p(1, 1)
	call C_poly2i(n , p)
	end

	subroutine poly(n, p)
	real p(1, 1)
	call C_poly(n , p)
	end

	subroutine polys(n, p)
	integer *2 p(1, 1)
	call C_polys(n , p)
	end

	subroutine polyi(n, p)
	integer p(1, 1)
	call C_polyi(n , p)
	end

	subroutine polf2(n, p)
	real p(1, 1)
	call C_polf2(n , p)
	end

	subroutine polf2s(n, p)
	integer *2 p(1, 1)
	call C_polf2s(n , p)
	end

	subroutine polf2i(n, p)
	integer p(1, 1)
	call C_polf2i(n , p)
	end

	subroutine polf(n, p)
	real p(1, 1)
	call C_polf(n , p)
	end

	subroutine polfs(n, p)
	integer *2 p(1, 1)
	call C_polfs(n , p)
	end

	subroutine polfi(n, p)
	integer p(1, 1)
	call C_polfi(n , p)
	end

	subroutine pmv(x, y, z)
	call C_pmv(x, y, z)
	end

	subroutine pmvs(x, y, z)
	integer *2 x, y, z
	call C_pmvs(x, y, z)
	end

	subroutine pmvi(x, y, z)
	integer x, y, z
	call C_pmvi(x, y, z)
	end

	subroutine pmv2(x, y)
	call C_pmv2(x, y)
	end

	subroutine pmv2s(x, y)
	integer *2 x, y
	call C_pmv2s(x, y)
	end

	subroutine pmv2i(x, y)
	integer x, y
	call C_pmv2i(x, y)
	end
	
	subroutine rpmv(x, y, z)
	call C_rpmv(x, y, z)
	end

	subroutine rpmvs(x, y, z)
	integer *2 x, y, z
	call C_rpmvs(x, y, z)
	end

	subroutine rpmvi(x, y, z)
	integer x, y, z
	call C_rpmvi(x, y, z)
	end

	subroutine rpmv2(x, y)
	call C_rpmv2(x, y)
	end

	subroutine rpmv2s(x, y)
	integer *2 x, y
	call C_rpmv2s(x, y)
	end

	subroutine rpmv2i(x, y)
	integer x, y
	call C_rpmv2i(x, y)
	end

	subroutine pdr(x, y, z)
	call C_pdr(x, y, z)
	end

	subroutine pdrs(x, y, z)
	integer *2 x, y, z
	call C_pdrs(x, y, z)
	end

	subroutine pdri(x, y, z)
	integer x, y, z
	call C_pdri(x, y, z)
	end

	subroutine pdr2(x, y)
	call C_pdr2(x, y)
	end

	subroutine pdr2s(x, y)
	integer *2 x, y
	call C_pdr2s(x, y)
	end

	subroutine pdr2i(x, y)
	integer x, y
	call C_pdr2i(x, y)
	end
	
	subroutine rpdr(x, y, z)
	call C_rpdr(x, y, z)
	end

	subroutine rpdrs(x, y, z)
	integer *2 x, y, z
	call C_rpdrs(x, y, z)
	end

	subroutine rpdri(x, y, z)
	integer x, y, z
	call C_rpdri(x, y, z)
	end

	subroutine rpdr2(x, y)
	call C_rpdr2(x, y)
	end

	subroutine rpdr2s(x, y)
	integer *2 x, y
	call C_rpdr2s(x, y)
	end

	subroutine rpdr2i(x, y)
	integer x, y
	call C_rpdr2i(x, y)
	end
	subroutine pclos
	call C_pclos
	end
