#include "vogl.h"

/*
 * polymode_
 */
void
polymode_(onoff)
	int	*onoff;
{
	polymode((long)*onoff);
}

/*
 * polymo_	(same as polymode_)
 */
void
polymo_(onoff)
	int	*onoff;
{
	polymode((long)*onoff);
}

/*
 * poly2_
 */
void
poly2_(n, parray)
	int	*n;
	float	parray[][2];
{
	poly2((long)*n, parray);
}

/*
 * poly2s_
 */
void
poly2s_(n, parray)
	int	*n;
	short	parray[][2];
{
	poly2s((long)*n, parray);
}

/*
 * poly2i_
 */
void
poly2i_(n, parray)
	int	*n;
	Icoord	parray[][2];
{
	poly2i((long)*n, parray);
}

/*
 * poly_
 */
void
poly_(n, parray)
	int	*n;
	float	parray[][3];
{
	poly((long)*n, parray);
}

/*
 * polys_
 */
void
polys_(n, parray)
	int	*n;
	short	parray[][3];
{
	polys((long)*n, parray);
}

/*
 * polyi_
 */
void
polyi_(n, parray)
	int	*n;
	Icoord	parray[][3];
{
	polyi((long)*n, parray);
}
/*
 * polf_
 */
void
polf_(n, parray)
	int	*n;
	float	parray[][3];
{
	polf((long)*n, parray);
}

/*
 * polfs_
 */
void
polfs_(n, parray)
	int	*n;
	short	parray[][3];
{
	polfs((long)*n, parray);
}

/*
 * polfi_
 */
void
polfi_(n, parray)
	int	*n;
	Icoord	parray[][3];
{
	polfi((long)*n, parray);
}

/*
 * polf2_
 */
void
polf2_(n, parray)
	int	*n;
	float	parray[][2];
{
	polf2((long)*n, parray);
}

/*
 * polf2s_
 */
void
polf2s_(n, parray)
	int	*n;
	short	parray[][2];
{
	polf2s((long)*n, parray);
}

/*
 * polf2i_
 */
void
polf2i_(n, parray)
	int	*n;
	int	parray[][2];
{
	polf2i((long)*n, parray);
}

/*
 * pmv_
 */
void
pmv_(x, y, z)
	float 	*x, *y, *z;
{
	pmv(*x, *y, *z);
}

/*
 * pmvs_
 */
void
pmvs_(x, y, z)
	short 	*x, *y, *z;
{
	pmv((float)*x, (float)*y, (float)*z);
}

/*
 * pmvi_
 */
void
pmvi_(x, y, z)
	int 	*x, *y, *z;
{
	pmv((float)*x, (float)*y, (float)*z);
}

/*
 * pmv2_
 */
void
pmv2_(x, y)
	float	*x, *y;
{
	pmv(*x, *y, 0.0);
}

/*
 * pmv2s_
 */
void
pmv2s_(x, y)
	short	*x, *y;
{
	pmv((float)*x, (float)*y, 0.0);
}

/*
 * pmv2i_
 */
void
pmv2i_(x, y)
	int	*x, *y;
{
	pmv((float)*x, (float)*y, 0.0);
}

/*
 * rpmv_
 */
void
rpmv_(dx, dy, dz)
	float	*dx, *dy, *dz;
{
	rpmv(*dx, *dy, *dz);
}

/*
 * rpmvs_
 */
void
rpmvs_(dx, dy, dz)
	short	*dx, *dy, *dz;
{
	rpmv((float)*dx, (float)*dy, (float)*dz);
}

/*
 * rpmvi_
 */
void
rpmvi_(dx, dy, dz)
	int	*dx, *dy, *dz;
{
	rpmv((float)*dx, (float)*dy, (float)*dz);
}

/*
 * rpmv2_
 */
void
rpmv2_(dx, dy)
	float	*dx, *dy;
{
	rpmv2(*dx, *dy);
}

/*
 * rpmv2s_
 */
void
rpmv2s_(dx, dy)
	short	*dx, *dy;
{
	rpmv2((float)*dx, (float)*dy);
}

/*
 * rpmv2i_
 */
void
rpmv2i_(dx, dy)
	int	*dx, *dy;
{
	rpmv2((float)*dx, (float)*dy);
}


/*
 * pdr_
 */
void
pdr_(x, y, z)
	float 	*x, *y, *z;
{
	pdr(*x, *y, *z);
}

/*
 * pdrs_
 */
void
pdrs_(x, y, z)
	short 	*x, *y, *z;
{
	pdr((float)*x, (float)*y, (float)*z);
}

/*
 * pdri_
 */
void
pdri_(x, y, z)
	int 	*x, *y, *z;
{
	pdr((float)*x, (float)*y, (float)*z);
}

/*
 * pdr2_
 */
void
pdr2_(x, y)
	float	*x, *y;
{
	pdr(*x, *y, 0.0);
}

/*
 * pdr2s_
 */
void
pdr2s_(x, y)
	short	*x, *y;
{
	pdr((float)*x, (float)*y, 0.0);
}

/*
 * pdr2i_
 */
void
pdr2i_(x, y)
	int	*x, *y;
{
	pdr((float)*x, (float)*y, 0.0);
}

/*
 * rpdr_
 */
void
rpdr_(dx, dy, dz)
	float	*dx, *dy, *dz;
{
	rpdr(*dx, *dy, *dz);
}

/*
 * rpdrs_
 */
void
rpdrs_(dx, dy, dz)
	short	*dx, *dy, *dz;
{
	rpdr((float)*dx, (float)*dy, (float)*dz);
}

/*
 * rpdri_
 */
void
rpdri_(dx, dy, dz)
	int	*dx, *dy, *dz;
{
	rpdr((float)*dx, (float)*dy, (float)*dz);
}

/*
 * rpdr2_
 */
void
rpdr2_(dx, dy)
	float	*dx, *dy;
{
	rpdr2(*dx, *dy);
}

/*
 * rpdr2s_
 */
void
rpdr2s_(dx, dy)
	short	*dx, *dy;
{
	rpdr2((float)*dx, (float)*dy);
}

/*
 * rpdr2i_
 */
void
rpdr2i_(dx, dy)
	int	*dx, *dy;
{
	rpdr2((float)*dx, (float)*dy);
}

/*
 * pclos_
 */
void
pclos_()
{
	pclos();
}

/*
 * backface_
 */
void
backface_(onoff)
	int	*onoff;
{
	backface(*onoff);
}

/*
 * backfa_
 */
void
backfa_(onoff)
	int	*onoff;
{
	backface(*onoff);
}

/*
 * frontface_
 */
void
frontface_(onoff)
	int	*onoff;
{
	frontface(*onoff);
}

/*
 * frontf_
 */
void
frontf_(onoff)
	int	*onoff;
{
	frontface(*onoff);
}
