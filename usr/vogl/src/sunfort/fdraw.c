#include "vogl.h"

/*
 * draw_
 */
void
draw_(x, y, z)
	float	*x, *y, *z;
{
	draw(*x, *y, *z);
}

/*
 * draws_
 */
void
draws_(x, y, z)
	short	*x, *y, *z;
{
	draw((float)*x, (float)*y, (float)*z);
}

/*
 * drawi_
 */
void
drawi_(x, y, z)
	int	*x, *y, *z;
{
	draw((float)*x, (float)*y, (float)*z);
}

/*
 * draw2_
 */
void
draw2_(x, y)
	float	*x, *y;
{
	draw(*x, *y, 0.0);
}

/*
 * draw2s_
 */
void
draw2s_(x, y)
	short	*x, *y;
{
	draw((float)*x, (float)*y, 0.0);
}

/*
 * draw2i_
 */
void
draw2i_(x, y)
	int	*x, *y;
{
	draw((float)*x, (float)*y, 0.0);
}

/*
 * rdr_
 */
void
rdr_(dx, dy, dz)
	float	*dx, *dy, *dz;
{
	rdr(*dx, *dy, *dz);
}

/*
 * rdrs_
 */
void
rdrs_(dx, dy, dz)
	short	*dx, *dy, *dz;
{
	rdr((float)*dx, (float)*dy, (float)*dz);
}

/*
 * rdri_
 */
void
rdri_(dx, dy, dz)
	int	*dx, *dy, *dz;
{
	rdr((float)*dx, (float)*dy, (float)*dz);
}

/*
 * rdr2_
 */
void
rdr2_(dx, dy)
	float	*dx, *dy;
{
	rdr2(*dx, *dy);
}

/*
 * rdr2s_
 */
void
rdr2s_(dx, dy)
	short	*dx, *dy;
{
	rdr2((float)*dx, (float)*dy);
}

/*
 * rdr2i_
 */
void
rdr2i_(dx, dy)
	int	*dx, *dy;
{
	rdr2((float)*dx, (float)*dy);
}
