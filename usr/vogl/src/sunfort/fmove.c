#include "vogl.h"

/*
 * move_
 */
void
move_(x, y, z)
	float 	*x, *y, *z;
{
	move(*x, *y, *z);
}

/*
 * moves_
 */
void
moves_(x, y, z)
	short 	*x, *y, *z;
{
	move((float)*x, (float)*y, (float)*z);
}

/*
 * movei_
 */
void
movei_(x, y, z)
	int 	*x, *y, *z;
{
	move((float)*x, (float)*y, (float)*z);
}
/*
 * move2_
 */
void
move2_(x, y)
	float	*x, *y;
{
	move(*x, *y, 0.0);
}

/*
 * move2s_
 */
void
move2s_(x, y)
	short	*x, *y;
{
	move((float)*x, (float)*y, 0.0);
}

/*
 * move2i_
 */
void
move2i_(x, y)
	int	*x, *y;
{
	move((float)*x, (float)*y, 0.0);
}

/*
 * rmv_
 */
void
rmv_(dx, dy, dz)
	float	*dx, *dy, *dz;
{
	rmv(*dx, *dy, *dz);
}

/*
 * rmvs_
 */
void
rmvs_(dx, dy, dz)
	short	*dx, *dy, *dz;
{
	rmv((float)*dx, (float)*dy, (float)*dz);
}

/*
 * rmvi_
 */
void
rmvi_(dx, dy, dz)
	int	*dx, *dy, *dz;
{
	rmv((float)*dx, (float)*dy, (float)*dz);
}

/*
 * rmv2_
 */
void
rmv2_(dx, dy)
	float	*dx, *dy;
{
	rmv2(*dx, *dy);
}

/*
 * rmv2s_
 */
void
rmv2s_(dx, dy)
	short	*dx, *dy;
{
	rmv2((float)*dx, (float)*dy);
}

/*
 * rmv2i_
 */
void
rmv2i_(dx, dy)
	int	*dx, *dy;
{
	rmv2((float)*dx, (float)*dy);
}
