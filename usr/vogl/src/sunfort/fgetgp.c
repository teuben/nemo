#include "vogl.h"

/*
 * getgp_
 */
void
getgp_(x, y, z)
	float	*x, *y, *z;
{
	getgp(x, y, z);
}

/*
 * getgpos_
 */
void
getgpos_(x, y, z, w)
	float	*x, *y, *z, *w;
{
	getgpos(x, y, z, w);
}

/*
 * getgpo_	(same as getgpos_)
 */
void
getgpo_(x, y, z, w)
	float	*x, *y, *z, *w;
{
	getgpos(x, y, z, w);
}
