#include "vogl.h"

/*
 * pnt_
 */
void
pnt_(x, y, z)
	float 	*x, *y, *z;
{
	pnt((Coord)*x, (Coord)*y, (Coord)*z);
}

/*
 * pnts_
 */
void
pnts_(x, y, z)
	short 	*x, *y, *z;
{
	pnt((Coord)*x, (Coord)*y, (Coord)*z);
}


/*
 * pnti_
 */
void
pnti_(x, y, z)
	int 	*x, *y, *z;
{
	pnt((Coord)*x, (Coord)*y, (Coord)*z);
}

/*
 * pnt2_
 */
void
pnt2_(x, y)
	float	*x, *y;
{
	pnt((Coord)*x, (Coord)*y, 0.0);
}

/*
 * pnt2s_
 */
void
pnt2s_(x, y)
	short	*x, *y;
{
	pnt((Coord)*x, (Coord)*y, 0.0);
}

/*
 * pnt2i_
 */
void
pnt2i_(x, y)
	int	*x, *y;
{
	pnt((Coord)*x, (Coord)*y, 0.0);
}
