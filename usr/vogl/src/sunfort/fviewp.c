#include "vogl.h"

/*
 * pushviewport_
 */
pushviewport_()
{
	pushviewport();
}

/*
 * pushvi_	(same as pushviewport_)
 */
pushvi_()
{
	pushviewport();
}

/*
 * popviewport_
 */
popviewport_()
{
	popviewport();
}

/*
 * popvie_	(same pushviewport_)
 */
popvie_()
{
	popviewport();
}

/*
 * viewport_
 */
viewport_(xlow, xhigh, ylow, yhigh)
	int	*xlow, *ylow, *xhigh, *yhigh;
{
	viewport(*xlow, *xhigh, *ylow, *yhigh);
}

/*
 * viewpo_	(same as viewport_)
 */
viewpo_(xlow, xhigh, ylow, yhigh)
	int	*xlow, *ylow, *xhigh, *yhigh;
{
	viewport(*xlow, *xhigh, *ylow, *yhigh);
}

/*
 * getviewport_
 */
getviewport_(left, right, bottom, top)
	short	*left, *right, *bottom, *top;
{
	getviewport(left, right, bottom, top);
}

/*
 * getvie_
 */
getvie_(left, right, bottom, top)
	short	*left, *right, *bottom, *top;
{
	getviewport(left, right, bottom, top);
}

