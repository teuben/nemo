
#include "vogl.h"

static	long	px = -1, py = -1, pxs = -1, pys = -1;

/*
 * prefposition
 *
 *	Specify a prefered position for a window that is
 *	under control of a window manager.
 *	Position is the location of the upper left corner.
 *	Should be called before ginit.
 */
void
prefposition(x1, x2, y1, y2) 
	long	x1, x2, y1, y2;
{
	if (x1 < 0 || x2 < 0)
		verror("prefposition: bad x value");

	if (y1 < 0 || y2 < 0)
		verror("prefposition: bad y value");


	px = x1;
	py = y1;
	pxs = x2 - x1;
	pys = y2 - y1;
	if (pxs <= 0 || pys <= 0)
		verror("prefposition: bad window size");
}

/*
 * prefsize
 *
 *	Specify the prefered size for a window under control of
 *	a window manager.
 *	Should be called before ginit.
 */
void
prefsize(x, y)
	long	x, y;
{
	if (x < 0)
		verror("prefsize: bad x value");

	if (y < 0)
		verror("prefsize: bad y value");

	pxs = x;
	pys = y;
}

/*
 * getprefposandsize
 *
 *	Returns the prefered position and size of a window under
 *	control of a window manager. (-1 for unset parameters)
 */
void
getprefposandsize(x, y, xs, ys)
	int	*x, *y, *xs, *ys;
{
	*x = px;
	*y = py;
	*xs = pxs;
	*ys = pys;
}

