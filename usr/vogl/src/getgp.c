#include "vogl.h"

/*
 * getgp
 *
 *	return the current (x, y, z) graphics position
 */
void
getgp(x, y, z)
	Coord	*x, *y, *z;
{
	*x = vdevice.cpW[V_X];
	*y = vdevice.cpW[V_Y];
	*z = vdevice.cpW[V_Z];
}

/*
 * getgpos
 *
 *	Get the current graphics position after transformation by
 *	the current matrix.
 */
void
getgpos(x, y, z, w)
	Coord	*x, *y, *z, *w;
{
	/* Make sure that it's all updated for current spot */

	multvector(vdevice.cpWtrans, vdevice.cpW, vdevice.transmat->m);
	*x = vdevice.cpWtrans[V_X];
	*y = vdevice.cpWtrans[V_Y];
	*z = vdevice.cpWtrans[V_Z];
	*w = vdevice.cpWtrans[V_W];
}

