#include "vogl.h"

static	float	Vcx, Vcy, Vsx, Vsy;

/*
 * calcW2Vcoeffs
 *
 *	Calculate the linear coeffs defining the mapping of world
 *	space to actual device space
 */
void
CalcW2Vcoeffs()
{
	Vcx = (float)(vdevice.maxVx + vdevice.minVx) * 0.5;
	Vcy = (float)(vdevice.maxVy + vdevice.minVy) * 0.5;

	Vsx = (float)(vdevice.maxVx - vdevice.minVx) * 0.5;
	Vsy = (float)(vdevice.maxVy - vdevice.minVy) * 0.5;
}

/*
 * WtoVx
 *
 * return the Screen X coordinate corresponding to world point 'p' 
 * (does the perspective division as well)
 */
int
WtoVx(p)
	float	p[];
{
	return((int)(p[0] * Vsx / p[3] + Vcx + 0.5));
}

/*
 * WtoVy
 *
 * return the Screen Y coordinate corresponding to world point 'p' 
 * (does the perspective division as well)
 */
int
WtoVy(p)
	float	p[];
{
	return((int)(p[1] * Vsy / p[3] + Vcy + 0.5));
}
