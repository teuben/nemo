#include "vogl.h"
/*
 * vega_draw
 *
 * For V/EGA
 * Draws a line from the current screen spot to x, y
 */

extern	unsigned	_cur_color;

vega_draw(x, y)
	register int	x, y;
{
	egaline(vdevice.cpVx, vdevice.sizeSy - vdevice.cpVy, x, vdevice.sizeSy - y, _cur_color);
	vdevice.cpVx = x;
	vdevice.cpVy = y;

	return(0);
}
