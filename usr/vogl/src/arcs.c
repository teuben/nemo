#include "vogl.h"

#ifdef	TC

extern	double	cos();
extern	double	sin();

#else 

#include <math.h>

#endif

static int	nsegs = 32;

/*
 * arcprecision
 *
 *	sets the number of segments in an arc or circle.
 *	- obsolete function.
 */
void
arcprecision(noseg)
	int	noseg;
{
	nsegs = noseg;
}

/*
 * circleprecision
 *
 *	sets the number of segments in an arc or circle.
 */
void
circleprecision(noseg)
	int	noseg;
{
	nsegs = noseg;
}

/*
 * arc
 *
 * draw an arc at a given location.  Precision of arc (# line segments)
 * is calculated from the value given to circleprecision.
 *
 */
void
arc(x, y, radius, sang, eang)
	Coord	x, y, radius;
	Angle	sang, eang;
{
	Token	*t;
	float	cx, cy, dx, dy;
	float	startang, endang, deltang, cosine, sine, angle;
	int	i, numsegs, sync;

	if (!vdevice.initialised)
		verror("arc: vogl not initialised");


	startang = (float)sang / 10.0;
	endang = (float)eang / 10.0;

	angle = startang * D2R;
	numsegs = (endang - startang) / 360.0 * nsegs + 0.5;
	deltang = (endang - startang) * D2R / numsegs;
	cosine = cos((double)deltang);
	sine = sin((double)deltang);

	if (vdevice.inobject) {
		t = newtokens(8);
		t[0].i = ARC;
		t[1].f = x;
		t[2].f = y;
		t[3].f = radius * cos((double)angle);
		t[4].f = radius * sin((double)angle);
		t[5].f = cosine;
		t[6].f = sine;
		t[7].i = numsegs;
		return;
	}

	if (sync = vdevice.sync)
		vdevice.sync = 0;

	/* calculates initial point on arc */

	cx = x + radius * cos((double)angle);
	cy = y + radius * sin((double)angle);
	move2(cx, cy);

	for (i = 0; i < numsegs; i++)  {
		dx = cx - x; 
		dy = cy - y;
		cx = x + dx * cosine - dy * sine;
		cy = y + dx * sine + dy * cosine;
		draw2(cx, cy);
	}

	if (sync) {
		vdevice.sync = 1;
		(*vdevice.dev.Vsync)();
	}
}

/*
 * arcs
 *
 * draw an arc at a given location.  (Expressed as short integers) 
 * Precision of arc (# line segments) is calculated from the value
 * given to circleprecision.
 *
 */
void
arcs(x, y, radius, sang, eang)
	Scoord	x, y, radius;
	Angle	sang, eang;
{
	arc((Coord)x, (Coord)y, (Coord)radius, sang, eang);
}

/*
 * arci
 *
 * draw an arc at a given location.  (Expressed as integers) 
 * Precision of arc (# line segments) is calculated from the value
 * given to circleprecision.
 *
 */
void
arci(x, y, radius, sang, eang)
	Icoord	x, y, radius;
	Angle	sang, eang;
{
	arc((Coord)x, (Coord)y, (Coord)radius, sang, eang);
}

/*
 * arcf
 *
 *	draw a filled sector in a given location. The number of line
 * segments in the arc of the segment is the same as in arc.
 */
void
arcf(x, y, radius, sang, eang)
	Coord	x, y, radius;
	Angle	sang, eang;
{
	Token	*t;
	float	cx, cy, dx, dy;
	float	deltang, cosine, sine, angle;
	int	i, numsegs;
	float	startang, endang;

	if (!vdevice.initialised)
		verror("arcf: vogl not initialised");

	startang = sang / 10.0;
	endang = eang / 10.0;

	angle = startang * D2R;
	numsegs = (endang - startang) / 360.0 * nsegs + 0.5;
	deltang = (endang - startang) * D2R / numsegs;
	cosine = cos((double)deltang);
	sine = sin((double)deltang);

	if (vdevice.inobject) {
		t = newtokens(8);
		t[0].i = ARCF;
		t[1].f = x;
		t[2].f = y;
		t[3].f = radius * cos((double)angle);
		t[4].f = radius * sin((double)angle);
		t[5].f = cosine;
		t[6].f = sine;
		t[7].i = numsegs;
		return;
	}

	pmv2(x, y);
			/* calculates initial point on arc */

	cx = x + radius * cos((double)angle);
	cy = y + radius * sin((double)angle);

	pdr2(cx, cy);

	for (i = 0; i < numsegs; i++)  {
		dx = cx - x; 
		dy = cy - y;
		cx = x + dx * cosine - dy * sine;
		cy = y + dx * sine + dy * cosine;
		pdr2(cx, cy);
	}

	pclos();
}

/*
 * arcfs
 *
 * draw a filled sector at a given location.  (Expressed as short integers) 
 * Precision of arc (# line segments) is calculated from the value
 * given to circleprecision.
 *
 */
void
arcfs(x, y, radius, sang, eang)
	Scoord	x, y, radius;
	Angle	sang, eang;
{
	arcf((Coord)x, (Coord)y, (Coord)radius, sang, eang);
}

/*
 * arcfi
 *
 * draw a filled sector at a given location.  (Expressed as integers) 
 * Precision of arc (# line segments) is calculated from the value
 * given to circleprecision.
 *
 */
void
arcfi(x, y, radius, sang, eang)
	Icoord	x, y, radius;
	Angle	sang, eang;
{
	arcf((Coord)x, (Coord)y, (Coord)radius, sang, eang);
}


/*
 * circ
 *
 * Draw a circle of given radius at given world coordinates. The number of
 * segments in the circle is the same as that of an arc.
 *
 */
void
circ(x, y, radius)
	Coord	x, y, radius;
{
	Token	*t;
	float	cx, cy, dx, dy;
	float	angle, cosine, sine;
	int	i, sync;

	if (!vdevice.initialised)
		verror("circ: vogl not initialised");

	angle = 2.0 * PI / nsegs;
	cosine = cos((double)angle);
	sine = sin((double)angle);

	if (vdevice.inobject) {
		t = newtokens(7);
		t[0].i = CIRCLE;
		t[1].f = x;
		t[2].f = y;
		t[3].f = radius;
		t[4].f = cosine;
		t[5].f = sine;
		t[6].i = nsegs;
		return;
	}

	cx = x + radius;
	cy = y;

	if (sync = vdevice.sync)
		vdevice.sync = 0;

	move2(cx, cy);
	for (i = 0; i < nsegs; i++) {
		dx = cx - x; 
		dy = cy - y;
		cx = x + dx * cosine - dy * sine;
		cy = y + dx * sine + dy * cosine;
		draw2(cx, cy);
	}

	if (sync) {
		vdevice.sync = 1;
		draw2(x + radius, y);
	}
}

/*
 * circs
 *
 * Draw a circle of given radius at given world coordinates expressed as
 * short integers. The number of segments in the circle is the same as that
 * of an arc.
 *
 */
void
circs(x, y, radius)
	Scoord	x, y, radius;
{
	circ((Coord)x, (Coord)y, (Coord)radius);
}


/*
 * circi
 *
 * Draw a circle of given radius at given world coordinates expressed as
 * integers. The number of segments in the circle is the same as that
 * of an arc.
 *
 */
void
circi(x, y, radius)
	Icoord	x, y, radius;
{
	circ((Coord)x, (Coord)y, (Coord)radius);
}

/*
 * circf
 *
 * Draw a filled circle of given radius at given world coordinates.
 * The number of segments in the circle is the same as that of an arc.
 *
 */
void
circf(x, y, radius)
	Coord	x, y, radius;
{
	Token	*t;
	float	cx, cy, dx, dy;
	float	angle, cosine, sine;
	int	i;

	if (!vdevice.initialised)
		verror("circf: vogl not initialised");

	angle = 2.0 * PI / nsegs;
	cosine = cos((double)angle);
	sine = sin((double)angle);

	if (vdevice.inobject) {
		t = newtokens(7);
		t[0].i = CIRCF;
		t[1].f = x;
		t[2].f = y;
		t[3].f = radius;
		t[4].f = cosine;
		t[5].f = sine;
		t[6].i = nsegs;
		return;
	}

	cx = x + radius;
	cy = y;

	pmv2(cx, cy);
	for (i = 0; i < nsegs; i++) {
		dx = cx - x; 
		dy = cy - y;
		cx = x + dx * cosine - dy * sine;
		cy = y + dx * sine + dy * cosine;
		pdr2(cx, cy);
	}

	pclos();
}

/*
 * circfs
 *
 * Draw a circle of given radius at given world coordinates expressed as
 * short integers. The number of segments in the circle is the same as that
 * of an arc.
 *
 */
void
circfs(x, y, radius)
	Scoord	x, y, radius;
{
	circf((Coord)x, (Coord)y, (Coord)radius);
}

/*
 * circfi
 *
 * Draw a circle of given radius at given world coordinates expressed as
 * integers. The number of segments in the circle is the same as that
 * of an arc.
 *
 */
void
circfi(x, y, radius)
	Icoord	x, y, radius;
{
	circf((Coord)x, (Coord)y, (Coord)radius);
}

