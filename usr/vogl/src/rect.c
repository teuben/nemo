#include "vogl.h"

/*
 * rect
 *
 * draw a rectangle given two opposite corners
 *
 */
void
rect(x1, y1, x2, y2)
	Coord 	x1, y1, x2, y2;
{
	if (!vdevice.initialised)
		verror("rect: vogl not initialised");

	move2(x1, y1);
	draw2(x2, y1);
	draw2(x2, y2);
	draw2(x1, y2);
	draw2(x1, y1);
}

/*
 * recti
 *
 * draw a rectangle given two opposite corners (expressed as integers)
 */
void
recti(x1, y1, x2, y2)
	Icoord 	x1, y1, x2, y2;
{
	rect((Coord)x1, (Coord)y1, (Coord)x2, (Coord)y2);
}

/*
 * rects
 *
 * draw a rectangle given two opposite corners (expressed as short integers)
 */
void
rects(x1, y1, x2, y2)
	Scoord 	x1, y1, x2, y2;
{
	rect((Coord)x1, (Coord)y1, (Coord)x2, (Coord)y2);
}

/*
 * rectf
 *
 * draw a filled rectangle given two opposite corners
 *
 */
void
rectf(x1, y1, x2, y2)
	Coord 	x1, y1, x2, y2;
{
	Token	*tok;

	if (!vdevice.initialised)
		verror("rect: vogl not initialised");

	if (vdevice.inobject) {
		tok = newtokens(5);
		tok[0].i = RECTF;
		tok[1].f = x1;
		tok[2].f = y1;
		tok[3].f = x2;
		tok[4].f = y2;
		return;
	}
		
	pmv2(x1, y1);
	pdr2(x2, y1);
	pdr2(x2, y2);
	pdr2(x1, y2);
	pdr2(x1, y1);
	pclos();
}


/*
 * rectfi
 *
 * draw a filled rectangle given two opposite corners (expressed as integers)
 */
void
rectfi(x1, y1, x2, y2)
	Icoord 	x1, y1, x2, y2;
{
	rectf((Coord)x1, (Coord)y1, (Coord)x2, (Coord)y2);
}

/*
 * rectfs
 *
 * draw a filled rectangle given two opposite corners (expressed as short
 * integers)
 */
void
rectfs(x1, y1, x2, y2)
	Scoord 	x1, y1, x2, y2;
{
	rectf((Coord)x1, (Coord)y1, (Coord)x2, (Coord)y2);
}

