#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

/*
 *	An array of points for a polygon
 */
static Coord	parray[][3] = {
	{-8.0, -8.0, 0.0},
	{-5.0, -8.0, 0.0},
	{-5.0, -5.0, 0.0},
	{-8.0, -5.0, 0.0}
};

/*
 * drawpoly
 *
 *	draw some polygons
 */
void
drawpoly()
{
	float	vec[3];
	short	val;

	color(YELLOW);

	/*
	 * Draw a polygon using poly, parray is our array of
	 * points and 4 is the number of points in it.
	 */
	poly(4L, parray);

	color(GREEN);

	/*
	 * Draw a 5 sided figure by using bgnpolygon, v3d, and endpolygon
	 */
	polymode(PYM_LINE);
	bgnpolygon();
		vec[0] = 0.0;
		vec[1] = 0.0;
		vec[2] = 0.0;
		v3f(vec);
		vec[0] = 3.0;
		vec[1] = 0.0;
		vec[2] = 0.0;
		v3f(vec);
		vec[0] = 3.0;
		vec[1] = 4.0;
		vec[2] = 0.0;
		v3f(vec);
		vec[0] = -1.0;
		vec[1] = 5.0;
		vec[2] = 0.0;
		v3f(vec);
		vec[0] = -2.0;
		vec[1] = 2.0;
		vec[2] = 0.0;
		v3f(vec);
	endpolygon();

	color(MAGENTA);

	/*
	 * draw a sector representing a 1/4 circle
	 */
	arc(1.5, -7.0, 3.0, 0, 900);

	move2(1.5, -7.0);
	draw2(1.5, -4.0);

	move2(1.5, -7.0);
	draw2(4.5, -7.0);

	qread(&val);
}

/*
 * drawpolyf
 *
 *	draw some filled polygons
 */
void
drawpolyf()
{
	short	val;

	color(YELLOW);

	polymode(PYM_FILL);

	/*
	 * Draw a polygon using poly, parray is our array of
	 * points and 4 is the number of points in it.
	 */
	polf(4L, parray);

	color(GREEN);

	/*
	 * Draw a filled 5 sided figure by using pmv, pdr and pclos.
	 */
	pmv(0.0, 0.0, 0.0);
		pdr(3.0, 0.0, 0.0);
		pdr(3.0, 4.0, 0.0);
		pdr(-1.0, 5.0, 0.0);
		pdr(-2.0, 2.0, 0.0);
	pclos();

	color(MAGENTA);

	/*
	 * draw a filled sector representing a 1/4 circle
	 */
	arcf(1.5, -7.0, 3.0, 0, 900);

	qread(&val);
}

/*
 * Using polygons, hatching, and filling.
 */
main()
{
	short	val;

	winopen("poly");

	unqdevice(INPUTCHANGE);
	qdevice(KEYBD);		/* enable keyboard */
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	color(BLACK);		/* clear to black */
	clear();

	/*
	 * world coordinates are now in the range -10 to 10
	 * in x, y, and z. Note that positive z is towards us.
	 */
	ortho(-10.0, 10.0, -10.0, 10.0, 10.0, -10.0);

	color(YELLOW);

	/*
	 * write out the string "Polygon from poly()" in the
	 * starting at (-8.0, -4.0) and scaled to be 4.0 units long,
	 * 0.5 units high.
	 */
	hfont("futura.m");
	hboxtext(-8.0, -4.0, 4.0, 0.5, "Polygon from poly()/ polf()");

	color(GREEN);

	/*
	 * write out a scaled string starting at (0.0, 6.0)
	 */
	hboxtext(0.0, 6.0, 4.5, 0.5, "Polygon from bgnpoly()/ endpoly()");
	hboxtext(0.0, 5.0, 4.5, 0.5, "             pmv()/ pdr()/ pclos()");

	color(MAGENTA);

	/*
	 * write out a scaled string starting at (0.0, 6.0)
	 */
	hboxtext(3.5, -3.5, 1.9, 0.5, "Arc/ Arcf");

	/*
	 * draw some wire frame polygons
	 */
	drawpoly();

	/*
	 *  rotate so the next polygons will appear in a different place.
	 */
	rot(20.0, 'x');
	rot(30.0, 'y');

	/*
	 * draw some filled polygons.
	 */
	drawpolyf();

	gexit();
}
