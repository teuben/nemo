#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

/*
 * curve basis types
 */
Matrix	bezier = {
	{-1.0,	3.0,	-3.0,	1.0},
	{3.0,	-6.0,	3.0,	0.0},
	{-3.0,	3.0,	0.0,	0.0},
	{1.0,	0.0,	0.0,	0.0} 
};

Matrix	cardinal = {
	{-0.5,	1.5,	-1.5,	0.5},
	{1.0,	-2.5,	2.0,	-0.5},
	{-0.5,	0.0,	0.5,	0.0},
	{0.0,	1.0,	0.0,	0.0}
};

Matrix	bspline = {
	{-1.0 / 6.0,	3.0 / 6.0,	-3.0 / 6.0,	1.0 / 6.0},
	{3.0 / 6.0,	-6.0 / 6.0,	3.0 / 6.0,	0.0},
	{-3.0 / 6.0,	0.0,		3.0 / 6.0,	0.0},
	{1.0 / 6.0,	4.0 / 6.0,	1.0 / 6.0,	0.0}    
};

/*
 * 	Geometry matrix to demonstrate basic spline segments
 */
float   geom1[4][3] = {
	{ -180.0, 10.0, 0.0 },
	{ -100.0, 110.0, 0.0 },
	{ -100.0, -90.0, 0.0 },
	{ 0.0, 50.0, 0.0 }
};

/*
 * 	Geometry matrix to demonstrate overlapping control points to
 *	produce continuous (Well, except for the bezier ones) curves
 *	from spline segments
 */
float	geom2[6][3] = {
	{ 200.0, 480.0, 0.0 },
	{ 380.0, 180.0, 0.0 },
	{ 250.0, 430.0, 0.0 },
	{ 100.0, 130.0, 0.0 },
	{ 50.0,  280.0, 0.0 },
	{ 150.0, 380.0, 0.0 }
};

/*
 * using curves
 */
main()
{
	char	dev[20];
	int	i;
	short	val;

	vinit("mswin");
	winopen("curves");

	qdevice(KEYBD);
	unqdevice(INPUTCHANGE);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	ortho2(-200.0, 400.0, -100.0, 500.0);

	color(BLACK);
	clear();

	color(YELLOW);

	/*
	 * label the control points in geom1
	 */
        for (i = 0; i < 4; i++) {
		cmov2(geom1[i][0], geom1[i][1]);
		sprintf(dev, "%d", i);
		charstr(dev);
	}
								 
	/*
	 * label the control points in geom2
	 */
	for (i = 0; i < 6; i++) {
		cmov2(geom2[i][0], geom2[i][1]);
		sprintf(dev, "%d", i);
		charstr(dev);
	}

	/*
	 * set the number of line segments appearing in each curve to 20
	 */
	curveprecision((short)20);

	/*
	 * copy the bezier basis matrix into the basis matrix stack and
	 * set the curve basis accordingly.
	 */
	defbasis((short)1, bezier);
	curvebasis((short)1);

	color(RED);

	/*
	 * draw a curve using the current basis matrix (bezier in this case)
	 * and the control points in geom1
	 */
	crv(geom1);

	cmov2(70.0, 60.0);
	charstr("Bezier Curve Segment");

	cmov2(-190.0, 450.0);
	charstr("Three overlapping Bezier Curves");

	/*
	 * crvn draws overlapping curve segments according to geom2, the
	 * number of curve segments drawn is three less than the number of
	 * points passed, assuming there are a least four points in the
	 * geometry matrix (in this case geom2). This call will draw 3
	 * overlapping curve segments in the current basis matrix - still
	 * bezier.
	 */
	crvn(6L, geom2);

	qread(&val);

	/*
	 * load in the cardinal basis matrix
	 */
	defbasis((short)1, cardinal);
	curvebasis((short)1);

	color(MAGENTA);

	cmov2(70.0, 10.0);
	charstr("Cardinal Curve Segment");

	/*
	 * plot out a curve segment using the cardinal basis matrix
	 */
	crv(geom1);

	cmov2(-190.0, 400.0);
	charstr("Three overlapping Cardinal Curves");

	/*
	 * now draw a bunch of them again.
	 */
	crvn(6L, geom2);

	qread(&val);

	/*
	 * change the basis matrix again
	 */
	defbasis((short)1, bspline);
	curvebasis((short)1);

	color(GREEN);

	cmov2(70.0, -40.0);
	charstr("Bspline Curve Segment");

	/*
	 * now draw our curve segment in the new basis...
	 */
	crv(geom1);

	cmov2(-190.0, 350.0);
	charstr("Three overlapping Bspline Curves");

	/*
	 * ...and do some overlapping ones
	 */
	crvn(6L, geom2);

	qread(&val);

	gexit();
}
