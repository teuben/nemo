/*
 * Demonstrate a rotating translating tetrahedron.
 */

#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

#ifndef TC
#include <math.h>
#else
extern double sin(), cos();
#endif

#define	TETRAHEDRON	1L
#define	NSIDES	3
#define	NFACES	4
#define	NPNTS	4

Coord	points[NPNTS][3] = {
	{-0.5, 0.866, -0.667},
	{-0.5, -0.866, -0.667},
	{ 1.0, 0.0, -0.667},
	{ 0.0, 0.0, 1.334}
};

int	faces[NFACES][NSIDES] = {
	{2, 1, 0},
	{0, 1, 3},
	{1, 2, 3},
	{2, 0, 3}
};

int	colface[NFACES] = {
		GREEN,
		YELLOW,
		CYAN,
		MAGENTA
};

main(argc, argv)
	int	argc;
	char	**argv;
{
	char	dev[20];
	int	i, but;
	int	rotval = 0, drotval = 2;
	float	R = 1.6, tx = 0.0, tz = R, zeye = 5.0;
	int	do_backface = 0;
	int	do_fill = 0;
	short	val;

	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-b"))
			do_backface = 1;
		if (!strcmp(argv[i], "-f"))
			do_fill = 1;
	}

	prefsize(400L, 400L);

	winopen("tetra");          /* set up device */

	qdevice(ESCKEY);
	qdevice(QKEY);
	unqdevice(INPUTCHANGE);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	doublebuffer();
	gconfig();


	polymode(PYM_LINE);
	if (do_fill)
		polymode(PYM_FILL);

	if (do_backface)
		backface(1);

	/*
	 * set up a perspective projection with a field of view of
	 * 40.0 degrees, aspect ratio of 1.0, near clipping plane 0.1,
	 * and the far clipping plane at 1000.0.
	 */
	perspective(400, 1.0, 0.001, 15.0);
	lookat(0.0, 0.0, zeye, 0.0, 0.0, 0.0, 0);


	/*
	 * Make a tetrahedron object
	 */

	maketetra();

	do {
		for (rotval = 0; rotval < 3600; rotval += drotval) {
			color(BLACK);
			clear();

			/*
			 * Rotate the whole scene...(this acumulates - hence
			 * drotval)
			 */
			rotate(drotval, 'x');
			rotate(drotval, 'z');

			color(RED);
			pushmatrix();
				rotate(900, 'x');
				circ(0.0, 0.0, R);
			popmatrix();

			color(BLUE);
			move(0.0, 0.0, 0.0);
			draw(tx, 0.0, tz);
			
			/*
			 * Remember! The order of the transformations is
			 * the reverse of what is specified here in between
			 * the pushmatrix and the popmatrix. These ones don't
			 * accumulate because of the push and pop.
			 */
			pushmatrix();
				translate(tx, 0.0, tz);
				rotate(rotval, 'x');
				rotate(rotval, 'y');
				rotate(rotval, 'z');
				scale(0.4, 0.4, 0.4);
				callobj(TETRAHEDRON);
			popmatrix();

			tz = R * cos((double)(rotval * 3.1415926535 / 180));
			tx = R * sin((double)(rotval * 3.1415926535 / 180));

			swapbuffers();

			if (qtest()) {
/*
				but = (int)qread(&val);
				fprintf(stderr, "but = %c (%d)\n", but, but);
*/
				gexit();
				exit(0);
			}
		}

	} while (1);
}

/*
 * maketetra
 *
 *	draw a tetrahedron
 */
maketetra()
{
	int	i, j;

	makeobj(TETRAHEDRON);

	for (i = 0; i < NFACES; i++) {
		color(colface[i]);
		bgnpolygon();
			for (j = 0; j < NSIDES; j++) 
				v3f(points[faces[i][j]]);
		endpolygon();
	}

	closeobj();
}
