
#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

Coord carray[][3] = { -1.0,  -1.0,   1.0, /* front */
		      1.0,  -1.0,   1.0,
		      1.0,   1.0,   1.0,
		     -1.0,   1.0,   1.0,
		     -1.0,  -1.0,  -1.0, /* rear */
		      1.0,  -1.0,  -1.0,
		      1.0,   1.0,  -1.0,
		     -1.0,   1.0,  -1.0
		};

int	nplanes;

/*
 * drawcube
 *
 *	draw the cube setting colours if available
 */
void
drawcube()
{
	if (nplanes > 1)
		color(RED);

	/* Front */
	pmv(carray[0][0], carray[0][1], carray[0][2]);
	pdr(carray[1][0], carray[1][1], carray[1][2]);
	pdr(carray[2][0], carray[2][1], carray[2][2]);
	pdr(carray[3][0], carray[3][1], carray[3][2]);
	pclos();
	
	if (nplanes > 1)
		color(GREEN);

	/* Back */
	pmv(carray[5][0], carray[5][1], carray[5][2]);
	pdr(carray[4][0], carray[4][1], carray[4][2]);
	pdr(carray[7][0], carray[7][1], carray[7][2]);
	pdr(carray[6][0], carray[6][1], carray[6][2]);
	pclos();

	if (nplanes > 1)
		color(YELLOW);

	/* Right side */
	pmv(carray[1][0], carray[1][1], carray[1][2]);
	pdr(carray[5][0], carray[5][1], carray[5][2]);
	pdr(carray[6][0], carray[6][1], carray[6][2]);
	pdr(carray[2][0], carray[2][1], carray[2][2]);
	pclos();

	if (nplanes > 1)
		color(BLUE);

	/* Left side */
	pmv(carray[0][0], carray[0][1], carray[0][2]);
	pdr(carray[3][0], carray[3][1], carray[3][2]);
	pdr(carray[7][0], carray[7][1], carray[7][2]);
	pdr(carray[4][0], carray[4][1], carray[4][2]);
	pclos();

	if (nplanes > 1)
		color(MAGENTA);

	/* Top */
	pmv(carray[2][0], carray[2][1], carray[2][2]);
	pdr(carray[6][0], carray[6][1], carray[6][2]);
	pdr(carray[7][0], carray[7][1], carray[7][2]);
	pdr(carray[3][0], carray[3][1], carray[3][2]);
	pclos();
	
	if (nplanes > 1)
		color(CYAN);

	/* Bottom */
	pmv(carray[0][0], carray[0][1], carray[0][2]);
	pdr(carray[4][0], carray[4][1], carray[4][2]);
	pdr(carray[5][0], carray[5][1], carray[5][2]);
	pdr(carray[1][0], carray[1][1], carray[1][2]);
	pclos();
}


main(argc, argv)
{
	float	t, dt = 0.2;
	int	r, dr = 100;
	short	val;

	prefsize(300L, 300L);

	vinit("mswin");
	winopen("cube");

	qdevice(KEYBD);
	unqdevice(INPUTCHANGE);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	nplanes = getplanes();

	color(BLACK);
	clear();

	window(-1.5, 1.5, -1.5, 1.5, 9.0, -5.0);
	lookat(0.0, 0.0, 12.0, 0.0, 0.0, 0.0, 0.0);

	backface(1);

	if (argc == 1)
		polymode(PYM_LINE);
	
	doublebuffer();
	gconfig();

	t = 0.0;

	do {
		for (r = 0; r < 3600; r += dr) {
			color(BLACK);
			clear();
			pushmatrix();
				translate(0.0, 0.0, t);
				rotate(r, 'y');
				rotate(r, 'z');
				rotate(r, 'x');
				color(WHITE);
				drawcube();
				if (nplanes == 1 && argc > 1) {
					polymode(PYM_LINE);
					color(0);
					drawcube();
				}
				if (argc > 1)
					polymode(PYM_FILL);

			popmatrix();

			t += dt;
			if (t > 3.0 || t < -18.0)
				dt = -dt;

			swapbuffers();

			if (qtest()) {
				qread(&val);
				gexit();
				exit(0);
			}
		}
	} while(1);
}

