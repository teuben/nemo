
#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

#define	CUBE_SIZE	200.0
#define	TRANS		25.0
#define	SCAL		0.1

main()
{
        char    *p;
	float	tdir = TRANS;
	float	scal = 1.0 + SCAL;
	int	but, nplanes;
	int	x, y, i, n;
	short	val;
	int	bf = 1;
	int	fill = 1;

	prefsize(500L, 500L);

	vinit("mswin");
	winopen("lcube");

	unqdevice(INPUTCHANGE);
	qdevice(SKEY);
	qdevice(XKEY);
	qdevice(YKEY);
	qdevice(ZKEY);
	qdevice(EQUALKEY);
	qdevice(MINUSKEY);
	qdevice(ESCKEY);
	qdevice(QKEY);
	qdevice(FKEY);
	qdevice(BKEY);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;
	

	window(-800.0, 800.0, -800.0, 800.0, -800.0, 800.0);
	lookat(0.0, 0.0, 1500.0, 0.0, 0.0, 0.0, 0);

	if ((nplanes = getplanes()) == 1)
		makecubes(0);

	makecubes(1);

	backface(1);
		
	doublebuffer();
	gconfig();

	/*
	 * Doublebuffer does a backbuffer(TRUE)....
	 */

	while(1) {
		x = 500 - (int)getvaluator(MOUSEX);
		y = 500 - (int)getvaluator(MOUSEY);
		x *= 3;
		y *= 3;
		pushmatrix();
			rotate(x, 'y');
			rotate(y, 'x');
			color(BLACK);
			clear();
			callobj((Object)3);
			if (nplanes == 1)
				callobj((Object)2);
		popmatrix();
		swapbuffers();

		if (qtest()) {
			but = qread(&val);
			but = qread(&val);	/* swallow up event */

			switch (but) {

			case SKEY:
				scale(scal, scal, scal);
				break;
			case XKEY:
				translate(tdir, 0.0, 0.0);
				break;
			case YKEY:
				translate(0.0, tdir, 0.0);
				break;
			case ZKEY:
				translate(0.0, 0.0, tdir);
				break;
			case MINUSKEY:
				tdir = -tdir;

				if (scal < 1.0)
					scal = 1.0 + SCAL;
				else
					scal = 1.0 - SCAL;

				break;
			case EQUALKEY:
				tdir = TRANS;
				break;
			case BKEY:
				bf = !bf;
				backface(bf);
				break;
			case FKEY:
				fill = !fill;
				if (fill)
					polymode(PYM_FILL);
				else
					polymode(PYM_LINE);
				break;
			case ESCKEY:
			case QKEY:
				gexit();
				exit(0);
			default:
				;
			}
		}
	}
}

makecubes(fill)
	int	fill;
{
	makeobj((Object)(fill + 2));
		if (!fill)
			color(BLACK);

		pushmatrix();
			translate(0.0, 0.0, CUBE_SIZE);
			if (fill) {
				color(RED);
				rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
			} else
				rect(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(CUBE_SIZE, 0.0, 0.0);
			rotate(900, 'y');
			if (fill) {
				color(GREEN);
				rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
			} else
				rect(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(0.0, 0.0, -CUBE_SIZE);
			rotate(1800, 'y');
			if (fill) {
				color(BLUE);
				rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
			} else 
				rect(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(-CUBE_SIZE, 0.0, 0.0);
			rotate(-900, 'y');
			if (fill) {
				color(CYAN);
				rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
			} else
				rect(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(0.0, CUBE_SIZE, 0.0);
			rotate(-900, 'x');
			if (fill) {
				color(MAGENTA);
				rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
			} else
				rect(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(0.0, -CUBE_SIZE, 0.0);
			rotate(900, 'x');
			if (fill) {
				color(YELLOW);
				rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
			} else
				rect(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

	closeobj();
}
