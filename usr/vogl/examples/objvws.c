#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

#define		CUBE		1L
#define		TOPLEFT		2L
#define		TOPRIGHT	3L
#define		BOTTOMLEFT	4L
#define		BOTTOMRIGHT	5L

/*
 * Demonstrate just how much you can put in an object
 */
main()
{
	void		makecube();

	Screencoord	minx, maxx, miny, maxy;
	short		val;

	winopen("objvws");
	unqdevice(INPUTCHANGE);
	qdevice(KEYBD);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	getviewport(&minx, &maxx, &miny, &maxy);


	pushviewport();

	hfont("futura.m");
	htextsize(0.5, 0.9);

	color(BLACK);
	clear();

	makecube();

	/*
	 * set up an object which draws in the top left of the screen.
	 */
	makeobj(TOPLEFT);
		viewport(minx, (maxx - minx) / 2, (maxy - miny) / 2, maxy);
		ortho2(-5.0, 5.0, -5.0, 5.0);

		color(RED);

		rect(-5.0, -5.0, 5.0, 5.0);

		perspective(400, 1.0, 0.1, 1000.0);
		lookat(5.0, 8.0, 5.0, 0.0, 0.0, 0.0, 0);

		callobj(CUBE);

		color(GREEN);

		move2(-4.5, -4.5);
		hcharstr("perspective/lookat");
	closeobj();

	/*
	 * now set up one which draws in the top right of the screen
	 */
	makeobj(TOPRIGHT);
		viewport((maxx - minx) / 2, maxx, (maxy - miny) / 2, maxy);
		ortho2(-5.0, 5.0, -5.0, 5.0);

		color(GREEN);

		rect(-5.0, -5.0, 5.0, 5.0);

		window(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0);
		lookat(5.0, 8.0, 5.0, 0.0, 0.0, 0.0, 0);

		callobj(CUBE);

		color(RED);

		move2(-4.5, -4.5);
		hcharstr("window/lookat");
	closeobj();

	/*
	 * try the bottom left
	 */
	makeobj(BOTTOMLEFT);
		viewport(minx, (maxx - minx) / 2, miny, (maxy - miny) / 2);
		ortho2(-5.0, 5.0, -5.0, 5.0);

		color(MAGENTA);

		rect(-5.0, -5.0, 5.0, 5.0);

		perspective(400, 1.0, 0.1, 1000.0);
		polarview(15.0, 300, 300, 300);

		callobj(CUBE);

		color(YELLOW);

		move2(-4.5, -4.5);
		hcharstr("perspective/polarview");
	closeobj();

	/*
	 * and the bottom right
	 */
	makeobj(BOTTOMRIGHT);
		viewport((maxx - minx) / 2, maxx, miny, (maxy - miny) / 2);
		ortho2(-5.0, 5.0, -5.0, 5.0);

		color(CYAN);

		rect(-5.0, -5.0, 5.0, 5.0);

		window(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0);
		polarview(8.0, -1800, -300, 1800);

		callobj(CUBE);

		color(BLUE);

		move2(-4.5, -4.5);
		hcharstr("window/polarview");
	closeobj();

	/*
	 * now draw them
	 */
	callobj(TOPLEFT);
	callobj(TOPRIGHT);
	callobj(BOTTOMLEFT);
	callobj(BOTTOMRIGHT);

	qread(&val);

	gexit();
}

/*
 * makecube
 *
 *	set up a cube
 */
void
makecube()
{
	void	side();

	makeobj(CUBE);

		/*
		 * The border around the cube
		 */
		rect(-5.0, -5.0, 10.0, 10.0);

		/*
		 * Make the cube from 4 squares
		 */
		pushmatrix();
			side();
			rotate(900, 'x');
			side();
			rotate(900, 'x');
			side();
			rotate(900, 'x');
			side();
		popmatrix();

	closeobj();
}

/*
 * side
 *
 *	define a face for the cube
 */
void
side()
{
	pushmatrix();
		translate(0.0, 0.0, 1.0);
		rect(-1.0, -1.0, 1.0, 1.0);
	popmatrix();
}

