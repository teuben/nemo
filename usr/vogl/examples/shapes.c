#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

/*
 * This program shows some of the simple primitives.
 */
main()
{
	Screencoord	minx, maxx, miny, maxy, edgelength;
	short		val;

	int     ls, lw;

	lw = 1;
	ls = 0xff00;

	winopen("shapes");
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	deflinestyle(34, ls);
	setlinestyle(34);
	linewidth(lw);

	/*
	 * the two lines below clear the screen to white if we have
	 * colours, on a monochrome device color is ignored so the
	 * screen will be cleared to its background color, normally black.
	 */
	color(BLACK);
	clear();

	/*
	 * set the screen to be 2.0 units wide and 2.0 units wide, with
	 * the drawable coordinates going from -1.0 to 1.0.
	 */
	ortho2(-1.0, 1.0, -1.0, 1.0);

	color(MAGENTA);

	/*
	 * okay, so we want to draw in the range -1 to 1, but we
	 * only want to draw in the top lefthand corner of the
	 * screen. The call to viewport allows us to do this. As
	 * viewport always takes screen coordinates, we need to
	 * call getviewport to found out how big our screen is
	 * at the moment. We use the values returned from getviewport
	 * calculate the positions for our new viewport. We note
	 * that on an Iris (0,0) is the bottom left pixel.
	 */
	getviewport(&minx, &maxx, &miny, &maxy);

	viewport(minx, (maxx - minx) / 2, (maxy - miny) / 2, maxy);

	cmov2(-0.9, -0.5);		/* write out a heading */
	charstr("rect");

	/*
	 * draw a rectangle around the points (-0.2, -0.2), (-0.2, 0.2),
	 * (0.3, 0.2), and (0.3, -0.2).
	 */
	rect(-0.2, -0.2, 0.3, 0.2);

	color(BLUE);

	/*
	 * now we want to draw in the top right corner of the screen,
	 * and we want to draw a circular circle so we must make sure
	 * our viewport is square (if it isn't we'll get an ellipse),
	 * so we calculate the shortest edge and set up a viewport which
	 * is in the top right region of the screen, but doesn't necessarilly
	 * occupy all the top right corner.
	 */
						/* find smallest edge */
	if (maxx - minx > maxy - miny)
		edgelength = (maxy - miny) / 2;
	else 
		edgelength = (maxx - minx) / 2;

						/* create a square viewport */

	viewport((maxx - minx) / 2, (maxx - minx) / 2 + edgelength, (maxy - miny) / 2, (maxy - miny) / 2 + edgelength);

	cmov2(-0.9, -0.5);
	charstr("circle");

	/*
	 * draw a circle of radius 0.4 around the point (0.0, 0.0)
	 */
	circ(0.0, 0.0, 0.4);

	color(GREEN);

	/*
	 * bottom left hand corner.
	 */
	viewport(minx, (maxx - minx) / 2, miny, (maxy - miny) / 2);

	cmov2(-0.9, -0.5);
	charstr("ellipse");

	/*
	 * To draw an ellipse we change the aspect ratio so it is no longer
	 * 1 and call circ. In this case we use ortho2 to make the square
	 * viewport appear to be higher than it is wide. Alternatively you
	 * could use arc to construct one.
	 *
	 * The call to pushmatrix saves the current viewing transformation.
	 * After the ortho2 has been done, we restore the current viewing
	 * transformation with a call to popmatrix. (Otherwise everything
	 * after the call to ortho would come out looking squashed as the
	 * world aspect ratio is no longer 1).
	 */
	pushmatrix();
		ortho2(-1.0, 1.0, -1.0, 2.0);
		circ(0.0, 0.5, 0.4);
	popmatrix();

	color(RED);

	/*
	 * bottom right hand corner
	 */
	viewport((maxx - minx) / 2, maxx, miny, (maxy - miny) / 2);

	cmov2(-0.9, -0.5);
	charstr("arc");

	/*
	 * draw an arc centered at (0.0, 0.0), radius of 0.4. 0.0 is the start
	 * angle and 90.0 is the end angle of the arc being drawn. So this
	 * draws a quarter circle - unless our viewport isn't square.
	 */
	arc(0.0, 0.0, 0.4, 0, 900);

	/*
	 * enable the keyboard
	 */
	qdevice(KEYBD);
	unqdevice(INPUTCHANGE);

	/*
	 * wait for the event
	 */
	qread(&val);

	gexit();
}
