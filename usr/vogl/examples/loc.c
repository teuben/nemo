#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

/*
 * a routine to demonstrate using locator.
 */
main()
{
	int		i, bt, act, nchars;
	short		data;
	Scoord		x, y, sx, sy;
	Screencoord	minx, maxx, miny, maxy;

	ginit();

	color(BLACK);
	clear();

	color(BLUE);

	getviewport(&minx, &maxx, &miny, &maxy);

	ortho2((Coord)minx, (Coord)maxx, (Coord)miny, (Coord)maxy);

	/*
	 * draw some axes
	 */
	move2s((Scoord)minx, (Scoord)((maxy - miny) / 2));
	draw2s((Scoord)maxx, (Scoord)((maxy - miny) / 2));

	move2s((Scoord)((maxx - minx) / 2), (Scoord)miny);
	draw2s((Scoord)((maxx - minx) / 2), (Scoord)maxy);

	color(GREEN);

	/*
	 * enable the left and middle mouse buttons
	 */
	unqdevice(INPUTCHANGE);
	qdevice(LEFTMOUSE);
	qdevice(MIDDLEMOUSE);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&data) != REDRAW)
		;

	act = 0;

	/*
	 * getvaluator tells us the valuator's value. In
	 * this case it's the X and Y positions of the mouse.
	 * Note: these come back to us in screen coordinates.
	 */
	while((bt = qread(&data)) != MIDDLEMOUSE) {
		sx = getvaluator(MOUSEX);
		sy = getvaluator(MOUSEY);
		if (bt == -1) {
			gexit();
			printf("No locator device found\n");
			exit(0);
		} else {
			if (act) {
				act = 0;
				move2s(sx, sy);
				draw2s(x, y);
			} else {
				act = 1;
				x = sx;
				y = sy;
			}
		}
		(void)qread(&data);	/* swallow the up event */
	}

	gexit();

}
