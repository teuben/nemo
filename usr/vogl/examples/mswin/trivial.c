#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

/*
 * the basic test program for a driver if we can draw a line and do
 * hardware text we are almost there!
 */
main(argc, argv)
	int	argc;
	char	**argv;
{
	short	val;

	vinit("mswin");
	winopen("trivial");

	if (argc == 2)
		font(atoi(argv[1]));		/* set font to argument */

	qdevice(KEYBD);			/* enable the keyboard */
	unqdevice(INPUTCHANGE);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	color(BLACK);			/* we want to clear in black */
	clear();			/* clear to current color */

	ortho2(-1.0, 1.0, -1.0, 1.0);	/* set up the coordinate system */

	color(GREEN);			/* set current color to green */

	move2(-1.0, 0.0);		/* draw a horizontal line at y = 0 */
	draw2(1.0, 0.0);

	qread(&val);			/* pause for some input */

	move2(0.0, 0.0);		/* draw a line along x = 0 */
	draw2(0.0, 1.0);

	cmov2(0.0, 0.0);		/* move to the middle of the screen */
	charstr("Hello");		/* draw "Hello" starting at the origin */
	qread(&val);			/* pause again */

	gexit();			/* set screen back to original state */
}
