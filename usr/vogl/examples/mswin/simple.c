#include <stdio.h>



#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

/*
 * A program showing basic line drawing, hardware text and (if applicable)
 * colour. We set the coordinate system to -1.0 to 1.0 in X and Y.
 */
main(ac, av)
	int	ac;
	char	**av;
{
	char	*p, tmp[2];
	float	cw, ch;
	short	val;

	prefposition(100L, 700L, 100L, 500L);
	vinit("mswin");
	winopen("simple");
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	if (ac == 2)
		font(atoi(av[1]));	/* change font to the argument */

	hfont("times.rb");
	htextsize(0.1, 0.2);

	color(BLACK);		/* set current color */
	clear();		/* clear screen to current color */

	ortho2(-1.0, 1.0, -1.0, 1.0);	/* set bounds for drawing */

	color(GREEN);
			/* 2 d move to start where we want drawstr to start */
	move2(-0.9, 0.9);

	hcharstr("A Simple Example");	/* draw string in current color */

	/*
	 * the next four lines draw the x 
	 */
	move2(0.0, 0.0);
	draw2(0.76, 0.76);
	move2(0.0, 0.76);
	draw2(0.76, 0.0);

	move2(0.0, 0.5);
	hcharstr("x done");
	hcharstr("Hi there");

	move(0.0, 0.1, -1.0);
	/*
	 * One character at a time...
	 */
	tmp[1] = '\0';
	for (p = "hello world"; *p != (char)NULL; p++) {
		tmp[0] = *p;
		hcharstr(tmp);     
	}

	/*
	 * the next five lines draw the square
	 */
	deflinestyle(34, 0xff00);
	setlinestyle(34);
	linewidth(3);

	move2(0.,0.);
	draw2(.76,0.);
	draw2(.76,.76);
	draw2(0.,.76);
	draw2(0.,0.);

	qdevice(KEYBD);		/* enable the keyboard */
	unqdevice(INPUTCHANGE);
	unqdevice(REDRAW);

	qread(&val);		/* wait for some input */

	gexit();		/* set the screen back to its original state */
}
