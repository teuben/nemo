#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

/*
 * drawgrid
 *
 *	draw a grid in the middle of the screen
 */
void
drawgrid()
{
	float	x;
	int	i;

	color(GREEN);

	rect(0.1, 0.4, 0.9, 0.6);

	x = 0.2;
	for (i = 0; i < 8; i++) {
		move2(x, 0.4);
		draw2(x, 0.6);
		x += 0.1;
	}
	move2(0.1, 0.5);
	draw2(0.9, 0.5);

	color(YELLOW);
}

/*
 * demonstrate some more features of text
 */
main(argc, argv)
	int	argc;
	char	**argv;
{
	int	i;
	float	x;
	short	val;
	
	winopen("moretxt");

	unqdevice(INPUTCHANGE);
	qdevice(KEYBD);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	if (argc == 2)
		hfont(argv[1]);
	else
		hfont("futura.l");

	htextsize(0.05, 0.05);

	ortho2(0.0, 1.0, 0.0, 1.0);

	color(BLACK);
	clear();

	drawgrid();

	/*
	 * show some scaled text on the grid (In the bottom part)
	 */
	hboxtext(0.1, 0.4, 0.8, 0.1, "{This is Some text] | $");

	qread(&val);

	color(BLACK);
	clear();

	drawgrid();

	/*
	 * centertext causes text to be centered around the current graphics
	 * position this is especially usefull if you want your text to come
	 * out centered on a line, or a character to be centered on a point
	 * in a graph. A non-zero argument turns centertext on.
	 *
	 * show a string centered on the center line
	 */
	hcentertext(1);

	hboxtext(0.5, 0.5, 0.8, 0.1, "{This is Some Centered text] | $");

	/*
	 * turn centertext off. We use an argument with the value zero.
	 */
	hcentertext(0);

	qread(&val);

	color(BLACK);
	clear();

	/*
	 * rotate the grid so that it is the same angle as the text after
	 * textang for text ang.
	 */
	pushmatrix();
		translate(0.5, 0.5, 0.0);
		rotate(900, 'z');
		translate(-0.5, -0.5, 0.0);

		drawgrid();
	popmatrix();

	/*
	 * turn on centered text again
	 */
	hcentertext(1);

	/*
	 * set the angle to 90.
	 */
	htextang(90.0);

	/*
	 * draw the string
	 */
	hboxtext(0.5, 0.5, 0.8, 0.1, "{This is Some Rotated Centered text] | $");

	/*
	 * turn off center text
	 */
	hcentertext(0);

	/*
	 * set text angle back to 90
	 */
	htextang(0.0);

	qread(&val);

	color(BLACK);
	clear();

	drawgrid();

	/*
	 * as all the software fonts are proportionally spaced we use
	 * the fixedwidth call to make each character take the same amount
	 * of horizontal space. As with centertext this is done by passing
	 * fixedwidth a non-zero argument.
	 */
	hfixedwidth(1);

	hboxtext(0.1, 0.5, 0.8, 0.1, "{This is Some Fixedwidth text] | $");

	qread(&val);

	color(BLACK);
	clear();

	drawgrid();

	/*
	 * now try centered and fixewidth at the same time
	 */
	hcentertext(1);

	color(RED);
	move2(0.5, 0.5);
	hcharstr("{This is Some Cent.Fixedwidth text] | $");

	hcentertext(0);
	
	qread(&val);
	color(BLACK);
	clear();

	drawgrid();

	color(RED);
	move2(0.9, 0.4);
	hrightjustify(1);
	hfixedwidth(0);
	hcharstr("{This is Some Right Justified text] | $");

	qread(&val);
	color(BLACK);
	clear();

	drawgrid();

	/*
	 * scale the text so tha a character is the size of a box in
	 * the grid.
	 */
	hboxfit(0.8, 0.1, 8);

	hfixedwidth(1);
	/*
	 * draw the two strings fixedwidth (it is still turned on)
	 */
	hleftjustify(1);
	move2(0.1, 0.4);
	hcharstr("ABCDefgh");

	move2(0.1, 0.5);
	hcharstr("IJKLmnop");

	qread(&val);

	gexit();
}
