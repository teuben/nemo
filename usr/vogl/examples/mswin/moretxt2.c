#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

/*
 * demonstrate still more features of text
 */
main(argc, argv)
	int	argc;
	char	**argv;
{
	short	val;
	
	vinit("mswin");
	winopen("moretxt2");

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

	htextsize(0.03, 0.04);

	ortho2(0.0, 1.0, 0.0, 1.0);


	color(RED);
	clear();


	drawstuff();

	/* Now do it all with the text rotated .... */

	htextang(45.0);
	drawstuff();

	htextang(160.0);
	drawstuff();

	htextang(270.0);
	drawstuff();

	/* Now with a single character */

	htextang(0.0);
	drawstuff2();

	htextang(45.0);
	drawstuff2();

	htextang(160.0);
	drawstuff2();

	htextang(270.0);
	drawstuff2();

	gexit();
}

drawstuff()
{
	short	val;

	color(BLACK);
	rectf(0.1, 0.1, 0.9, 0.9);
	color(WHITE);
	move2(0.1, 0.5);
	draw2(0.9, 0.5);
	move2(0.5, 0.1);
	draw2(0.5, 0.9);

	color(GREEN);
	move2(0.5, 0.5);
	hleftjustify(1);
	hcharstr("This is Left Justified text");

	qread(&val);

	color(BLACK);
	rectf(0.1, 0.1, 0.9, 0.9);
	color(WHITE);
	move2(0.1, 0.5);
	draw2(0.9, 0.5);
	move2(0.5, 0.1);
	draw2(0.5, 0.9);

	color(YELLOW);
	move2(0.5, 0.5);
	hcentertext(1);
	hcharstr("This is Centered text");
	hcentertext(0);

	qread(&val);

	color(BLACK);
	rectf(0.1, 0.1, 0.9, 0.9);
	color(WHITE);
	move2(0.1, 0.5);
	draw2(0.9, 0.5);
	move2(0.5, 0.1);
	draw2(0.5, 0.9);

	color(MAGENTA);
	move2(0.5, 0.5);
	hrightjustify(1);
	hcharstr("This is Right Justified text");
	hrightjustify(0);

	qread(&val);
}

drawstuff2()
{
	short	val;

	color(BLACK);
	rectf(0.1, 0.1, 0.9, 0.9);
	color(WHITE);
	move2(0.1, 0.5);
	draw2(0.9, 0.5);
	move2(0.5, 0.1);
	draw2(0.5, 0.9);

	color(GREEN);
	move2(0.5, 0.5);
	hleftjustify(1);
	hdrawchar('B');

	qread(&val);

	color(BLACK);
	rectf(0.1, 0.1, 0.9, 0.9);
	color(WHITE);
	move2(0.1, 0.5);
	draw2(0.9, 0.5);
	move2(0.5, 0.1);
	draw2(0.5, 0.9);

	color(YELLOW);
	move2(0.5, 0.5);
	hcentertext(1);
	hdrawchar('B');
	hcentertext(0);

	qread(&val);

	color(BLACK);
	rectf(0.1, 0.1, 0.9, 0.9);
	color(WHITE);
	move2(0.1, 0.5);
	draw2(0.9, 0.5);
	move2(0.5, 0.1);
	draw2(0.5, 0.9);

	color(MAGENTA);
	move2(0.5, 0.5);
	hrightjustify(1);
	hdrawchar('B');
	hrightjustify(0);

	qread(&val);
}
