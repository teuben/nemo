#include <stdio.h>
#ifdef SGI
#include <gl.h>
#include <device.h>
#else
#include "vogl.h"
#include "vodevice.h"
#endif

/*
 *	displays every character in a hershey font at 64 characters
 * per screen. Note: this program reads the binary format as created
 * by h2v.
 */
main(ac, av)
	int	ac;
	char	**av;
{
	char	dev[50];
	int	i, nchars;
	float	x, y;
	short	val;

	if (ac != 2) {
		fprintf(stderr, "fdisp: usage fdisp fontname\n");
		exit(1);
	}

	winopen("fdisp");
	ortho2(-1.0, 1.0, -1.0, 1.0);
	qdevice(KEYBD);
	color(BLACK);
	clear();

	color(GREEN);

	hfont(av[1]);

	nchars = hnumchars();

	htextsize(0.2, 0.2);

	x = -0.94;
	y = 0.77;
	for (i = 0; i < nchars; i++) {
		move2(x, y);
		hdrawchar(' ' + i);
		x += 0.25;
		if (x > 0.86) {
			y -= 0.25;
			if (y < -1.1) {
				qread(&val);
				color(BLACK);
				clear();
				color(GREEN);
				y = 0.77;
			}
			x = -0.94;
		}
	}

	qread(&val);

	gexit();
}

