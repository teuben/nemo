#include <stdio.h>

#ifdef SGI
#include <gl.h>
#include <device.h>
#else
#include "vogl.h"
#include "vodevice.h"
#endif

extern int	getcharacter();

#define	XCOORD(x)	((int)(x) - (int)'R')
#define	YCOORD(y)	((int)'R' - (int)(y))	/* invert as in tv coords */

/*
 * newpage
 *
 *	draw up a new page with title, boxes, etc..
 */
newpage(fname, pageno)
	char	*fname;
	int	pageno;
{
	char	str[100];

	hcentertext(0);

	htextsize(0.05, 0.07);

	move2(-0.91, 0.9);
	hcharstr("Hershey Character File: ");
	hcharstr(fname);

	move2(0.45, 0.9);
	sprintf(str, "Page No: %d", pageno);
	hcharstr(str);

	htextsize(0.03, 0.03);

	hcentertext(1);
}

/*
 * display the hershey data set in the input file.
 */
main(ac, av)
	int	ac;
	char	**av;
{
	FILE	*fp;
	int	charno, numpairs, page;
	char    c, device[20], buf[1000], *p, str[100];
	float	x, y, ox, oy;
	short	val;

	if (ac < 2) {
		fprintf(stderr, "hdisp: usage hdisp datafile\n");
		exit(1);
	}

	if ((fp = fopen(av[1], "r")) == NULL) {
		fprintf(stderr, "hdisp: unable to open file %s\n", av[1]);
		exit(1);
	}

	hfont("times.r");

	winopen("hdisp");
	qdevice(KEYBD);
	ortho2(-1.0, 1.0, -1.0, 1.0);

	color(BLACK);
	clear();

	ox = -0.8;
	oy = 0.75;

	page = 1;

	color(WHITE);

	newpage(av[1], page);

	while (getcharacter(fp, &charno, &numpairs, buf)) {

		if (buf[2] != 0) {
			p = &buf[2];		/* skip the width bytes */

			x = XCOORD(*p++) / 280.0;
			y = YCOORD(*p++) / 280.0;
			move2((Coord)(ox + x), (Coord)(oy + y));

			while (*p != 0) {
				if (*p == ' ') {
					p += 2;
					x = XCOORD(*p++) / 280.0;
					y = YCOORD(*p++) / 280.0;
					move2((Coord)(ox + x), (Coord)(oy + y));
				} else {
					x = XCOORD(*p++) / 280.0;
					y = YCOORD(*p++) / 280.0;
					draw2((Coord)(ox + x), (Coord)(oy + y));
				}
			}
		}

		move2((Coord)ox, (Coord)(oy - 0.11));
		sprintf(str, "(%d)", charno);
		hcharstr(str);

		ox += 0.22;
		if (ox > 0.9) {
			oy -= 0.22;
			ox = -0.8;
		}
		if (oy < -0.9) {
			oy = 0.75;
			qread(&val);

			color(BLACK);
			clear();

			color(WHITE);

			newpage(av[1], ++page);
		}
	}

	qread(&val);

	gexit();
}
