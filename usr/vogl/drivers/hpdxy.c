#undef VOGLE

/*
 *	DXY & HPGL Driver for vogle/vogl.
 */
#include <stdio.h>
#ifdef VOGLE
#include "vogle.h"
#else
#include "vogl.h"
#endif

extern FILE	*_voutfile();

static int	plotlstx, plotlsty;	/* position of last draw */

static FILE	*fp;

#define	P_RESET		0
#define	P_MOVE		1
#define	P_DRAW		2
#define	P_TXTSIZE	3
#define	P_BEGTXT	4
#define	P_ENDTXT	5
#define	P_PEN		6

/*
 * basic commands for hpgl
 */
/* 
 * Changed to the delimiters to commas and removed spaces (for older plotter).
 * (mike@penguin.gatech.edu) 
 */
static char	*hpgl[] = {
	"DF;\n",
	"PU%d,%d;\n",
	"PD%d,%d;\n",
	"SI%.4f,%.4f;\n",
	"LB",
	"\003\n",
	"SP%d;\n"
};

/*
 * basic commands for dxy
 */
static char	*dxy[] = {
	"",
	"M %d,%d\n",
	"D %d,%d\n",
	"S %d\n",
	"P",
	"\n",
	"J %d\n"
};

static char	**plotcmds;

/*
 * noop
 *
 *      do nothing but return-1
 */
static int
noop()
{
	return(-1);
}
/*
 * HPGL_common_init()
 *
 * Performs the common parts of HPGL initialization.
 */
static int
HPGL_common_init(minx, maxx, miny, maxy)
	int	minx, maxx, miny, maxy;
{
	vdevice.depth = 4;

	fp = _voutfile();

	/*
	 * The next line is for serial lines if you need to set modes
	 */
	fprintf(fp, "\033.(;\033.I81;;17;\033.N;19:IN;");

	/*
	 * Cause scaling to be 0 to maxX maxY.
	 */
	fprintf(fp, "IP%d,%d,%d,%d;", minx, miny, maxx, maxy);
	fprintf(fp, "SC0,%d,0,%d;", vdevice.sizeX, vdevice.sizeY);

	plotcmds = hpgl;
	plotlstx = -1111111;
	plotlsty = -1111111;

	return(1);
}

/*
 * HPGL_init
 *
 *	set up hp plotter. Returns 1 on success.
 */

static int
HPGL_A4_init()
{
	/*
	 * A4 paper
	 */
	vdevice.sizeX = vdevice.sizeY = 7320;

	vdevice.sizeSx = 10200;
	vdevice.sizeSy = 7320;

	/* 
	 * Changed to 7000 (from 7721) as noted by Michael J. Gourlay 
	 * (mike@penguin.gatech.edu) 
	 */
	return(HPGL_common_init(-7000, 7000, -7000, 7000));
}

static int
HPGL_A3_init()
{
	/*
	 * A3 paper
	 */
	vdevice.sizeX = vdevice.sizeY = 10560;

	vdevice.sizeSx = 14720; 
	vdevice.sizeSy = 10560; 

	return(HPGL_common_init(-10000, 10000, -10000, 10000));
}

static int
HPGL_A2_init()
{
	/*
	 * A2 paper
	 */
	vdevice.sizeX = vdevice.sizeY = 13440;

	vdevice.sizeSx = 18734;
	vdevice.sizeSy = 13440;

	return(HPGL_common_init(-13000, 13000, -13000, 13000));
}

static int
HPGL_A1_init()
{
	/*
	 * A1 paper
	 */
	vdevice.sizeX = vdevice.sizeY = 21360;

	vdevice.sizeSx = 29774;
	vdevice.sizeSy = 21360;

	return(HPGL_common_init(-21000, 21000, -21000, 21000));
}


/*
 * DXY_init
 *
 *	set up dxy plotter. Returns 1 on success.
 */
static int
DXY_init()
{
	fp = _voutfile();

	vdevice.sizeX = vdevice.sizeY = 1920; 

	vdevice.sizeSx = 2668; 
	vdevice.sizeSy = 1920; 

	plotcmds = dxy;
	plotlstx = -1;
	plotlsty = -1;

	fprintf(fp, plotcmds[P_RESET]);

	return(1);
}

/*
 * PLOT_draw
 *
 *	print the commands to draw a line from the current graphics position
 * to (x, y).
 */
static
PLOT_draw(x, y)
	int	x, y;
{
	if (plotlstx != vdevice.cpVx || plotlsty != vdevice.cpVy)
		fprintf(fp, plotcmds[P_MOVE], vdevice.cpVx, vdevice.cpVy);

	fprintf(fp, plotcmds[P_DRAW], x, y);
	plotlstx = x;
	plotlsty = y;
}

/*
 * PLOT_exit
 *
 *	exit from vogle printing the command to put away the pen and flush
 * the buffer.
 */
static
PLOT_exit()
{
	fprintf(fp, plotcmds[P_PEN], 0);
	fprintf(fp, "\033.)");
	fflush(fp);

	if (fp != stdout)
		fclose(fp);
}

/*
 * PLOT_color
 *
 *	change the current pen number.
 */
static
PLOT_color(i)
	int	i;
{ 
	fprintf(fp, plotcmds[P_PEN], i);
}

/*
 * HPGL_font
 *
 *	load in large or small
 */
static int
HPGL_font(font)
	char	*font;
{
	if (strcmp(font, "small") == 0) {
		vdevice.hwidth = 97.01;	/* Size in plotter resolution units */
		vdevice.hheight = vdevice.hwidth * 2.0;
		fprintf(fp, plotcmds[P_TXTSIZE], 0.16, 0.32);
	} else if (strcmp(font, "large") == 0) {
		vdevice.hwidth = 145.5;
		vdevice.hheight = vdevice.hwidth * 2.0;
		fprintf(fp, plotcmds[P_TXTSIZE], 0.24, 0.48);
	} else 
		return(0);

	return(1);
}

/*
 * DXY_font
 *
 *	load in large or small.
 */
static int
DXY_font(font)
	char	*font;
{
	if (strcmp(font, "small") == 0) {
		vdevice.hwidth = 24.25;
		vdevice.hheight = vdevice.hwidth * 2.0;
		fprintf(fp, plotcmds[P_TXTSIZE], 3);
	} else if (strcmp(font, "large") == 0) {
		vdevice.hwidth = 36.375;
		vdevice.hheight = vdevice.hwidth * 2.0;
		fprintf(fp, plotcmds[P_TXTSIZE], 5);
	} else 
		return(0);

	return(1);
}

/*
 * PLOT_char
 *
 *	draw a character.
 */
static
PLOT_char(c)
	char	c;
{
	int	dy, dx;

	if (plotlstx != vdevice.cpVx || plotlsty != vdevice.cpVy)
		fprintf(fp, plotcmds[P_MOVE], vdevice.cpVx, vdevice.cpVy);

	fprintf(fp, plotcmds[P_BEGTXT]);

	fprintf(fp, "%c", c);

	fprintf(fp, plotcmds[P_ENDTXT]);

	plotlstx = plotlsty = -1111111;
}

/*
 * PLOT_string
 *
 *	output a string.
 */
static
PLOT_string(s)
	char	*s;
{
	int		dy, dx;

	if (plotlstx != vdevice.cpVx || plotlsty != vdevice.cpVy)
		fprintf(fp, plotcmds[P_MOVE], vdevice.cpVx, vdevice.cpVy);

	fprintf(fp, plotcmds[P_BEGTXT]);

	fputs(s, fp);

	fprintf(fp, plotcmds[P_ENDTXT]);

	plotlstx = plotlsty = -1111111;
}

/*
 * PLOT_fill
 *
 *      "fill" a polygon
 */
static
PLOT_fill(n, x, y)
	int     n, x[], y[];
{
	int     i;

	if (plotlstx != x[0] || plotlsty != y[0])
		fprintf(fp, plotcmds[P_MOVE], x[0], y[0]);

	for (i = 1; i < n; i++)
		fprintf(fp, plotcmds[P_DRAW], x[i], y[i]);

	fprintf(fp, plotcmds[P_DRAW], x[0], y[0]);

	plotlstx = vdevice.cpVx = x[n - 1];
	plotlsty = vdevice.cpVy = y[n - 1];
}

static DevEntry hpgldev = {
	"hpgl",
	"large",
	"small",
	noop,
	PLOT_char,
	noop,
	noop,
	PLOT_color,
	PLOT_draw,
	PLOT_exit,
	PLOT_fill,
	HPGL_font,
	noop,
	noop,
	HPGL_A2_init,
	noop,
	noop,
#ifndef VOGLE
	noop,
	noop,
#endif
	PLOT_string,
	noop,
	noop
};

/*
 * _HPGL_devcpy
 *
 *	copy the HPGL device into vdevice.dev.
 */
_HPGL_A2_devcpy()
{
/* if you don't have structure assignment ...
	char    *dev, *tdev, *edev;

	dev = (char *)&hpgldev;
	tdev = (char *)&vdevice.dev;
	edev = dev + sizeof(Device);

	while (dev != edev)
		*tdev++ = *dev++;
*/
	vdevice.dev = hpgldev;
}

_HPGL_A3_devcpy()
{
	vdevice.dev = hpgldev;
	vdevice.dev.Vinit = HPGL_A3_init;
}

_HPGL_A4_devcpy()
{
	vdevice.dev = hpgldev;
	vdevice.dev.Vinit = HPGL_A4_init;
}

_HPGL_A1_devcpy()
{
	vdevice.dev = hpgldev;
	vdevice.dev.Vinit = HPGL_A1_init;
}

static DevEntry dxydev = {
	"dxy",
	"large",
	"small",
	noop,
	PLOT_char,
	noop,
	noop,
	PLOT_color,
	PLOT_draw,
	PLOT_exit,
	PLOT_fill,
	DXY_font,
	noop,
	noop,
	DXY_init,
	noop,
	noop,
#ifndef VOGLE
	noop,
	noop,
#endif
	PLOT_string,
	noop,
	noop
};

/*
 * _DXY_devcpy
 *
 *	copy the DXY device into vdevice.dev.
 */
_DXY_devcpy()
{
	vdevice.dev = dxydev;
}
