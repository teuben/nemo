/*
 * This file is directly from the VOGLE library. It's simply had a few
 * gratuitous name changes and some comments added.
 */
#ifdef SGI
#include <gl.h>
#include <device.h>
#else
#include "vogl.h"
#include "vodevice.h"
#endif

#include <stdio.h>
#include <fcntl.h>
#ifdef TC
extern double sin();
extern double cos();
#else
#include <math.h>
#endif
#include <string.h>

#ifndef PATH_SIZE
#define	PATH_SIZE	256
#endif

#ifdef PC
#ifndef FONTLIB
#define	FONTLIB		"c:\\lib\\hershey\\"
#endif
#else
#ifndef FONTLIB
#define	FONTLIB		"/usr/local/lib/hershey/"
#endif
#endif

#define	ABS(a)		((a) < 0 ? -(a) : (a))
#define	MAX(a, b)	((a) < (b) ? (b) : (a))
#define	XCOORD(x)	(int)((x) - 'R')
#define	YCOORD(y)	(int)('R' - (y))

#ifndef PI
#define PI 	3.14159265358979
#endif

#ifndef D2R
#define D2R	(PI / 180.0)
#endif

#define LEFT		0	/* The default */
#define CENTERED	1
#define RIGHT		2

static	float	tcos = 1.0, tsin = 0.0;			/* For rotations */
static	float	SCSIZEX = 1.0, SCSIZEY = 1.0;		/* Scale factors */
static	int	Justify = LEFT;
static	int	Fixedwidth = 0;			/* Some flags	 */
static	short	nchars;					/* No. in font	 */
int	hLoaded = 0;

static	hershfont();

static	char	old_font[PATH_SIZE] = "";	/* So we don't have to reload it */
static	char	fpath[PATH_SIZE] = "";	


static	struct	{
	int	as;	/* Max ascender of a character in this font */
	int	dec;	/* Max decender of a character in this font */
	int	mw;	/* Max width of a character in this font */
	char	*p;	/* All the vectors in the font */
	char	**ind;	/* Pointers to where the chars start in p */
} ftab;

float	hstrlength();	/* Length of a set of Hershey characters  */
float	hgetfontwidth(); /* The width of this font */
float	hgetfontheight(); /* The height of this font */
extern	char	*hallocate();
extern	char	*getenv();

/*
 * hfont
 * 	loads in a hershey font.
 */
void
hfont(name)
	char	*name;
{
	/*
	 * check we aren't loading the same font twice in a row
	 */
	if (*name == '/') {
		if (strcmp(strrchr(name, '/') + 1, old_font) == 0)
			return;

	} else if (strcmp(name, old_font) == 0)
		return;

	/*
	 * Try and load it
	 */
	if (!hershfont(name)) {
		fprintf(stderr, "hershlib: problem reading font file '%s'.\n", name);
		exit(1);
	}

	/*
	 * Save the name of it...
	 */
	if (*name == '/')
		strcpy(old_font, strrchr(name, '/') + 1);
	else 
		strcpy(old_font, name);

}

/*
 * hnumchars
 *
 *	Return the number of characters in the currently loaded hershey font.
 */
int
hnumchars()
{

	check_loaded("hnumchars");

	return((int)nchars);
}

/*
 * hsetpath
 *
 *	Set the path of the directory to look for fonts.
 */
void
hsetpath(path)
	char	*path;
{
	int	l;

	strcpy(fpath, path);

	l = strlen(fpath);

#ifdef PC
	if (fpath[l] != '\\')
		strcat(fpath, "\\");
#else
	if (fpath[l] != '/')
		strcat(fpath, "/");
#endif
}

/*
 * hershfont
 *
 * Load in a hershey font. First try the environment, then the font library,
 * if that fails try the current directory, otherwise return 0.
 */
static int
hershfont(fontname)
	char	*fontname;
{
	FILE	*fp, *fopen();
	int	i, j;
	short	nvects, n;
	char	*flib;
	char	path[PATH_SIZE];	
	
	if (fpath[0] != '\0') {
		strcpy(path, fpath);
		strcat(path, fontname);
	} else if ((flib = getenv("HFONTLIB")) != (char *)NULL) {
		strcpy(path, flib);
#ifdef PC
		strcat(path, "\\");
#else
		strcat(path, "/");
#endif
		strcat(path, fontname);
	} else if ((flib = getenv("VFONTLIB")) != (char *)NULL) {
		strcpy(path, flib);	/* To be compatible with VOGLE */
#ifdef PC
		strcat(path, "\\");
#else
		strcat(path, "/");
#endif
		strcat(path, fontname);
	} else {
		strcpy(path, FONTLIB);
#ifdef PC
		strcat(path, "\\");
#else
		strcat(path, "/");
#endif
		strcat(path, fontname);
	}

#ifdef PC
	if ((fp = fopen(path, "r+b")) == (FILE *)NULL) 
		if ((fp = fopen(fontname, "r+b")) == (FILE *)NULL) {
#else
	if ((fp = fopen(path, "r")) == (FILE *)NULL)
		if ((fp = fopen(fontname, "r")) == (FILE *)NULL) {
#endif
			fprintf(stderr, "hershlib: Can't open Hershey fontfile '%s' or './%s'.\n", path, fontname);
			exit(1);
		}

	if (fread(&nchars, sizeof(nchars), 1, fp) != 1)
		return (0);

#ifdef DEBUG
	printf("nchars = %d\n", nchars);
#endif

	if (fread(&nvects, sizeof(nvects), 1, fp) != 1)
		return(0);

#ifdef DEBUG
	printf("nvects = %d\n", nvects);
#endif
	if (fread(&n, sizeof(n), 1,  fp) != 1)
		return(0);

#ifdef DEBUG
	printf("ftab.as = %d\n", n);
#endif
	ftab.as = (int)n;

	if (fread(&n, sizeof(n), 1, fp) != 1)
		return(0);

#ifdef DEBUG
	printf("ftab.dec = %d\n", n);
#endif
	ftab.dec = (int)n;

	if (fread(&n, sizeof(n), 1, fp) != 1)
		return(0);

#ifdef DEBUG
	printf("ftab.mw = %d\n", n);
#endif

	ftab.mw = (int)n;

	/*
	 *  Allocate space for it all....
	 */
	if (hLoaded) {
		if (ftab.ind[0])
			free(ftab.ind[0]);
		if (ftab.ind)
			free(ftab.ind);
		hLoaded = 0;
	}

	ftab.ind = (char **)hallocate(sizeof(char *)*(nchars + 1));

	ftab.p = (char *)hallocate((unsigned)(2 * nvects));

	/*
	 *  As we read in each character, figure out what ind should be
	 */

	for (i = 0; i < nchars; i++) {
		if (fread(&n , sizeof(n), 1, fp) != 1)
			return(0);

		if (fread(ftab.p, 1, (unsigned)n, fp) != (unsigned)n)
			return(0);

		ftab.ind[i] = ftab.p;
		ftab.p += n;
	}

	ftab.ind[nchars] = ftab.p;	/* To Terminate the last one */

	fclose(fp);
	hLoaded = 1;
	return(1);
}

/*
 * hgetcharsize
 *
 *	get the width and height of a single character. At the moment, for
 * the hershey characters, the height returned is always that of the
 * difference between the maximun descender and ascender.
 *
 */
void
hgetcharsize(c, width, height)
	char	c;
	float	*width, *height;
{
	check_loaded("hgetcharsize");
	
	*height = (float)(ftab.as - ftab.dec) * SCSIZEY;

	if (Fixedwidth)
		*width = (float)ftab.mw * SCSIZEX;
	else
		*width = (float)(ftab.ind[c - 32][1] - ftab.ind[c - 32][0]) * SCSIZEX;
}

/*
 * hdrawchar
 *
 * Display a character from the currently loaded font.
 */
void
hdrawchar(c)
	int	c;
{
	char	*p, *e;
	int	Move, i, x, y, xt, yt;
	float	xs, ys, xp, yp, tmp, xtmp, ytmp;
#ifdef VOGL
	int	sync;
#endif

	check_loaded("hdrawchar");

	if ((i = c - 32) < 0)
		i = 0;
	if (i >= nchars)
		i = nchars - 1;

	Move = 1;

	xt = yt = 0;
	if (Justify == LEFT) {
		xt = (Fixedwidth ? -ftab.mw / 2 : XCOORD(ftab.ind[i][0]));
		yt = ftab.dec;
	} else if (Justify == RIGHT) {
		xt = (Fixedwidth ?  ftab.mw / 2 : -XCOORD(ftab.ind[i][0]));
		yt = ftab.dec;
	}

	e = ftab.ind[i + 1];
	p = ftab.ind[i] + 2;

	xtmp = ytmp = 0.0;

#ifdef VOGL
	if (sync = vdevice.sync)
		vdevice.sync = 0;
#endif

	while(p < e) {
		x = XCOORD((int)*p++);
		y = YCOORD((int)*p++);
		if (x != -50) {			/* means move */
			xp = (float)(x - xt) * SCSIZEX;
			yp = (float)(y - yt) * SCSIZEY;
			tmp = xp;
			xp = tcos*tmp - tsin*yp;
			yp = tsin*tmp + tcos*yp;
			xs = xp - xtmp;
			ys = yp - ytmp;
			xtmp = xp;
			ytmp = yp;
			if (Move) {
				Move = 0;
				rmv((Coord)xs, (Coord)ys, (Coord)0.0);
			} else
				rdr((Coord)xs, (Coord)ys, (Coord)0.0);

		} else {
			Move = 1;
		}
	}
	/*
	 * Move to right hand of character.
	   BLARK
	 */
	
	tmp = Fixedwidth ? (float)ftab.mw : (float)(ftab.ind[i][1] - ftab.ind[i][0]);

	tmp *= SCSIZEX;
	xs = tcos * tmp - xtmp;
	ys = tsin * tmp - ytmp;
	rmv((Coord)xs, (Coord)ys, 0.0);

#ifdef VOGL
	if (sync) {
		vdevice.sync = 1;
		(*vdevice.dev.Vsync)();
	}
#endif
}

/*
 * htextsize
 *
 * set software character scaling values 
 *
 * Note: Only changes software char size. Should be called
 * after a font has been loaded.
 *
 */
void
htextsize(width, height)
	float	width, height;
{
	float	a;

	check_loaded("htextsize");

	a = (float)MAX((int)ftab.mw, (int)(ftab.as - ftab.dec));

	SCSIZEX = width / ABS(a);
	SCSIZEY = height / ABS(a);
}

/*
 * hgetfontwidth
 *
 * Return the maximum Width of the current font.
 *
 */
float
hgetfontwidth()
{
	check_loaded("hgetfontwidth");

	return((float)(SCSIZEX * MAX((int)ftab.mw, (int)(ftab.as - ftab.dec))));
}

/* 
 * hgetfontheight
 *
 * Return the maximum Height of the current font
 */
float 
hgetfontheight()
{
	check_loaded("hgetfontheight");

	return((float)(SCSIZEY * MAX((int)ftab.mw, (int)(ftab.as - ftab.dec))));
}

/*
 * hgetfontsize
 *
 * Get the current character size in user coords.
 * For software Hershey fonts, the character width is that of
 * a the widest character and the height the height of the tallest.
 *
 */
void
hgetfontsize(cw, ch)
	float 	*cw, *ch;
{
	check_loaded("hgetfontsize");

	*cw = hgetfontwidth();
	*ch = hgetfontheight();
}

/*
 * hgetdecender
 *
 *	Return the maximum decender of the current font.
 * 	(In world coords).
 */
float
hgetdecender()
{
	check_loaded("hgetdecender");

	return((float)ftab.dec * SCSIZEY);
}

/*
 * hgetascender
 *
 *	Return the maximum assender of the current font.
 * 	(In world coords).
 */
float
hgetascender()
{
	check_loaded("hgetascender");

	return((float)ftab.as * SCSIZEY);
}

/*
 * hcharstr
 *
 * Draw a string from the current pen position.
 *
 */
void
hcharstr(string)
	char 	*string;
{
	float	width, height, cx, cy, cz;
	char	*str = string, c;
	int	oldJustify;
#ifdef VOGL
	int	sync;
#endif

	check_loaded("hcharstr");

	height = hgetfontheight();
	width = hstrlength(string);

	cx = cy = 0.0;
	if (Justify == CENTERED) {
		height /= 2.0;
		width /= 2.0;
		cx =  height * tsin - width * tcos;
		cy = -height * tcos - width * tsin;
	} else if (Justify == RIGHT) {
		height = 0.0;
		cx =  height * tsin - width * tcos;
		cy = -height * tcos - width * tsin;
	}
	rmv((Coord)cx, (Coord)cy, (Coord)0.0);


	/*
	 * For the duration of hershey strings, turn off
	 * "Justify" as we have already compensated
	 * for it in hcharstr()
	 */
	oldJustify = Justify;
	Justify = 0;

	/*
	 * Now display each character
	 *
	 */
#ifdef VOGL
	if (sync = vdevice.sync)
		vdevice.sync = 0;
#endif

	while (c = *str++)
		hdrawchar(c);

#ifdef VOGL
	if (sync) {
		vdevice.sync = 1;
		(*vdevice.dev.Vsync)();
	}
#endif
	
	Justify = oldJustify;

}

/*
 * istrlength
 *
 * Find out the length of a string in raw "Hershey coordinates".
 */
static	int
istrlength(s)
	char	*s;
{
	char	c;
	int	i, len = 0;
	
	if (Fixedwidth)
		return((int)(strlen(s) * ftab.mw));
	else {
		while (c = *s++) {
			if ((i = (int)c - 32) < 0 || i >= nchars)
				i = nchars - 1;

			len += (ftab.ind[i][1] - ftab.ind[i][0]);
		}
		return (len);
	}
}

/*
 * hstrlength
 *
 * Find out the length (in world coords) of a string.
 *
 */
float
hstrlength(s)
	char	*s;
{
	check_loaded("hstrlength");

	return((float)(istrlength(s) * SCSIZEX));
}

/*
 * hboxtext
 *
 * Draw text so it fits in a "box" - note only works with hershey text
 */
void
hboxtext(x, y, l, h, s)
	float	x, y, l, h;
	char	*s;
{
	float	oscsizex, oscsizey, cz;

	check_loaded("hboxtext");

	oscsizex = SCSIZEX;
	oscsizey = SCSIZEY;
	/*
	 * set width so string length is the same a "l" 
	 */
	SCSIZEX = l / (float)istrlength(s);

	/* 
	 * set character height so it's the same as "h" 
	 */
	SCSIZEY = h / (float)(ftab.as - ftab.dec);

	move2(x, y);
	hcharstr(s);

	SCSIZEX = oscsizex;
	SCSIZEY = oscsizey;
}

/*
 * hboxfit
 *
 * Set up the scales etc for text so that a string of "nchars" characters
 * of the maximum width in the font fits in a box.
 */
void
hboxfit(l, h, nchars)
	float	l, h;
	int	nchars;
{
	check_loaded("hboxfit");

	SCSIZEX = l / (float)(nchars * ftab.mw);
	SCSIZEY = h / (float)(ftab.as - ftab.dec);
}

/*
 * The following can be set without a font actually being loaded....
 */
/*
 * hcenter
 *
 *	Turns centering of text on or off
 */
void
hcentertext(onoff)
	int	onoff;
{
	Justify = (onoff != 0 ? CENTERED : LEFT);
}

/*
 * hrightjustify
 *
 *	Right Justifies the text.
 */
void
hrightjustify(onoff)
	int	onoff;
{
	Justify = (onoff != 0 ? RIGHT : LEFT);
}

/*
 * hleftjustify
 *
 *	Left Justifies the text. (the default).
 */
void
hleftjustify(onoff)
	int	onoff;
{
	Justify = (onoff != 0 ? LEFT : RIGHT);
}

/*
 * fixedwidth
 *
 *	Turns fixedwidth text on or off
 */
void
hfixedwidth(onoff)
	int	onoff;
{
	Fixedwidth = onoff;
}

/*
 * htextang
 *
 * set software character angle in degrees
 *
 * strings will be written along a line 'ang' degrees from the 
 * horizontal screen direction
 *
 * Note: only changes software character angle
 *
 */
void
htextang(ang)
	float	ang;
{

	tcos = cos((double)(ang * D2R));
	tsin = sin((double)(ang * D2R));
}
