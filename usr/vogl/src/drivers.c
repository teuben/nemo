#include <stdio.h>
#include "vogl.h"
#include "vodevice.h"

extern char	*getenv();

struct vdev	vdevice;

static FILE	*fp;   /*  = stdout; */

static int	allocated = 0;

void		gexit();

/* device-independent function routines */

/*
 * voutput
 *
 *	redirect output - only for postscript, hpgl (this is not a feature)
 */
void
voutput(path)
	char	*path;
{
	char	buf[128];

	if ((fp = fopen(path, "w")) == (FILE *)NULL) {
		sprintf(buf, "voutput: couldn't open %s", path);
		verror(buf);
	}
}

/*
 * _voutfile
 *
 *	return a pointer to the current output file - designed for internal
 * use only.
 */
FILE *
_voutfile()
{
	return(fp);
}

/*
 * verror
 *
 *	print an error on the graphics device, and then exit. Only called
 * for fatal errors. We assume that stderr is always there.
 *
 */
void
verror(str)
	char	*str;
{
#ifdef MSWIN
	mswin_verror(str);
	if (vdevice.initialised)
		gexit();
#else
	if (vdevice.initialised)
		gexit();

	fprintf(stderr, "%s\n", str);
#endif
	exit(1);
}

void
viniterror(str)
	char	*str;
{
	fprintf(stderr, "%s: vogl not initialised\n", str);
	exit(1);
}

/*
 * gexit
 *
 *	exit the vogl/vogle system
 *
 */
void
gexit()
{
	if (!vdevice.initialised)
		verror("gexit: vogl not initialised");

	(*vdevice.dev.Vexit)();

	vdevice.devname = (char *)NULL;
	vdevice.initialised = 0;
	fp = stdout;
}

/*
 * getdevice
 *
 *	get the appropriate device table structure
 */
static void
getdevice(device)
	char	*device;
{
	char	buf[100];
#ifdef SUN
	if (strncmp(device, "sun", 3) == 0)
		_SUN_devcpy();
	else
#endif
#ifdef PIXRECT
	if (strncmp(device, "pixrect", 7) == 0)
		_PIXRECT_devcpy();
	else
#endif
#ifdef X11
	if (strncmp(device, "X11", 3) == 0)
		_X11_devcpy();
	else
#endif
#ifdef DECX11
	if (strncmp(device, "decX11", 6) == 0)
		_DECX11_devcpy();
	else
#endif
#ifdef NeXT
	if (strncmp(device, "NeXT", 4) == 0)
		_NeXT_devcpy();
	else
#endif
#ifdef POSTSCRIPT
	if (strncmp(device, "postscript", 10) == 0) {
		_PS_devcpy();
	} else
	if (strncmp(device, "ppostscript", 11) == 0) {
		_PSP_devcpy();
	} else
	if (strncmp(device, "cps", 3) == 0) {
		_CPS_devcpy();
	} else
	if (strncmp(device, "pcps", 4) == 0) {
		_PCPS_devcpy();
	} else
#endif
#ifdef HPGL
	if (strncmp(device, "hpgla1", 6) == 0)
		_HPGL_A1_devcpy();
	else if (strncmp(device, "hpgla3", 6) == 0)
		_HPGL_A3_devcpy();
	else if (strncmp(device, "hpgla4", 6) == 0)
		_HPGL_A4_devcpy();
	else if (strncmp(device, "hpgla2", 6) == 0 || strncmp(device, "hpgl", 4) == 0)
		_HPGL_A2_devcpy();
	else
#endif
#ifdef DXY
	if (strncmp(device, "dxy", 3) == 0)
		_DXY_devcpy();
	else
#endif
#ifdef TEK
	if (strncmp(device, "tek", 3) == 0)
		_TEK_devcpy();
	else
#endif
#ifdef GRX
	if (strncmp(device, "grx", 3) == 0)
		_grx_devcpy();
	else
#endif
#ifdef HERCULES
	if (strncmp(device, "hercules", 8) == 0)
		_hgc_devcpy();
	else
#endif
#ifdef MSWIN
	if (strncmp(device, "mswin", 5) == 0)
		_mswin_devcpy();
	else
#endif
#ifdef CGA
	if (strncmp(device, "cga", 3) == 0)
		_cga_devcpy();
	else
#endif
#ifdef EGA
	if (strncmp(device, "ega", 3) == 0)
		_ega_devcpy();
	else
#endif
#ifdef VGA
	if (strncmp(device, "vga", 3) == 0)
		_vga_devcpy();
	else
#endif
#ifdef SIGMA
	if (strncmp(device, "sigma", 5) == 0)
		_sigma_devcpy();
	else
#endif
	{
		if (*device == 0)
			sprintf(buf, "vogl: expected the enviroment variable VDEVICE to be set to the desired device.\n");
		else
			sprintf(buf, "vogl: %s is an invalid device type\n", device);
#ifdef MSWIN
		mswin_verror(buf);
#else
		fputs(buf, stderr);
		fprintf(stderr, "The devices compiled into this library are:\n");
#endif
#ifdef SUN
		fprintf(stderr, "sun\n");
#endif
#ifdef PIXRECT
		fprintf(stderr, "pixrect\n");
#endif
#ifdef X11
		fprintf(stderr, "X11\n");
#endif
#ifdef DECX11
		fprintf(stderr, "decX11\n");
#endif
#ifdef NeXT
		fprintf(stderr, "NeXT\n");
#endif
#ifdef POSTSCRIPT
		fprintf(stderr, "postscript\n");
		fprintf(stderr, "ppostscript\n");
		fprintf(stderr, "cps\n");
		fprintf(stderr, "pcps\n");
#endif
#ifdef HPGL
		fprintf(stderr, "hpgla1\n");
		fprintf(stderr, "hpgla2 (or hpgl)\n");
		fprintf(stderr, "hpgla3\n");
		fprintf(stderr, "hpgla4\n");
#endif
#ifdef DXY
		fprintf(stderr, "dxy\n");
#endif
#ifdef TEK
		fprintf(stderr, "tek\n");
#endif
#ifdef HERCULES
		fprintf(stderr, "hercules\n");
#endif
#ifdef CGA
		fprintf(stderr, "cga\n");
#endif
#ifdef EGA
		fprintf(stderr, "ega\n");
#endif
#ifdef VGA
		fprintf(stderr, "vga\n");
#endif
#ifdef SIGMA
		fprintf(stderr, "sigma\n");
#endif
#ifdef GRX
		fprintf(stderr, "grx\n");
#endif
		exit(1);
	}
}

/*
 * vinit
 *
 * 	Just set the device name. ginit and winopen are basically
 * the same as the VOGLE the vinit function.
 *
 */
void
vinit(device)
	char	*device;
{
	vdevice.devname = device;
}

/*
 * winopen
 *
 *	use the more modern winopen call (this really calls ginit),
 * we use the title if we can
 */
long
winopen(title)
	char	*title;
{
	vdevice.wintitle = title;

	ginit();

	return(1L);
}

/*
 * ginit
 *
 *	by default we check the environment variable, if nothing
 * is set we use the value passed to us by the vinit call.
 */
void
ginit()
{
	char	*dev;
	int	i;

	if (vdevice.devname == (char *)NULL) {
		if ((dev = getenv("VDEVICE")) == (char *)NULL)
			getdevice("");
		else
			getdevice(dev);
	} else 
		getdevice(vdevice.devname);

	if (vdevice.initialised)
		gexit();

	if (!allocated) {
		allocated = 1;
		deflinestyle(0, 0xffff);
		vdevice.transmat = (Mstack *)vallocate(sizeof(Mstack));
		vdevice.transmat->back = (Mstack *)NULL;
		vdevice.attr = (Astack *)vallocate(sizeof(Astack));
		vdevice.attr->back = (Astack *)NULL;
		vdevice.viewport = (Vstack *)vallocate(sizeof(Vstack));
		vdevice.viewport->back = (Vstack *)NULL;
		vdevice.bases = (Matrix *)vallocate(sizeof(Matrix) * 10);
		vdevice.enabled = (char *)vallocate(MAXDEVTABSIZE);
	}

	for (i = 0; i < MAXDEVTABSIZE; i++)
		vdevice.enabled[i] = 0;

	/* NOTE:
	 * There is a slight behaviour change from previous versions of VOGL
	 * if you define FIRST_REDRAW... you always get a REDRAW event as the
	 * first event in the queue.
	 */
#define FIRST_REDRAW 1
#ifdef FIRST_REDRAW
	/*
	 * Arrange for a REDRAW to be the first thing in the queue...
         * (winopen always enters a REDRAW in the queue on a real SGI).
	 */
	vdevice.alreadyread = TRUE;
	vdevice.data = 0;
	vdevice.devno = REDRAW;
	vdevice.enabled[REDRAW / 8] |= (1 << (REDRAW & 0x7));
#else
	vdevice.alreadyread = FALSE;
	vdevice.data = 0;
	vdevice.devno = 0;
#endif
	vdevice.kbdmode = vdevice.mouseevents = vdevice.kbdevents = 0;

	vdevice.concave = 0;
	vdevice.clipoff = 0;
	vdevice.sync = 1;
	vdevice.cpW[V_W] = 1.0;			/* never changes */

	vdevice.maxfontnum = 2;

	vdevice.attr->a.lw = 1;
	vdevice.attr->a.fontnum = 0;
	vdevice.attr->a.mode = 0;
	vdevice.attr->a.backface = 0;

	if ((*vdevice.dev.Vinit)()) {
		vdevice.initialised = 1;

		vdevice.inobject = 0;
		vdevice.inpolygon = 0;

		viewport((Screencoord)0, (Screencoord)vdevice.sizeSx - 1,
			(Screencoord)0, (Screencoord)vdevice.sizeSy - 1);

		ortho2(0.0, (Coord)(vdevice.sizeSx - 1), 0.0, (Coord)(vdevice.sizeSy - 1));

		move(0.0, 0.0, 0.0);

		font(0);	/* set up default font */

	} else {
		fprintf(stderr, "vogl: error while setting up device\n");
		exit(1);
	}

	setlinestyle(0);
}

/*
 * gconfig
 *
 *	thankfully a noop.
 */
void
gconfig()
{
}

/*
 * Hacky new device changing routines...
 */
#define	DEVSTACK	8
static	VDevice	vdevstk[DEVSTACK];
static	int	vdevindx = 0;

void
pushdev(device)
	char	*device;
{
	/*
	 * Save the old vdevice structure
	 */
	pushattributes();
	pushviewport();

	if (vdevindx < DEVSTACK)
		vdevstk[vdevindx++] = vdevice;
	else
		verror("vogl: pushdev: Device stack overflow");

	vdevice.initialised = 0;

	getdev(device);

	(*vdevice.dev.Vinit)();

	vdevice.initialised = 1;

	popviewport();
	popattributes();
}

void
popdev()
{
	/*
	 * Restore the old vdevice structure
	 */
	pushattributes();
	pushviewport();

	(*vdevice.dev.Vexit)();
	if (vdevindx > 0)
		vdevice = vdevstk[--vdevindx];
	else
		verror("vogl: popdev: Device stack underflow");

	popviewport();
	popattributes();
}
/*
 * vnewdev
 *
 * reinitialize vogl to use a new device but don't change any
 * global attributes like the window and viewport settings.
 */
void
vnewdev(device)
	char	*device;
{
	if (!vdevice.initialised)
		verror("vnewdev: vogl not initialised\n");

	pushviewport();	

	(*vdevice.dev.Vexit)();

	vdevice.initialised = 0;

	getdevice(device);

	(*vdevice.dev.Vinit)();

	vdevice.initialised = 1;

	/*
	 * Need to update font for this device...
	 */
	font(vdevice.attr->a.fontnum);

	popviewport();
}

/*
 * vgetdev
 *
 *	Returns the name of the current vogl device 
 *	in the buffer buf. Also returns a pointer to
 *	the start of buf.
 */
char	*
vgetdev(buf)
	char	*buf;
{
	/*
	 * Note no exit if not initialized here - so that gexit
	 * can be called before printing the name.
	 */
	if (vdevice.dev.devname)
		strcpy(buf, vdevice.dev.devname);
	else
		strcpy(buf, "(no device)");

	return(&buf[0]);
}

/*
 * getvaluator
 *
 *	similar to the VOGLE locator only it returns either x (MOUSEX) or y (MOUSEY).
 */
long
getvaluator(dev)
	Device	dev;
{
	int	a, b, c;

	if (!vdevice.initialised)
		verror("getvaluator: vogl not initialised");

	c = (*vdevice.dev.Vlocator)(&a, &b);

	if (c != -1) {
		if (dev == MOUSEX)
			return((long)a);
		else 
			return((long)b);
	}

	return(-1);
}

/*
 * getbutton
 *
 *	returns the up (or down) state of a button. 1 means down, 0 up,
 * -1 invalid.
 */
Boolean
getbutton(dev)
	Device	dev;
{
	int	a, b, c;

	if (dev < 256) {
		c = (*vdevice.dev.Vcheckkey)();
		if (c >= 'a' && c <= 'z')
			c = c - 'a' + 'A';
		if (c == dev)
			return(1);
		return(0);
	} else if (dev < 261) {
		c = (*vdevice.dev.Vlocator)(&a, &b);
		if (c & 0x01 && dev == MOUSE3)
			return(1);
		if (c & 0x02 && dev == MOUSE2)
			return(1);
		if (c & 0x04 && dev == MOUSE1)
			return(1);
		return(0);
	}

	return(-1);
}

/*
 * Get the values of the valuators in devs and put them into vals
 */
void 
getdev(n, devs, vals)
	long n;
	Device devs[];
	short vals[];
{
	int	i;

	for( i=0; i < n; i++)
		vals[i] = (short)getvaluator(devs[i]);
}


/*
 * clear
 *
 *	clears the screen to the current colour, excepting devices
 * like a laser printer where it flushes the page.
 *
 */
void
clear()
{
	Token	*tok;

	if (!vdevice.initialised)
		verror("clear: vogl not initialised");

	if (vdevice.inobject) {
		tok = newtokens(1);
		tok->i = CLEAR;

		return;
	}

	(*vdevice.dev.Vclear)();
}

/*
 * colorf
 *
 *	set the current colour to colour index given by
 * the rounded value of f.
 *
 */
void
colorf(f)
	float	f;
{
	color((int)(f + 0.5));
}

/*
 * color
 *
 *	set the current colour to colour index number i.
 *
 */
void
color(i)
	int	i;
{
	Token	*tok;

	if (!vdevice.initialised)
		verror("color: vogl not initialised");

	if (vdevice.inobject) {
		tok = newtokens(2);

		tok[0].i = COLOR;
		tok[1].i = i;
		return;
	}

	vdevice.attr->a.color = i;
	(*vdevice.dev.Vcolor)(i);
}

long
getcolor()
{
	return((long)vdevice.attr->a.color);
}

/*
 * mapcolor
 *
 *	set the color of index i.
 */
void
mapcolor(i, r, g, b)
	Colorindex	i;
	short		r, g, b;
{
	Token	*tok;

	if (!vdevice.initialised)
		verror("mapcolor: vogl not initialised");

	if (vdevice.inobject) {
		tok = newtokens(5);

		tok[0].i = MAPCOLOR;
		tok[1].i = i;
		tok[2].i = r;
		tok[3].i = g;
		tok[4].i = b;

		return;
	}

	(*vdevice.dev.Vmapcolor)(i, r, g, b);
}

/*
 * getplanes
 *
 *	Returns the number if bit planes on a device.
 */
long
getplanes()
{
	if (!vdevice.initialised)
		verror("getdepth: vogl not initialised\n");

	return((long)vdevice.depth);
}

/*
 * reshapeviewport
 *
 *	Simply sets the viewport to the size of the current window
 */
void
reshapeviewport()
{
	viewport(0, vdevice.sizeSx - 1, 0, vdevice.sizeSy - 1);
}

/*
 * winconstraints
 *		- does nothing
 */
void
winconstraints()
{
}

/*
 * keepaspect
 *		- does nothing
 */
void
keepaspect()
{
}

/*
 * shademodel
 *		- does nothing
 */
void
shademodel(model)
	long	model;
{
}

/*
 * getgdesc
 *
 *	Inquire about some stuff....
 */
long
getgdesc(inq)
	long	inq;
{	
	/*
	 * How can we know before the device is inited??
	 */

	switch (inq) {
	case GD_XPMAX:
		if (vdevice.initialised)
			return((long)vdevice.sizeSx);
		else
			return(500L);	/* A bullshit number */
	case GD_YPMAX:
		if (vdevice.initialised)
			return((long)vdevice.sizeSy);
		else
			return(500L);
	default:
		return(-1L);
	}
}

/*
 * foregound
 * 		Dummy - does nothing.
 */
void
foreground()
{
}

/*
 * vsetflush
 *
 * Controls flushing of the display - we can get considerable
 * Speed up's under X11 using this...
 */
void
vsetflush(yn)
	int	yn;
{
	vdevice.sync = yn;
}

/*
 * vflush
 *
 * Explicitly call the device flushing routine...
 * This is enabled for object so that you can force an update
 * in the middle of an object, as objects have flushing off
 * while they are drawn anyway.
 */
void
vflush()
{
	Token	*tok;

	if (!vdevice.initialised)
		verror("vflush: vogl not initialised");

	if (vdevice.inobject) {
		tok = newtokens(1);
		tok->i = VFLUSH;

		return;
	}

	(*vdevice.dev.Vsync)();
}


/* 
 * getorigin
 *
 *	Returns the origin of the window. This is a dummy.
 */
void
getorigin(x, y)
	long	*x, *y;
{
	*x = *y = 0;
}

/*
 * getsize
 *
 *	Returns the approximate size of the window (some window managers
 *	stuff around with your borders).
 */
void
getsize(x, y)
	long	*x, *y;
{
	*x = (long)vdevice.sizeSx;
	*y = (long)vdevice.sizeSy;
}

winattach()
{}

winset()
{}

