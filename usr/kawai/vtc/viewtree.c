#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "vtc.h"

#if USEX11
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif /* USEX11 */

#if USEGD
#include <gd.h>
#include <gdfontg.h>
#include <gdfonts.h>
#include <gdfontmb.h>
#endif /* USEGD */

#include <vtclocal.h>

#define BORDER (2)
#define WIDTH (512)
#define HEIGHT (512)
/*
#define WIDTH (700)
#define HEIGHT (700)
*/

typedef enum
{
    BLACK = 0,
    WHITE,
    YELLOW,
    GREEN,
}
Color;
typedef enum
{
    CELLONLY = 0,
    BODYONLY,
    BODY_AND_CELL,
}
Showmode;
typedef enum
{
    LISTBODYONLY = 0,
    LISTONLY,
    LIST_AND_BODY,
    BODYONLY2,
}
Showmode2;

typedef struct polyh_t
{
    int n; /* # of vertices */
    double (*v3)[3]; /* 3D vertice coordinate */
    int (*v2)[2]; /* 2D vertice coordinate */
    int ns; /* # of sides */
    int (*s)[2]; /* side start & end vertice */
} Polyh;

typedef struct
{
    int nb; /* # of bodies */
    double (*bpos)[3]; /* body positions */
    int nc; /* # of nodes */
    double (*cpos)[3]; /* node positions. (0,0,0) for root */
    double *csize; /* side len */
} Snapshot;


typedef struct imgcontext_t
{
#if USEX11
    /* X related */
    Display *dpy;
    Window w;
    Window root;
    int screen;
    XEvent event;
    GC gcaxis;
    GC gcline;
    GC gcline2;
    GC gcclear;
    unsigned long black, white;
#endif /* USEX11 */
    /* point of view */
    double theta;
    double phi;
    double zoom;
    Polyh axis;
    Showmode showmode;
} ImgContext;

static int will_quit;
static int will_dump = 0;

static void initXContext(ImgContext *ic);
static void readData(char *fname, Polyh *p);
static void mainLoop0(ImgContext *ic, Polyh *bp, Polyh *cp);
static void exposeEvent0(ImgContext *ic, Polyh *bp, Polyh *cp);
static void mainLoop1(ImgContext *ic, Polyh *ballp, Polyh *bip, Polyh *bjp, Polyh *cip, Polyh *cjp);
static void exposeEvent1(ImgContext *ic, Polyh *ballp, Polyh *bip, Polyh *bjp, Polyh *cip, Polyh *cjp);
static void keypressEvent(ImgContext *ic);
static void projectPolyh(ImgContext *ic, Polyh *p, double zoom);
#if USEX11
static void displayAxis(ImgContext *ic, Polyh *p, GC gc);
static void displayPolyh(ImgContext *ic, Polyh *p, GC gc);
#endif /* USEX11 */

static void clearDrawable(ImgContext *ic);
static void body2polyh(Snapshot *t, Polyh *p);
static void treenode2polyh(Snapshot *t, Polyh *p);
static void freepolyh(Polyh *p);

#if USEGD
static void dumpimage0(ImgContext *ic, Polyh *bp, Polyh *cp);
static void dumpimage1(ImgContext *ic, Polyh *ballp, Polyh *bip, Polyh *bjp, Polyh *cip, Polyh *cjp);
static void dumpAxis(gdImagePtr im, ImgContext *ic, Polyh *p, int c);
static void dumpPolyh(gdImagePtr im, ImgContext *ic, Polyh *p, int c);
#endif /* USEGD */


static int ilen, jlen;

#if USEGD

#if 1 /* high quality image */

#define IMG_WIDTH (768)
#define IMG_HEIGHT (768)
#define IMG_NCOLORS (256)
void
vtc_plotstar(char *fname, Nbodyinfo *nb, char *msg, double scale, double center[2], double cmass,
	      double *xheavy, double *xmin)
{
    FILE *fp, *fpw;
    double pixelvalmax, pixelvalmin, logmax, logmin;
    int pe, i, i0, j, jmax, k;
    int n = nb->n;
    double *buf = NULL;
    int nparticle;
    int bufsize;
    double r[3];
    double x, y, w;
    gdImagePtr im;
    int ix, iy;
    static double pixelval[IMG_HEIGHT][IMG_WIDTH];
    static int color[IMG_NCOLORS];
    char textbuf[255];

    /* project mass to a 2D plane */
    for (ix = 0; ix < IMG_WIDTH; ix ++) {
	for (iy = 0; iy < IMG_HEIGHT; iy ++) {
	    pixelval[iy][ix] = 0.0;
	}
    }
    i0 = 0;
    for (i = 0; i < n; i++) {
	if (i > i0) {
	    fprintf(stderr, "%d particles done\n", i);
	    i0 += 1000000;
	}
	r[0] = nb->x[i][0] - center[0];
	r[1] = nb->x[i][1] - center[1];
	r[2] = nb->x[i][2];
	x = (r[0]*scale*0.5+0.5)*IMG_WIDTH;
	y = (r[1]*scale*0.5+0.5)*IMG_HEIGHT;
	if (0 <= x && x < IMG_WIDTH && 0 <= y && y < IMG_HEIGHT) {
	    if (nb->m[i] < cmass) {
		pixelval[(int)y][(int)x] += nb->m[i];
	    }
	}
    }

    /* assign color to each pixel */
    im = gdImageCreate(IMG_WIDTH, IMG_HEIGHT);
    for (i = 0; i < IMG_NCOLORS/2; i++) {
	color[i] = gdImageColorAllocate(im, 0, 0, i);
    }
    for (i = 0; i < IMG_NCOLORS/2; i++) {
	color[i+IMG_NCOLORS/2] = gdImageColorAllocate(im, i*2, i*2, i+IMG_NCOLORS/2);
    }
    gdImageFill(im, IMG_WIDTH, IMG_HEIGHT, color[0]);
    pixelvalmax = 0.0;
    pixelvalmin = 1e30;
    for (ix = 0; ix < IMG_WIDTH; ix ++) {
	for (iy = 0; iy < IMG_HEIGHT; iy ++) {
	    if (pixelvalmax < pixelval[iy][ix]) {
		pixelvalmax = pixelval[iy][ix];
	    }
	    if (pixelval[iy][ix] > 0.0 && pixelvalmin > pixelval[iy][ix]) {
		pixelvalmin = pixelval[iy][ix];
	    }
	}
    }
    fprintf(stderr, "pixelval max: %f  min: %f\n", pixelvalmax, pixelvalmin);
    logmax = log(pixelvalmax);
    logmin = log(pixelvalmin);
    fprintf(stderr, "logmax: %f  logmin: %f\n", logmax, logmin);
    for (ix = 0; ix < IMG_WIDTH; ix ++) {
	for (iy = 0; iy < IMG_HEIGHT; iy ++) {
	    int cindex;
	    double p = pixelval[iy][ix];
	    if (p > 0.0) {
		cindex = (log(p)-logmin)/(logmax-logmin)*IMG_NCOLORS;
	    }
	    else {
		cindex = 0;
	    }
	    cindex = (cindex < IMG_NCOLORS ? cindex : IMG_NCOLORS-1);
	    gdImageSetPixel(im, ix, iy, color[cindex]);
	}
    }

    /* point the position of the heavy particle closest to the center */
    if (xheavy != NULL) {
	r[0] = xheavy[0] - center[0];
	r[1] = xheavy[1] - center[1];
	r[2] = xheavy[2];
	ix = (r[0]*scale*0.5+0.5)*IMG_WIDTH;
	iy = (r[1]*scale*0.5+0.5)*IMG_HEIGHT;
	gdImageLine(im, ix-2, iy, ix+2, iy, IMG_NCOLORS-1);
	gdImageLine(im, ix, iy-2, ix, iy+2, IMG_NCOLORS-1);

	r[0] = xmin[0] - center[0];
	r[1] = xmin[1] - center[1];
	r[2] = xmin[2];
	ix = (r[0]*scale*0.5+0.5)*IMG_WIDTH;
	iy = (r[1]*scale*0.5+0.5)*IMG_HEIGHT;
	gdImageLine(im, ix-2, iy-2, ix+2, iy+2, IMG_NCOLORS-1);
	gdImageLine(im, ix-2, iy+2, ix+2, iy-2, IMG_NCOLORS-1);
    }

    /* draw a gauge */
    w = 0.1;
    sprintf(textbuf, "%4.2f Mpc", 2.0*w/scale);
    gdImageLine(im,
		IMG_WIDTH*0.1, IMG_HEIGHT-45+7,
		IMG_WIDTH*(0.1+w), IMG_HEIGHT-45+7,
		color[IMG_NCOLORS-1]);
    gdImageLine(im,
		IMG_WIDTH*0.1, IMG_HEIGHT-45-2+7,
		IMG_WIDTH*0.1, IMG_HEIGHT-45+2+7,
		color[IMG_NCOLORS-1]);
    gdImageLine(im,
		IMG_WIDTH*(0.1+w), IMG_HEIGHT-45-2+7,
		IMG_WIDTH*(0.1+w), IMG_HEIGHT-45+2+7,
		color[IMG_NCOLORS-1]);
    gdImageString(im,
		  gdFontMediumBold,
		  IMG_WIDTH*(0.1+w)+10, IMG_HEIGHT-45,
		  textbuf, color[IMG_NCOLORS-1]);

    /* put a text in msg */
    gdImageString(im,
		  gdFontMediumBold,
		  IMG_WIDTH*3.0/4.0, IMG_HEIGHT-45,
		  msg, color[IMG_NCOLORS-1]);

    fpw = fopen(fname, "w");
    if (!fpw) {
	perror(__FILE__ " write_image");
	exit (2);
    }
    gdImageGif(im, fpw);
    gdImageDestroy(im);
    fclose(fpw);
}

#else /* poor image */

void
vtc_plotstar(double (*pos)[3], int n, double time, double ratio, char *fname)
{
    int i, ii;
    gdImagePtr im;
    int black, yellow, red;
    double x, y;
    FILE *fp;
    static int imgno = 0;

    ratio *= 4.0;
    im = gdImageCreate(WIDTH, HEIGHT);
    black = gdImageColorAllocate(im, 0, 0, 0);
    yellow = gdImageColorAllocate(im, 255, 255, 0);
    red = gdImageColorAllocate(im, 255, 0, 0);
    gdImageFill(im, WIDTH, HEIGHT, black);

#if 1 /* plot all */
    for (i = 0; i < n; i++) {
	x = (pos[i][0]*ratio*0.5+0.5)*WIDTH;
	y = (pos[i][1]*ratio*0.5+0.5)*HEIGHT;
	gdImageSetPixel(im, x, y, yellow);
    }
#elif 1 /* plot -5.0<z<5.0 slice */
    for (i = 0; i < n; i++) {
	if (pos[i][2] < -5.0 || 5.0 < pos[i][2]) continue;
	x = (pos[i][0]*ratio*0.5+0.5)*WIDTH;
	y = (pos[i][1]*ratio*0.5+0.5)*HEIGHT;
	gdImageSetPixel(im, x, y, yellow);
    }
#else /* plot randomlt chosen 100000 particles */
    for (i = 0, ii = 0; i < n; i++) {
	if (i > 100000 && i%(n/100000) != 1)  continue;
	x = (pos[ii][0]*ratio*0.5+0.5)*WIDTH;
	y = (pos[ii][1]*ratio*0.5+0.5)*HEIGHT;
	gdImageSetPixel(im, x, y, yellow);
	ii++;
    }
#endif

    imgno++;
    fp = fopen(fname, "w");
    if (!fp) {
	fprintf(stderr, "fopen() failed for `%s'\n", fname);
	exit (2);
    }
    gdImageGif(im, fp);
    gdImageDestroy(im);
    fclose(fp);
}
#endif /* image quality */

#endif /* USEGD */

#if USEX11

void
vtc_viewlist(int nall, double (*ballpos)[3], /* all bodies */
	     int n, double (*bpos)[3], /* bodies in the interaction list */
	     int nc, double (*cpos)[3], double *csize, /* cells in the interaction list */
	     int ni) /* bodies in the acceleration box */
{
    static int firstcall = 1;
    static Polyh ballpolyh, bipolyh, bjpolyh, cipolyh, cjpolyh;
    static Snapshot s;
    static ImgContext imgcontext;

    will_quit = 0;
    if (firstcall) {
	firstcall = 0;
	initXContext(&imgcontext);
	imgcontext.showmode = CELLONLY;
    }
    s.nb = nall;
    s.bpos = ballpos;
    body2polyh(&s, &ballpolyh);
    s.nb = n-ni;
    s.bpos = bpos+ni;
    body2polyh(&s, &bjpolyh);
    s.nb = ni;
    s.bpos = bpos;
    body2polyh(&s, &bipolyh);
    s.nc = nc-1;
    s.cpos = cpos+1;
    s.csize = csize+1;
    treenode2polyh(&s, &cjpolyh);
    s.nc = 1;
    s.cpos = cpos;
    s.csize = csize;
    treenode2polyh(&s, &cipolyh);
    ilen = ni;
    jlen = n+nc-1;
    Cfprintf(stderr, "ni: %d list len: %d\n" ,ni, n+nc-1);
    exposeEvent1(&imgcontext, &ballpolyh, &bipolyh, &bjpolyh, &cipolyh, &cjpolyh);
    mainLoop1(&imgcontext, &ballpolyh, &bipolyh, &bjpolyh, &cipolyh, &cjpolyh);
    freepolyh(&ballpolyh);
    freepolyh(&bipolyh);
    freepolyh(&bjpolyh);
    freepolyh(&cipolyh);
    freepolyh(&cjpolyh);
}

void
vtc_viewtree(int n, double (*bpos)[3], int nc, double (*cpos)[3], double *csize)
{
    static int firstcall = 1;
    static Polyh bpolyh, cpolyh;
    static Snapshot s;
    static ImgContext imgcontext;

    will_quit = 0;
    if (firstcall) {
	firstcall = 0;
	initXContext(&imgcontext);
	imgcontext.showmode = CELLONLY;
    }
    s.nb = n;
    s.bpos = bpos;
    s.nc = nc;
    s.cpos = cpos;
    s.csize = csize;
    body2polyh(&s, &bpolyh);
    treenode2polyh(&s, &cpolyh);
    exposeEvent0(&imgcontext, &bpolyh, &cpolyh);
    mainLoop0(&imgcontext, &bpolyh, &cpolyh);
    freepolyh(&bpolyh);
    freepolyh(&cpolyh);
}

static void
initXContext(ImgContext *ic)
{
    double za = 0.2;

    ic->dpy = XOpenDisplay(NULL);
    ic->root=DefaultRootWindow(ic->dpy);
    ic->screen = DefaultScreen(ic->dpy);
    ic->white = WhitePixel (ic->dpy, ic->screen);
    ic->black = BlackPixel (ic->dpy, ic->screen);
    ic->w = XCreateSimpleWindow(ic->dpy, ic->root, 100, 100,
			    WIDTH, HEIGHT, BORDER,
			    ic->white, ic->black);
    ic->gcline = XCreateGC(ic->dpy, ic->w, 0, NULL);
    ic->gcline2 = XCreateGC(ic->dpy, ic->w, 0, NULL);
    ic->gcclear = XCreateGC(ic->dpy, ic->w, 0, NULL);
    ic->gcaxis  = XCreateGC(ic->dpy, ic->w, 0, NULL);
    XSetBackground(ic->dpy, ic->gcline, ic->black);
    XSetForeground(ic->dpy, ic->gcline, 0xffff00);
    XSetBackground(ic->dpy, ic->gcline2, ic->black);
    XSetForeground(ic->dpy, ic->gcline2, 0x0f0f00);
    XSetBackground(ic->dpy, ic->gcaxis, ic->black);
    XSetForeground(ic->dpy, ic->gcaxis, ic->white);
    XSetBackground(ic->dpy, ic->gcclear, ic->black);
    XSetForeground(ic->dpy, ic->gcclear, ic->black);
    XSetLineAttributes(ic->dpy, ic->gcaxis,
		       1, LineOnOffDash, 0, 0);
    XSetLineAttributes(ic->dpy, ic->gcline,
		       1, LineSolid, 0, 0);
    XSelectInput(ic->dpy, ic->w, ExposureMask|KeyPressMask);
    XMapWindow(ic->dpy, ic->w);
    XFlush(ic->dpy);

    /* init point of view */
    ic->theta = 0.0;
    ic->phi = 0.0;
    ic->zoom = 200.0;

    /* init axis object */
    ic->axis.n = 6;
    ic->axis.ns = 3;
    ic->axis.v3 = (double (*)[3])calloc(sizeof(double) * 3, 6);
    ic->axis.v2 = (int (*)[2])calloc(sizeof(int) * 2, 6);
    ic->axis.s = (int (*)[2])calloc(sizeof(int) * 2, 3);
    ic->axis.v3[0][0] = za;
    ic->axis.v3[0][1] = 0.0;
    ic->axis.v3[0][2] = 0.0;
    ic->axis.v3[1][0] = -za;
    ic->axis.v3[1][1] = 0.0;
    ic->axis.v3[1][2] = 0.0;
    ic->axis.v3[2][0] = 0.0;
    ic->axis.v3[2][1] = za;
    ic->axis.v3[2][2] = 0.0;
    ic->axis.v3[3][0] = 0.0;
    ic->axis.v3[3][1] = -za;
    ic->axis.v3[3][2] = 0.0;
    ic->axis.v3[4][0] = 0.0;
    ic->axis.v3[4][1] = 0.0;
    ic->axis.v3[4][2] = za;
    ic->axis.v3[5][0] = 0.0;
    ic->axis.v3[5][1] = 0.0;
    ic->axis.v3[5][2] = -za;
    ic->axis.s[0][0] = 0;
    ic->axis.s[0][1] = 1;
    ic->axis.s[1][0] = 2;
    ic->axis.s[1][1] = 3;
    ic->axis.s[2][0] = 4;
    ic->axis.s[2][1] = 5;
}

static void
body2polyh(Snapshot *t, Polyh *p)
{
    int n, i, k;
    double (*v)[3];
    int (*s)[2];

    p->n = n = t->nb;

    /* set v[n][3] */
    p->v3 = v = (double (*)[3])calloc(sizeof(double) * 3, n);
    p->v2 = (int (*)[2])calloc(sizeof(int) * 2, n);
    for (i = 0; i < n; i++) {
	for (k = 0; k < 3; k++) {
	    p->v3[i][k] = t->bpos[i][k];
	}
    }
    p->ns = 0;
}


static void
treenode2polyh(Snapshot *t, Polyh *p)
{
    int n, ns, i, j;
    double (*v)[3];
    int (*s)[2];

    p->n = n = t->nc*8;

    /* set v[n][3] */
    p->v3 = v = (double (*)[3])calloc(sizeof(double) * 3, n);
    p->v2 = (int (*)[2])calloc(sizeof(int) * 2, n);
    for (i = 0; i < n; i += 8)
    {
	double c[3];
	double w, w2;

	c[0] = t->cpos[i/8][0];
	c[1] = t->cpos[i/8][1];
	c[2] = t->cpos[i/8][2];
	w = t->csize[i/8];
	w2 = w/2.0;
/*
	fprintf(stderr, "w: %f w2: %f\n", w, w2);
	*/
	for (j = 0; j < 8; j++) {
	    v[i+j][0] = c[0]+w2*((j>>2)%2>0 ? 1.0:-1.0);
	    v[i+j][1] = c[1]+w2*((j>>1)%2>0 ? 1.0:-1.0);
	    v[i+j][2] = c[2]+w2*((j>>0)%2>0 ? 1.0:-1.0);
	}
    }
    /* set s[ns][2] */
    p->ns = ns = t->nc*12;
    p->s = s = (int (*)[2])calloc(sizeof(int) * 2, ns);
    for (i = 0, j = 0; i < ns; i += 12, j += 8) {
	s[i+0][0] = j+0;
	s[i+0][1] = j+1;
	s[i+1][0] = j+1;
	s[i+1][1] = j+3;
	s[i+2][0] = j+3;
	s[i+2][1] = j+2;
	s[i+3][0] = j+2;
	s[i+3][1] = j+0;

	s[i+4][0] = j+4;
	s[i+4][1] = j+5;
	s[i+5][0] = j+5;
	s[i+5][1] = j+7;
	s[i+6][0] = j+7;
	s[i+6][1] = j+6;
	s[i+7][0] = j+6;
	s[i+7][1] = j+4;

	s[i+8][0] = j+0;
	s[i+8][1] = j+4;
	s[i+9][0] = j+1;
	s[i+9][1] = j+5;
	s[i+10][0] = j+3;
	s[i+10][1] = j+7;
	s[i+11][0] = j+2;
	s[i+11][1] = j+6;
    }
}

static void
mainLoop0(ImgContext *ic, Polyh *bp, Polyh *cp)
{
    while (!will_quit)
    {
	XNextEvent(ic->dpy, &(ic->event));
	switch (ic->event.type)
	{
	case Expose:
	    exposeEvent0(ic, bp, cp);
	    break;

	case KeyPress:
	    keypressEvent(ic);
	    if (will_dump) {
		will_dump = 0;
#if USEGD
		dumpimage0(ic, bp, cp);
#endif /* USEGD */
	    }
	    exposeEvent0(ic, bp, cp);
	    break;

	default:
	    break;
	}
    }
}

static void
mainLoop1(ImgContext *ic, Polyh *ballp, Polyh *bip, Polyh *bjp, Polyh *cpi, Polyh *cpj)
{
    while (!will_quit)
    {
	XNextEvent(ic->dpy, &(ic->event));
	switch (ic->event.type)
	{
	case Expose:
	    exposeEvent1(ic, ballp, bip, bjp, cpi, cpj);
	    break;

	case KeyPress:
	    keypressEvent(ic);
	    if (will_dump) {
		will_dump = 0;
#if USEGD
		dumpimage1(ic, ballp, bip, bjp, cpi, cpj);
#endif /* USEGD */
	    }
	    exposeEvent1(ic, ballp, bip, bjp, cpi, cpj);
	    break;

	default:
	    break;
	}
    }
}

#if USEGD

static void
dumpimage0(ImgContext *ic, Polyh *bp, Polyh *cp)
{
    gdImagePtr im;
    int bg, fg0, fg1, fg2;
    FILE *fp;
    char fname[255];
    static int imgno = 0;

    im = gdImageCreate(WIDTH, HEIGHT);
    bg = gdImageColorAllocate(im, 0, 0, 0);
    fg0 = gdImageColorAllocate(im, 255, 255, 255);
    fg1 = gdImageColorAllocate(im, 255, 255, 0);
    fg2 = gdImageColorAllocate(im, 0, 255, 0);
    gdImageFill(im, WIDTH, HEIGHT, bg);

    dumpAxis(im, ic, &(ic->axis), fg0);
    ic->showmode = ic->showmode % 3;
    switch (ic->showmode) {
    case BODYONLY:
	dumpPolyh(im, ic, bp, fg1);
	break;
    case CELLONLY:
	dumpPolyh(im, ic, cp, fg1);
	break;
    case BODY_AND_CELL:
	dumpPolyh(im, ic, bp, fg1);
	dumpPolyh(im, ic, cp, fg1);
	break;
    default:
	fprintf(stderr, "unknown showmode: %d\n",
		ic->showmode);
	break;
    }

    sprintf(fname, "img%03d.gif", imgno);
    fp = fopen(fname, "w");
    if (!fp) {
	fprintf(stderr, "fopen() failed for `%s'\n", fname);
	exit (2);
    }
    gdImageGif(im, fp);
    gdImageDestroy(im);
    fclose(fp);
    imgno++;
}

static void
dumpimage1(ImgContext *ic, Polyh *ballp, Polyh *bip, Polyh *bjp,
	   Polyh *cip, Polyh *cjp)
{
    gdImagePtr im;
    int bg, fg0, fg1, fg2;
    FILE *fp;
    char fname[255];
    char buf[128];
    static int imgno = 0;

    im = gdImageCreate(WIDTH, HEIGHT);
    bg = gdImageColorAllocate(im, 0, 0, 0);
    fg0 = gdImageColorAllocate(im, 255, 255, 255);
    fg1 = gdImageColorAllocate(im, 255, 255, 0);
    fg2 = gdImageColorAllocate(im, 0, 255, 0);
    gdImageFill(im, WIDTH, HEIGHT, bg);

    sprintf(buf, "i-particle: %d", (int)ilen);
    gdImageString(im,
		  gdFontMediumBold,
		  WIDTH*3.0/4.0, HEIGHT-45,
		  buf, fg0);
    sprintf(buf, "j-particle: %d", (int)jlen);
    gdImageString(im,
		  gdFontMediumBold,
		  WIDTH*3.0/4.0, HEIGHT-30,
		  buf, fg0);

    projectPolyh(ic, &(ic->axis), 250);
    projectPolyh(ic, ballp, ic->zoom);
    projectPolyh(ic, bip, ic->zoom);
    projectPolyh(ic, bjp, ic->zoom);
    projectPolyh(ic, cip, ic->zoom);
    projectPolyh(ic, cjp, ic->zoom);
    
    dumpAxis(im, ic, &(ic->axis), fg0);
    ic->showmode = ic->showmode % 4;
    switch (ic->showmode) {
    case LISTBODYONLY:
	dumpPolyh(im, ic, bjp, fg1);
	dumpPolyh(im, ic, cip, fg2);
	dumpPolyh(im, ic, bip, fg2);
	break;
    case LISTONLY:
	dumpPolyh(im, ic, cjp, fg1);
	dumpPolyh(im, ic, bjp, fg1);
	dumpPolyh(im, ic, cip, fg2);
	dumpPolyh(im, ic, bip, fg2);
	break;
    case LIST_AND_BODY:
	dumpPolyh(im, ic, cjp, fg1);
	dumpPolyh(im, ic, ballp, fg1);
	dumpPolyh(im, ic, cip, fg2);
	dumpPolyh(im, ic, bip, fg2);
	break;
    case BODYONLY2:
	dumpPolyh(im, ic, ballp, fg1);
	dumpPolyh(im, ic, bip, fg2);
	break;
    default:
	fprintf(stderr, "unknown showmode: %d\n",
		ic->showmode);
	break;
    }

    sprintf(fname, "img%03d.gif", imgno);
    fp = fopen(fname, "w");
    if (!fp) {
	fprintf(stderr, "fopen() failed for `%s'\n", fname);
	exit (2);
    }
    gdImageGif(im, fp);
    gdImageDestroy(im);
    fclose(fp);
    imgno++;
}

#endif /* USEGD */

static void
exposeEvent0(ImgContext *ic, Polyh *bp, Polyh *cp)
{
    clearDrawable(ic);

    projectPolyh(ic, &(ic->axis), WIDTH*0.5);
    projectPolyh(ic, bp, ic->zoom);
    projectPolyh(ic, cp, ic->zoom);
    
    displayAxis(ic, &(ic->axis), ic->gcaxis);
    ic->showmode = ic->showmode % 3;
    switch (ic->showmode) {
    case BODYONLY:
	displayPolyh(ic, bp, ic->gcline);
	break;
    case CELLONLY:
	displayPolyh(ic, cp, ic->gcline);
	break;
    case BODY_AND_CELL:
	displayPolyh(ic, bp, ic->gcline);
	displayPolyh(ic, cp, ic->gcline);
	break;
    default:
	fprintf(stderr, "unknown showmode: %d\n",
		ic->showmode);
	break;
    }
}

static void
exposeEvent1(ImgContext *ic, Polyh *ballp, Polyh *bip, Polyh *bjp,
	     Polyh *cip, Polyh *cjp)
{
    char buf[128];

    clearDrawable(ic);
    sprintf(buf, "i-particle: %d", (int)ilen);
    XDrawString(ic->dpy, ic->w, ic->gcaxis,
		WIDTH*3.0/4.0, HEIGHT-45,
		buf, strlen(buf));
    sprintf(buf, "j-particle: %d", (int)jlen);
    XDrawString(ic->dpy, ic->w, ic->gcaxis,
		WIDTH*3.0/4.0, HEIGHT-30,
		buf, strlen(buf));

    projectPolyh(ic, &(ic->axis), 250);
    projectPolyh(ic, ballp, ic->zoom);
    projectPolyh(ic, bip, ic->zoom);
    projectPolyh(ic, bjp, ic->zoom);
    projectPolyh(ic, cip, ic->zoom);
    projectPolyh(ic, cjp, ic->zoom);
    
    displayAxis(ic, &(ic->axis), ic->gcaxis);
    ic->showmode = ic->showmode % 4;
    switch (ic->showmode) {
    case LISTBODYONLY:
	displayPolyh(ic, bjp, ic->gcline);
	displayPolyh(ic, cip, ic->gcline2);
	displayPolyh(ic, bip, ic->gcline2);
	break;
    case LISTONLY:
	displayPolyh(ic, cjp, ic->gcline);
	displayPolyh(ic, bjp, ic->gcline);
	displayPolyh(ic, cip, ic->gcline2);
	displayPolyh(ic, bip, ic->gcline2);
	break;
    case LIST_AND_BODY:
	displayPolyh(ic, cjp, ic->gcline);
	displayPolyh(ic, ballp, ic->gcline);
	displayPolyh(ic, cip, ic->gcline2);
	displayPolyh(ic, bip, ic->gcline2);
	break;
    case BODYONLY2:
	displayPolyh(ic, ballp, ic->gcline);
	displayPolyh(ic, bip, ic->gcline2);
	break;
    default:
	fprintf(stderr, "unknown showmode: %d\n",
		ic->showmode);
	break;
    }
}

static void
keypressEvent(ImgContext *ic)
{
    char key[20];
    KeySym keysym;
    XComposeStatus xcom;

    XLookupString((XKeyEvent *)&(ic->event), key, sizeof(key), &keysym, &xcom);

    switch (key[0])
    {
    case '+':
	ic->zoom *= 1.1;
	break;
    case '-':
	ic->zoom /= 1.1;
	break;
    case 'j':
	ic->phi +=6.0;
	ic->phi = (int)(ic->phi)%360;
	break;
    case 'k':
	ic->phi -=6.0;
	ic->phi = (int)(ic->phi)%360;
	break;
    case 'h':
	ic->theta +=6.0;
	ic->theta = (int)(ic->theta)%360;
	break;
    case 'l':
    case ' ':
	ic->theta -=6.0;
	ic->theta = (int)(ic->theta)%360;
	break;
    case 's':
	ic->showmode++;
	break;
    case 'd':
	will_dump = 1;
	break;
    case 'q':
	will_quit = 1;
	break;
    default:
	break;
    }
}

static void
projectPolyh(ImgContext *ic, Polyh *p, double zoom)
{
    int i;
    double x0, x1, x2;
    double sint = sin(M_PI * (ic->theta) / 180.0);
    double cost = cos(M_PI * (ic->theta) / 180.0);
    double sinp = sin(M_PI * (ic->phi) / 180.0);
    double cosp = cos(M_PI * (ic->phi) / 180.0);

    for (i = 0; i < p->n; i++) {
	x0 = p->v3[i][0];
	x1 = p->v3[i][1];
	x2 = p->v3[i][2];

	p->v2[i][0] = (x0 * cost - x2 * sint)
	    * zoom + WIDTH/2;
	p->v2[i][1] = -(x0 * sint * sinp + x1 * cosp + x2 * cost * sinp)
	    * zoom + HEIGHT/2;
    }
}

#if USEGD

static void
dumpAxis(gdImagePtr im, ImgContext *ic, Polyh *p, int c)
{
    int offx = WIDTH*0.3;
    int offy = -HEIGHT*0.4;
    int i;
    char buf[128];

    sprintf(buf, "zoom: %d", (int)ic->zoom);
    gdImageString(im,
		  gdFontMediumBold,
		  10, HEIGHT-45,
		  buf, c);
    sprintf(buf, "theta: %d", (int)ic->theta);
    gdImageString(im,
		  gdFontMediumBold,
		  10, HEIGHT-30,
		  buf, c);

    sprintf(buf, "phi: %d", (int)ic->phi);
    gdImageString(im,
		  gdFontMediumBold,
		  10, HEIGHT-15,
		  buf, c);

    gdImageString(im,
		  gdFontMediumBold,
		  ic->axis.v2[0][0]-offx, ic->axis.v2[0][1]-offy,
		  "x", c);
    gdImageString(im,
		  gdFontMediumBold,
		  ic->axis.v2[2][0]-offx, ic->axis.v2[2][1]-offy,
		  "y", c);
    gdImageString(im,
		  gdFontMediumBold,
		  ic->axis.v2[4][0]-offx, ic->axis.v2[4][1]-offy,
		  "z", c);
    for (i = 0; i < p->ns; i++)
    {
	gdImageLine(im,
		    p->v2[p->s[i][0]][0]-offx, p->v2[p->s[i][0]][1]-offy,
		    p->v2[p->s[i][1]][0]-offx, p->v2[p->s[i][1]][1]-offy, c);
    }
}

static void
dumpPolyh(gdImagePtr im, ImgContext *ic, Polyh *p, int c)
{
    int i;
    double zoom = ic->zoom;
    int z = 20.0*zoom/WIDTH+1.0;
    if (z < 2) {
	z = 2;
    }
    if (p->ns == 0) { /* particle data */
	for (i = 0; i < p->n; i++) {
#if 0
	    gdImageArc(im, p->v2[i][0]-z*0.5, p->v2[i][1]-z*0.5,
		       z, z,
		       0, 360, c);
	    gdImageFill(im, p->v2[i][0], p->v2[i][1], c);
#else
	    gdImageSetPixel(im, p->v2[i][0]-1, p->v2[i][1]-1, c);
	    gdImageSetPixel(im, p->v2[i][0]-0, p->v2[i][1]-1, c);
	    gdImageSetPixel(im, p->v2[i][0]+1, p->v2[i][1]-1, c);
	    gdImageSetPixel(im, p->v2[i][0]-1, p->v2[i][1]-0, c);
	    gdImageSetPixel(im, p->v2[i][0]-0, p->v2[i][1]-0, c);
	    gdImageSetPixel(im, p->v2[i][0]+1, p->v2[i][1]-0, c);
	    gdImageSetPixel(im, p->v2[i][0]-1, p->v2[i][1]+1, c);
	    gdImageSetPixel(im, p->v2[i][0]-0, p->v2[i][1]+1, c);
	    gdImageSetPixel(im, p->v2[i][0]+1, p->v2[i][1]+1, c);
	    if (c == GREEN) {
		int ii;
		for (ii = 0; ii < 4; ii++) {
		    gdImageSetPixel(im, p->v2[i][0]-2, p->v2[i][1]+ii-2, c);
		    gdImageSetPixel(im, p->v2[i][0]+2, p->v2[i][1]+ii-2, c);
		}
		for (ii = 0; ii < 3; ii++) {
		    gdImageSetPixel(im, p->v2[i][0]+ii-1, p->v2[i][1]+2, c);
		    gdImageSetPixel(im, p->v2[i][0]+ii-1, p->v2[i][1]-2, c);
		}
	    }

#endif
	}
    }
    else { /* polyhedron data */
	for (i = 0; i < p->ns; i++) {
	    gdImageLine(im,
			p->v2[p->s[i][0]][0], p->v2[p->s[i][0]][1],
			p->v2[p->s[i][1]][0], p->v2[p->s[i][1]][1],
			c);
	    if (c == GREEN) {
		gdImageLine(im,
			    p->v2[p->s[i][0]][0]+1, p->v2[p->s[i][0]][1]+1,
			    p->v2[p->s[i][1]][0]+1, p->v2[p->s[i][1]][1]+1,
			    c);
	    }
	}
    }
}

#endif /* USEGD */

static void
displayAxis(ImgContext *ic, Polyh *p, GC gc)
{
    int offx = WIDTH*0.3;
    int offy = -HEIGHT*0.4;
    int i;
    char buf[128];

    sprintf(buf, "zoom: %d", (int)ic->zoom);
    XDrawString(ic->dpy, ic->w, ic->gcaxis,
		10, HEIGHT-45,
		buf, strlen(buf));

    sprintf(buf, "theta: %d", (int)ic->theta);
    XDrawString(ic->dpy, ic->w, ic->gcaxis,
		10, HEIGHT-30,
		buf, strlen(buf));

    sprintf(buf, "phi: %d", (int)ic->phi);
    XDrawString(ic->dpy, ic->w, ic->gcaxis,
		10, HEIGHT-15,
		buf, strlen(buf));

    XDrawString(ic->dpy, ic->w, ic->gcaxis,
		ic->axis.v2[0][0]-offx, ic->axis.v2[0][1]-offy,
		"x", strlen("x"));
    XDrawString(ic->dpy, ic->w, ic->gcaxis,
		ic->axis.v2[2][0]-offx, ic->axis.v2[2][1]-offy,
		"y", strlen("y"));
    XDrawString(ic->dpy, ic->w, ic->gcaxis,
		ic->axis.v2[4][0]-offx, ic->axis.v2[4][1]-offy,
		"z", strlen("z"));

    for (i = 0; i < p->ns; i++)
    {
	XDrawLine(ic->dpy, ic->w, gc,
		  p->v2[p->s[i][0]][0]-offx, p->v2[p->s[i][0]][1]-offy,
		  p->v2[p->s[i][1]][0]-offx, p->v2[p->s[i][1]][1]-offy);
    }
}

static void
displayPolyh(ImgContext *ic, Polyh *p, GC gc)
{
    int i;
    double zoom = ic->zoom;
    int z = 20.0*zoom/WIDTH+1.0;
    if (z < 2) {
	z = 2;
    }
    if (p->ns == 0) { /* particle data */
	for (i = 0; i < p->n; i++) {
	    XFillArc(ic->dpy, ic->w, gc,
		     p->v2[i][0]-z*0.5, p->v2[i][1]-z*0.5,
		     z, z, /* radius */
		     0*64, 360*64);
	}
    }
    else { /* polyhedron data */
	for (i = 0; i < p->ns; i++) {
	    XDrawLine(ic->dpy, ic->w, gc,
		      p->v2[p->s[i][0]][0], p->v2[p->s[i][0]][1],
		      p->v2[p->s[i][1]][0], p->v2[p->s[i][1]][1]);
	}
    }
}

static void
clearDrawable(ImgContext *ic)
{
    XFillRectangle(ic->dpy, ic->w, ic->gcclear,
		   0, 0,
		   WIDTH, HEIGHT);
}

static void
freepolyh(Polyh *p)
{
    free(p->v3);
    free(p->v2);
    if (p->ns != 0) {
	free(p->s);
    }
}
#endif /* USEX11 */
