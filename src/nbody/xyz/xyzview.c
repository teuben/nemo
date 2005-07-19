/*
 * XYZVIEW.C: interactively display a file of x,y,z,color points.
 *
 *      xx-xxx-xx  V1.0  Created                        Joshua Barnes
 *       8-mar-92  V1.9  some version                             JEB
 *      21-jan-93  V2.0  VOGL implementation (tested: X windows)  PJT
 *                    a  some mods to get it compile on our SGI   
 *			 added movie features 			  PJT
 *	30-mar-97  V2.1  readied for NEMO export, optional -DOLDVOGL
 *	 4-apr-97     a  FOURKEY is also a movie toggle		  pjt
 *	22-aug-00     b  ansi cc				  pjt
 *	 4-sep-00  V2.2  select some point always visible	  pjt
 *                       with orbit=
 *      19-may-00  V2.3  if input is SnapShot, autoconvert to PointData
 *      10-dec-03  V3.0  added movie=t|f option
 *
 *
 *  Bugs: VOGL has an event stack which is only one deep, reason
 *        to have to hit a key twice to get the requested action done
 *  Todo: when switching modes in VOGL, using the 1,2,3 keys, it would
 *        be good to remember the cursor position and restore it when
 * 	  you get back to that mode - this prevent annoying jumping
 *	  around in rotation angles or zoom depth (1,2) in particular
 *	  i.e. not to force an absolute cursor positioning?
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		 /* NEMO  -- now also includes ZENO's vectdefs.h */
#include <filestruct.h>

#if defined(OLD_VOGL)
# include <vogl.h> 		 /* probably need to add -I$(VOGLINC) or so */
# include <vodevice.h>
#elif defined(MESA)
# include <GL/gl.h>
#else
# include <gl.h>                 /* if VOGL: add: -I$NEMOINC/vogl to CFLAGS */
# include <device.h>
#endif

#if 0
#ifndef VOGL
# include <fmclient.h>           /* VOGL doesn't have font manager (-lfm) */
#endif
#endif

#ifdef VOGL
# ifndef XMAXSCREEN 
#  define XMAXSCREEN  511          /* 1279 on SGI */
# endif
# ifndef YMAXSCREEN 
#  define YMAXSCREEN  511         /* 1023 on SGI */
# endif
#endif

#define XMAXSCREEN  768
#define YMAXSCREEN  768

string defv[] = {
    "in=???\n			Point data to display",
    "times=all\n		Range of times to display",
    "scale=2.0\n		Scale factor for data",
    "nfast=2048\n		Max pnts to draw in fast mode",
    "showbox=true\n		Draw a cube, scale units high",
    "maxpoint=1024\n		Process at least this many",
    "colormap=\n		File with color table (binary)",
    "fullscreen=false\n		Use full screen (no border)",
    "position=\n		If given, set screen position",
    "aspect=\n			If given, set aspect ratio",
    "noborder=false\n		Disable drawing window border",
    "textsize=20\n		Size of text in points",
    "ident=\n			ID message for display",
    "viewfile=view.dat\n	Output viewing parameters",
    "maxframe=500\n             Maximum frames for movie storage",
    "orbit=\n                   Select this point/star always visible",
    "movie=f\n                  If true, start loading all frames",
    "VERSION=3.0a\n		17-jul-05 PJT",
    NULL,
};

#ifndef XYZVEL
string usage = "Display 3-D point data";
#else
string usage = "Display 3-D velocity data";
#endif

/* local routines : */

bool get_xyz(), handle_que();

extern string *burststring();
extern int xstrlen();

static bool Qmovie;

nemo_main()
{
    Qmovie = getbparam("movie");
    init_display();
    if (hasvalue("colormap")) setcolors(getparam("colormap"));
    if (!get_xyz(TRUE)) error("no data in file");
    plot_points();
    while (handle_que()) plot_points();
}

init_display()
{
    string *posstrs, *aspstrs;

    prefsize(512L,512L);	/* VOGL as well as GL */

    if (getbparam("fullscreen")) {
	if (hasvalue("position") || hasvalue ("aspect"))
	    error("fullscreen precludes position or aspect");
	prefposition(0, XMAXSCREEN, 0, YMAXSCREEN);
    }
    if (hasvalue("position")) {
	if (hasvalue("aspect"))
	    error("position precludes aspect");
	posstrs = burststring(getparam("position"), ",:");
	if (xstrlen(posstrs, sizeof(string)) != 5)	/* 5 includes NULL  */
	    error("position must have form xmin,xmax,ymin,ymax");
	prefposition(atoi(posstrs[0]), atoi(posstrs[1]),
		     atoi(posstrs[2]), atoi(posstrs[3]));
    }
    else if (hasvalue("aspect")) {
	aspstrs = burststring(getparam("aspect"), ",:");
	if (xstrlen(aspstrs, sizeof(string)) != 3)
	    error("%s: aspect must have form xratio,yratio\n",
		  getargv0());
	keepaspect(atoi(aspstrs[0]), atoi(aspstrs[1]));
    }
    if (getbparam("noborder") || getbparam("fullscreen"))
	noborder();
    foreground();
    winopen(getargv0());
    doublebuffer();
    gconfig();
#ifdef VOGL
    qdevice(ESCKEY);					/* exit program	    */
    qdevice(SPACEKEY);					/* next frame	    */
    qdevice(RETKEY);					/* save view params */
    qdevice(REDRAW);
    qdevice(MOUSEX);
    qdevice(MOUSEY);
    qdevice(LEFTMOUSE);
    qdevice(MIDDLEMOUSE);
    qdevice(RIGHTMOUSE);
    qdevice(HKEY);  /* help */
    qdevice(QKEY);  /* quit */
    qdevice(DKEY);  /* debug */
    qdevice(JKEY);  /* previous frame */
    qdevice(KKEY);  /* next frame */
    qdevice(LKEY);  /* load all frames */
    qdevice(MKEY);  /* movie toggle */
    qdevice(BKEY);  /* toggle viewing the box */
    qdevice(ZEROKEY);    /* as if deactivated all buttons */
    qdevice(ONEKEY);     /* as if activated left button */
    qdevice(TWOKEY);     /* as if activated middle button */
    qdevice(THREEKEY);   /* as if activated right button */
    qdevice(FOURKEY);    /* also movie toggle */
#endif
}

#define MAXCOL  4096

int ncolors = 8;	/* default SGI white color */

setcolors(map)
string map;
{
    stream cstr;
    int i, ir, ig, ib;
    real red[MAXCOL], green[MAXCOL], blue[MAXCOL];

    cstr = stropen(map, "r");
    get_history(cstr);
    get_data(cstr, "Ncolors", IntType, &ncolors, 0);
    get_data_coerced(cstr, "Red", RealType, red, ncolors, 0);
    get_data_coerced(cstr, "Green", RealType, green, ncolors, 0);
    get_data_coerced(cstr, "Blue", RealType, blue, ncolors, 0);
    strclose(cstr);
    for (i = 0; i < ncolors; i++) {
	ir = 255 * MAX(0.0, MIN(red[i], 1.0));
	ig = 255 * MAX(0.0, MIN(green[i], 1.0));
	ib = 255 * MAX(0.0, MIN(blue[i], 1.0));
	mapcolor(i, ir, ig, ib);
    }
}

stream instr = NULL;

real tpoint = 0.0;

int npoint = 0, maxpoint = 0;
int    *ctab = NULL;
vector *ptab = NULL;
vector *vtab = NULL;

typedef struct frame {
    int npoint;
    real tpoint;
    int *ctab;
    vector *ptab;
    vector *vtab;
} frame;


int maxframe=0;
int nframe=0, iframe=0;
int orbit;
frame *saved;


#define TIMEFUZZ  0.001

bool get_xyz(save)
bool save;          /* save it too ? */
{
    bool done;
    string times;
    char cmd[128];

    if (instr == NULL) {         /* first time around, open file */
	instr = stropen(getparam("in"), "r");
	get_history(instr);
	if (get_tag_ok(instr,"SnapShot")) {     /* ieck, convert via pipe */
            warning("stropen CANNOT YET read from popen");
	    strclose(instr);            /* note this program never strclose's */
            sprintf(cmd,"snapxyz %s -",getparam("in"));
            instr = stropen(cmd,"rp");
	}

        if (hasvalue("orbit"))
            orbit = getiparam("orbit");
        else
            orbit = -1;
    }
    done = FALSE;
    times = getparam("times");
    while (! done && get_tag_ok(instr, "PointData")) {
	get_set(instr, "PointData");
	get_data_coerced(instr, "Tpoint", RealType, &tpoint, 0);
	get_data(instr, "Npoint", IntType, &npoint, 0);
	if (maxpoint == 0)
	    maxpoint = MAX(npoint, getiparam("maxpoint"));
	if (npoint > maxpoint)
	    error("npoint = %d exceeds maxpoint", npoint);
	if (get_tag_ok(instr, "Color")) {
	    if (ctab == NULL) {
                dprintf(0,"[Allocating ctab for %d particles]\n",maxpoint);
		ctab = (int *) allocate(maxpoint * sizeof(int));
            }
	    get_data(instr, "Color", IntType, ctab, npoint, 0);
	}
	if (streq(times, "all") || within(tpoint, times, TIMEFUZZ)) {
	    if (ptab == NULL) {
                dprintf(0,"[Allocating ptab for %d particles]\n",maxpoint);
		ptab = (vector *) allocate(maxpoint * sizeof(vector));
            }
	    get_data_coerced(instr, "Position", RealType, ptab, npoint, 3, 0);
#ifdef XYZVEL
	    if (vtab == NULL) {
                dprintf(0,"[Allocating vtab for %d particles]\n",maxpoint);
		vtab = (vector *) allocate(maxpoint * sizeof(vector));
            }
	    get_data_coerced(instr, "Velocity", RealType, vtab, npoint, 3, 0);
#endif
	    done = TRUE;
	    dprintf(0, "[get_xyz: time = %.2f]\n", tpoint);
	}
	get_tes(instr, "PointData");
    }
    if (done && ctab == NULL) error("no color data provided");

    if (save && done) {
        if (maxframe==0) {
            maxframe = getiparam("maxframe");
            saved = (frame *) allocate(maxframe*sizeof(frame));
        }
        if (nframe < maxframe) {
           dprintf(0,"Saving %d points in frame %d\n",npoint,nframe+1);

           saved[nframe].npoint = npoint;
           saved[nframe].tpoint = tpoint;

            /* WARNING: should do npoint+maxframe or reallocate */

           saved[nframe].ctab = (int *) allocate(maxpoint * sizeof(int));
           memcpy(saved[nframe].ctab,ctab,npoint*sizeof(int));

           saved[nframe].ptab = (vector *) allocate(maxpoint * sizeof(vector));
           memcpy(saved[nframe].ptab,ptab,npoint*sizeof(vector));
#ifdef XYZVEL
           saved[nframe].vtab = (vector *) allocate(maxpoint * sizeof(vector));
           memcpy(saved[nframe].vtab,vtab,npoint*sizeof(vector));
#endif
           nframe++;
           iframe = nframe;
        } else
           warning("No more stored frames; increase maxframe=%d",maxframe);
    }    

    return done;
}

int  field;				/* field of view in 0.1 degrees	    */
real aspect;				/* aspect ratio of display	    */
real zview;				/* viewing distance		    */
int  xrot, yrot, zrot;			/* viewing angles in 0.1 degrees    */
real boxscale;				/* size of cube to draw             */
int  pskip;				/* sampling factor for points	    */

#ifdef XYZVEL
 real vscale;				/* scaling for velocity vectors     */
 int  vskip;				/* sampling for velocity vectors    */
#endif

plot_points()
{
    int i;

    set_view();
    color(0);
    clear();
    color(ncolors - 1);				/* assume last is white     */
    if (boxscale > 0.0)
	drawbox();
    dprintf(1,"plot_points: npoint=%d\n",npoint);
    for (i = 0; i < npoint; i += pskip) {
	color(ctab[i]);
	pnt(ptab[i][0], ptab[i][1], ptab[i][2]);
    }
#ifdef XYZVEL
    for (i = 0; i < npoint; i += pskip*vskip) {
	color(ctab[i]);
	move(ptab[i][0], ptab[i][1], ptab[i][2]);
	draw(ptab[i][0] + vscale*vtab[i][0],
	     ptab[i][1] + vscale*vtab[i][1],
	     ptab[i][2] + vscale*vtab[i][2]);
    }
#endif
    end_view();
    if (pskip == 1 && getiparam("textsize") != 0)
	plottext();
    swapbuffers();
}

set_view()
{
    set_transform();
    pushmatrix();
    perspective(field, aspect, 0.01 * zview, 100.0 * zview);
    polarview(zview, xrot, yrot, zrot);
}

drawbox()
{
    pushmatrix();
    drawside();
    rotate(900, 'x');
    drawside();
    rotate(900, 'x');
    drawside();
    rotate(900, 'x');
    drawside();
    popmatrix();
}

drawside()
{
    pushmatrix();
    translate(0.0, 0.0, boxscale);
    rect(-boxscale, -boxscale, boxscale, boxscale);
    popmatrix();
}

end_view()
{
    popmatrix();
}

#ifdef VOGL
plottext()
{
    int tsize, tlen;
    long xmax, ymax;
    char buf[64];

    tsize = ABS(getiparam("textsize"));
    sprintf(buf, "%.3f", tpoint);
    tlen = strwidth(buf);
    getsize(&xmax, &ymax);
    color(ncolors-1);			/* text in white color */
    if (getiparam("textsize") > 0) {
	cmov2i(xmax - tlen - tsize/2, ymax - 3*tsize/2);
        charstr(buf);
	cmov2i(tsize, ymax - 3*tsize/2);
	charstr(getparam("ident"));
    } else {
        rotate(900,'z');
	cmov2i(xmax - 3*tsize/2, tlen + tsize/2);
	charstr(buf);
	cmov2i(xmax - 3*tsize/2, ymax - tsize);
	charstr(getparam("ident"));
        rotate(900,'z');
    }	
}
#else
fmfonthandle textfont = NULL;

plottext()
{
    int tsize, tlen;
    long xmax, ymax;
    char buf[64];

    tsize = ABS(getiparam("textsize"));
    if (textfont == NULL) {
	fminit();
	textfont = fmfindfont("Times-Roman");
	if (textfont == NULL) error("font not found");
	textfont = fmscalefont(textfont, (double) tsize);
	fmsetfont(textfont);
    }
    sprintf(buf, "%.2f", tpoint);
    tlen = fmgetstrwidth(textfont, buf);
    getsize(&xmax, &ymax);
    cpack(ncolors - 1);
    if (getiparam("textsize") > 0) {
	cmov2i(xmax - tlen - tsize/2, ymax - 3*tsize/2);
	fmprstr(buf);
	cmov2i(tsize, ymax - 3*tsize/2);
	fmprstr(getparam("ident"));
    } else {
	fmrotatepagematrix(-90.0);
	cmov2i(xmax - 3*tsize/2, tlen + tsize/2);
	fmprstr(buf);
	cmov2i(xmax - 3*tsize/2, ymax - tsize);
	fmprstr(getparam("ident"));
	fmrotatepagematrix( 90.0);
    }	
}
#endif

bool ldown = FALSE, mdown = FALSE, rdown = FALSE, showbox, movie=FALSE;
int xcur, ycur;
long xorig, yorig, xsize, ysize;

set_transform()
{
    permanent real zscale = 0.0;
    int nfast;
    bool fastdraw;

    if (zscale == 0.0) {			/* init first time around */
	zscale = getdparam("scale");
	getsize(&xsize, &ysize);
	field = 600;
	aspect = (real) xsize / (real) ysize;
	zview = 2.5 * zscale;
	xrot = 0;
	yrot = 0;
	zrot = 0;
#ifdef XYZVEL
	vscale = 0.0;
	vskip = 1;
#endif
	showbox = getbparam("showbox");
    }
    if (!movie) {

    if (ldown) {
	xrot = (3600 * (xcur - xorig - xsize/2)) / xsize;
	yrot = (2000 * (ycur - yorig - ysize/2)) / ysize;
    }
    if (mdown) {
	zview = (5.0*zscale * (xcur - xorig)) / xsize;
	field = (1200 * (ycur - yorig)) / ysize;
    }
#ifdef XYZVEL
    if (rdown) {
	vscale = 0.025 * (zscale * (xcur - xorig)) / xsize;
	vskip = 1 << (8 * MAX(ycur - yorig, 0)) / ysize;
    }
#endif
    fastdraw = ldown || mdown || rdown;
    if (fastdraw & showbox)
	boxscale = 0.5 * zscale;
    else
	boxscale = 0.0;
    if (fastdraw) {
	nfast = getiparam("nfast");
	pskip = npoint / MIN(nfast, npoint);
	if (pskip * nfast < npoint)
	    pskip++;
    } else
	pskip = 1;

    } else {  /* movie mode */
	if (mdown)
            iframe = ((xcur-xorig)*nframe)/xsize;
        if (iframe>=0 && iframe<nframe) loadframe();
    }
}

bool handle_que()
{
    bool newmouse, exitflag;
    short dev, val;

    if (Qmovie)
      loadkey();

    getorigin(&xorig, &yorig);      /* why do this again and again */
    getsize(&xsize, &ysize);        /* aren't they fixed ?? */
    do {
	newmouse = exitflag = FALSE;
	switch (dev = qread(&val)) {
          case QKEY:
	  case ESCKEY:
	    exitflag = TRUE;
	    qreset();
	    break;
	  case SPACEKEY:    	    spacekey();     break;
	  case RETKEY:              returnkey();    break;
	  case REDRAW:              reshapeviewport();  break;
	  case MOUSEX:		/* doesn't seem to work in VOGL ? */
	    xcur = val;
	    newmouse = ! (ldown || mdown || rdown);
	    break;
	  case MOUSEY:		/* doesn't seem to work in VOGL ? */
	    ycur = val;
	    newmouse = ! (ldown || mdown || rdown);
	    break;
	  case LEFTMOUSE:	    lbutton();  break;
	  case MIDDLEMOUSE:	    mbutton();  break;
	  case RIGHTMOUSE:          rbutton();  break;
          case HKEY:      help();             break;
          case DKEY:      debug();            break;
          case JKEY:      previouskey();      break;
          case KKEY:      nextkey();          break;
          case LKEY:      loadkey();          break;
          case FOURKEY:
          case MKEY:      moviekey();         break;
          case BKEY:      boxkey();	      break;
          case ZEROKEY:   zerokey();          break;
          case ONEKEY:    leftkey();          break;
          case TWOKEY:    middlekey();        break;
          case THREEKEY:  rightkey();         break;
	  default:
	    break;
	}
    } while (qtest() || newmouse);
    mouse();    			/* VOGL only */
    return (! exitflag);
}

mouse()     /* VOGL only */
{   
#ifdef VOGL
    xcur = getvaluator(MOUSEX);
    ycur = getvaluator(MOUSEY);
#endif
}

bell(msg)
char *msg;
{
#ifdef VOGL
    printf("%s: mouse: ",msg);
    printf("%s", ldown ? "*" : "0");    
    printf("%s", mdown ? "*" : "0");    
    printf("%s", rdown ? "*" : "0");
    printf("\n");
    fflush(stdout);
#endif
}

debug()	    /* VOGL only */
{
    permanent bool keydown = FALSE;
    keydown = !keydown;
    if (!keydown) return;

    printf("________________________________________________________\n");
    printf("Program: %s\n",getargv0());
    printf("Mouse: xcur = %d ycur=%d\n",xcur,ycur);    
    printf("       ldown = %d, mdown = %d, rdown = %d\n",ldown,mdown,rdown);
    printf("Field of view = %f degrees\n",field*0.1);
    printf("Aspect ratio  = %f\n",aspect);
    printf("Viewing distance = %f\n",zview);
    printf("Viewing angles   = %f %f %f\n",xrot*0.1, yrot*0.1, zrot*0.1);
    printf("Size of cube     = %f\n",boxscale);
    printf("Sampling factor  = %d\n",pskip);
#ifdef XYZVEL
    printf("Vel scaling      = %f\n",vscale);
    printf("Vel sampling     = %d\n",vskip);
#endif
    printf("________________________________________________________\n");
}

zerokey()
{
    ldown = FALSE;   mdown=FALSE;    rdown=FALSE;
}
leftkey()
{
    ldown = TRUE;    mdown=FALSE;    rdown=FALSE;
}
middlekey()
{
    ldown = FALSE;   mdown=TRUE;     rdown=FALSE;
}
rightkey()
{
    ldown = FALSE;   mdown=FALSE;    rdown=TRUE;
}

lbutton()
{
    permanent int xlast = 0, ylast = 0;

    mouse();
    ldown = ! ldown;
    if (movie) {
	if (iframe>0) iframe--;
        return;
    }
    if (ldown) {
	xcur = xlast + xorig + xsize/2;
	ycur = ylast + yorig + ysize/2;
	setvaluator(MOUSEX, xcur, 0, XMAXSCREEN);
	setvaluator(MOUSEY, ycur, 0, YMAXSCREEN);
    } else {
	xlast = xcur - xorig - xsize/2;
	ylast = ycur - yorig - ysize/2;
    }
}

mbutton()
{
    permanent int xlast = 0, ylast = 0;
    
    mouse();
    mdown = ! mdown;
    if (mdown) {
	xcur = xlast + xorig + xsize/2;
	ycur = ylast + yorig + ysize/2;
	setvaluator(MOUSEX, xcur, 0, XMAXSCREEN);
	setvaluator(MOUSEY, ycur, 0, YMAXSCREEN);
    } else {
	xlast = xcur - xorig - xsize/2;
	ylast = ycur - yorig - ysize/2;
    }
}

rbutton()
{
    permanent int xlast = 0, ylast = 0;

    mouse();
    rdown = ! rdown;
    if (movie) {
	if (iframe<nframe-1) iframe++;
        return;
    }
#ifdef XYZVEL
    if (rdown) {
	xcur = xlast + xorig;
	ycur = ylast + yorig;
	setvaluator(MOUSEX, xcur, 0, XMAXSCREEN);
	setvaluator(MOUSEY, ycur, 0, YMAXSCREEN);
    } else {
	xlast = xcur - xorig;
	ylast = ycur - yorig;
    }
#endif
}

returnkey()
{
    permanent bool retdown = FALSE;
    permanent stream vstr = NULL;
    permanent int framecnt = 0;
    Matrix vmat;

    retdown = ! retdown;
    if (retdown) {
	if (vstr == NULL) {
	    vstr = stropen(getparam("viewfile"), "w!");
	    put_history(vstr);
	}
	set_view();
	getmatrix(vmat);
	end_view();
	put_set(vstr, "XyzView");
	put_data(vstr, "Xsize", IntType, &xsize, 0);
	put_data(vstr, "Ysize", IntType, &ysize, 0);
	put_data(vstr, "Tview", RealType, &tpoint, 0);
	put_data(vstr, "ViewMatrix", RealType, vmat, 4, 4, 0);
	put_tes(vstr, "XyzView");
	fflush(vstr);
	dprintf(0, "[%s: view %d at time %.3f]\n",
		getargv0(), ++framecnt, tpoint);
    }
}

spacekey()
{
    permanent bool spacedown = FALSE;

    spacedown = ! spacedown;
    if (spacedown)
	if (! get_xyz(TRUE))
	    ringbell();
}

loadkey() {
    int i, j;

    while (get_xyz(TRUE))
        plot_points();
    if (orbit >= 0 && iframe > 0) {
        dprintf(0,"Loading orbit for %d extra particles (%d)\n",iframe,orbit);
        if (orbit >= npoint) error("no such particle number %d",orbit);


        if (maxpoint < npoint + iframe + 1) {
            error("not enough space to hold buffers for extra orbit: %d < %d",
                maxpoint, npoint+iframe+1);
        }

        for (i=0; i<iframe; i++) {
            saved[i].npoint = npoint + iframe;
            for (j=0; j<iframe; j++) {
                saved[i].ctab[npoint+j] =  saved[j].ctab[orbit];
                saved[i].ctab[npoint+j] =  2;
                SETV(saved[i].ptab[npoint+j], saved[j].ptab[orbit]);
#ifdef XYZVEL
                SETV(saved[i].vtab[npoint+j], saved[j].vtab[orbit]);
#endif
            }
        }
        npoint += iframe;
        dprintf(0,"npoint now %d\n",npoint);
    }
    iframe = 0;
    loadframe();
}

nextkey()
{
    permanent bool spacedown = FALSE;

    spacedown = ! spacedown;
    if (spacedown) {
        if (iframe < nframe) {
            iframe++;
            loadframe();
        } else if (! get_xyz(TRUE))
	    ringbell();
    }
}

previouskey()
{
    permanent bool spacedown = FALSE;

    spacedown = ! spacedown;
    if (spacedown) {
        if (iframe>0) {
            iframe--;
            loadframe();
        }
    }

}

loadframe()
{
    permanent int lastframe = -1;

    if (lastframe == iframe) return;
    dprintf(1,"Loading frame %d\n",iframe);
    if (iframe < 0 || iframe >=  nframe) return;
    npoint = saved[iframe].npoint;
    tpoint = saved[iframe].tpoint;
    ctab = saved[iframe].ctab;
    ptab = saved[iframe].ptab;
#ifdef XYZVEL
    vtab = saved[iframe].vtab;
#endif
    lastframe = iframe;
}

moviekey()
{
    permanent bool keydown = FALSE;
    keydown = !keydown;
    if (keydown) {
	movie = !movie;
        if (nframe<2) {
            warning("You have not loaded enough frames; use K or L");
            movie = FALSE;
            return;
        }
        printf("MOVIE-mode (nframe=%d) is now: %s\n", nframe, 
                                              movie ? "on" : "off");
    }
}

boxkey()
{
    permanent bool keydown = FALSE;
    keydown = !keydown;
    if (keydown) {
	/* figure this out later - toggle box display */
    } 
}

#define p printf
help()
{
    permanent bool keydown = FALSE;

    keydown = !keydown;
    if (!keydown) return;

    p("ESC          quit                Q       quit\n");
    p("SPACE        next frame          H       help\n");
    p("RETURN       write viewfile      D       debug\n");
    p("LEFT MOUSE   (xrot,yrot)         J       get next frame and store\n");
    p("MIDDLE MOUSE (zview,field)       K       previous frame\n");
    p("RIGHT MOUSE  (vscale,vskip)      L       load all next frames\n");
    p("                                 M       toggle movie mode\n");
    p("In movie mode: LEFT/MIDDLE/MOUSE control the viewed frame\n");
    
}


/* stubs for things VOGL doesn't have */
#if defined(VOGL)
ringbell() 
{
warning("At end of data, use ESC to quit");
}

extern struct vdev     vdevice;

noborder() {}

setvaluator(dev,j,k,l)
Device dev;
int j,k,l;
{
    /* dummy stub: VOGL doesn't have it */
}


				/* some old VOGL's don't have this ? */
				/* 1.3.0 is ok, VOGL5D is ok too */
#ifdef OLD_VOGL

getsize(x,y)
int *x, *y;
{
    *x = XMAXSCREEN+1;
    *y = YMAXSCREEN+1;
}

getorigin(x,y)
int *x, *y;
{
    *x = 0;
    *y = 0;
}
#endif

#endif
