/*
 * YAPP: Yet Another Plotting Package.
 * Joshua Barnes  Sep 1986  I. A. S.  Princeton, NJ.
 * This version works with the SunCore graphics library.
 *
 * Mods by Peter Teuben 
 *    11-nov-86  open mouse device for interrogation (see e.g. snapplot)
 *		 automatic detection window/fullscreen
 *		 still nu auto-BW/CG detection
 *		 switch between temporary/retained segments 
 *		 allow old (retained) segments to be saved = overwrite mode
 *    17-Jun-87  new hacking PJT
 */

#include <usercore.h>
#include <stdinc.h>

#define BWScreen     0   	/*  0=color    1=BW          */
#define RETAINED     1          /*  0=use temporary segments  1=retained  */
	/*  temporary are much faster and can still be redrawn/scaled 
	    via mouse device in window mode
	 */
	 
#if BWScreen
int bw2dd();	
local struct vwsurf rawsurf = DEFAULT_VWSURF(bw2dd);	 /* BW full screen */
int pixwindd();	
local struct vwsurf winsurf = DEFAULT_VWSURF(pixwindd);  /* BW window      */
#else
int cg2dd();	
local struct vwsurf rawsurf = DEFAULT_VWSURF(cg2dd);	 /* CG full screen */
int cgpixwindd();	
local struct vwsurf winsurf = DEFAULT_VWSURF(cgpixwindd);/* CG window      */
#endif


static double dxymax;	/* size of user window */

static struct vwsurf *surface;	     /* this will be the surface we're using */

#if RETAINED
static int  keep_old_segments=0;      /*  1=TRUE, keep'm   0=FALSE, delete'm */
static int  segment_name=0;           /* name of current segment,  0=unused */
#endif

/*
 * PLINIT: initalize the plotting package.
 */

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;		/* output device name (ignored) */
double xmin, xmax;	/* user plotting area */
double ymin, ymax;
{
    float width, height;
    float gx1, gx2, gy1, gy2;
    float x1, x2, y1, y2;
    struct vwsurf *get_surface();
        
    printf ("TEST version PJT for YAPP.C\n");

    initialize_core(BASIC, SYNCHRONOUS, TWOD);

    surface = get_surface();
    initialize_view_surface(surface, FALSE);
    select_view_surface(surface);
/*
    inquire_ndc_space_2(&width,&height);
    printf ("NDC_space= %f %f\n",width,height);
    inquire_viewport_2(&gx1,&gx2,&gy1,&gy2);
    printf ("VIEWPORT = %f %f %f %f\n",gx1,gx2,gy1,gy2);
    inquire_window(&x1,&x2,&y1,&y2);
    printf ("WINDOW = %f %f %f %f \n",x1,x2,y1,y2);
 */
    if (ymax - ymin < xmax - xmin) {
	dxymax = xmax - xmin;
	width = 1.0;
	height = (ymax - ymin) / dxymax;
    } else {
	dxymax = ymax - ymin;
	width = (xmax - xmin) / dxymax;
	height = 1.0;
    }
/*
    printf ("NDC_space= %f %f\n",width,height);
 */
    set_ndc_space_2(width, height);
    set_viewport_2(0.0, width, 0.0, height);
    set_window(xmin, xmax, ymin, ymax);
    set_font(ROMAN);
    set_charprecision(CHARACTER);
    initialize_device(KEYBOARD, 1);
    initialize_device(BUTTON,1);
    initialize_device(BUTTON,2);
    initialize_device(BUTTON,3);        
    set_echo(KEYBOARD, 1, 0);
    set_echo_surface(BUTTON,1,surface);
    set_echo_surface(BUTTON,2,surface);
    set_echo_surface(BUTTON,3,surface);        
#if RETAINED
    create_retained_segment(++segment_name);
#else
    create_temporary_segment();
#endif
}

struct vwsurf *get_surface()	/*  returns a pointer to correct surface */
{
#if BWScreen
	printf ("BW SUN");
#else
	printf ("CG SUN");
#endif

	if (getenv("WINDOW_ME")) {
		printf (" in WINDOW mode\n");
		return (&winsurf);
	}
	else {
		printf (" in FULL SCREEN mode\n");
		return (&rawsurf);
	}

}		


/*
 * PLSWAP: does nothing.
 */

plswap() { }

/*
 * PLXSCALE, PLYSCALE: transform from usr to plotting coordinates.
 * At present, these do nothing (identity transformation).
 */

double plxscale(x, y)
double x, y;		/* user coordinates */
{
    return (x);
}

double plyscale(x, y)
double x, y;		/* user coordinates */
{
    return (y);
}

/*
 * PLLTYPE: does nothing.
 */

plltype(lwid, lpat)
int lwid;		/* line width */
int lpat;		/* line pattern */
{
    if (lwid > 0)
	set_linewidth(0.1 * (lwid - 1));
    if (lpat > 0)
	set_linestyle(lpat - 1);
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

plline(x, y)
double x, y;		/* user coordinates */
{
    line_abs_2(x, y);
}

plmove(x, y)
double x, y;		/* user coordinates */
{
    move_abs_2(x, y);
}

plpoint(x, y)
double x, y;		/* user coordinates */
{
    move_abs_2(x , y );
    line_abs_2(x , y );
}

/*
 * PLJUST: specify justification of strings and numbers.
 * Imports: jus: 
 */

static int textjust = -1;

pljust(jus)
int jus;		/* -1, 0, 1 for left, mid, right just */
{
    textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));
}

/*
 * PLTEXT: plot a text string.
 */

pltext(msg, x, y, hgt, ang)
char *msg;		/* message to plot, with NULL termination */
double x, y;		/* user coordinates (modified by justification) */
double hgt;		/* height of characters in user coordinates */
double ang;		/* angle of message, counterclockwise in degrees */
{
    double c, cos(), s, sin();
    float dx, dy;

    c = cos(ang / 57.296);
    s = sin(ang / 57.296);
    set_charpath_2(c, s);
    set_charup_2(- s, c);
    set_charsize(0.75 * hgt, hgt);
    inquire_text_extent_2(msg, &dx, &dy);
    move_abs_2(x - 0.55 * (1 + textjust) * dx,
               y - 0.55 * (1 + textjust) * dy);
				/* note fudge factor of 1.1 in offsets */
    text(msg);
}


/*
 * PLCIRCLE, PLCROSS, PLBOX: plot various symbols.
 */

plcircle(x, y, r)
double x, y;
double r;
{
    int npnts, i;
    double theta, sin(), cos();

    npnts = MAX(2400 * r / dxymax, 6.0);
    plmove(x + r, y);
    for (i = 1; i <= npnts; i++) {
	theta = TWO_PI * ((double) i) / ((double) npnts);
	plline(x + r * cos(theta), y + r * sin(theta));
    }
}

plcross(x, y, s)
double x, y;
double s;
{
    if (s > 0.0) {
	plmove(x - s, y);
	plline(x + s, y);
	plmove(x, y - s);
	plline(x, y + s);
    } else {
	s = s / 1.4142;
	plmove(x - s, y - s);
	plline(x + s, y + s);
	plmove(x - s, y + s);
	plline(x + s, y - s);
    }
}

plbox(x, y, s)
double x, y;
double s;
{
    if (s > 0.0) {
	plmove(x - s, y - s);
	plline(x + s, y - s);
	plline(x + s, y + s);
	plline(x - s, y + s);
	plline(x - s, y - s);
    } else {
	s = s * 1.4142;
	plmove(x - s, y);
	plline(x, y - s);
	plline(x + s, y);
	plline(x, y + s);
	plline(x - s, y);
    }
}

/*
 * PLFLUSH: output any pending graphics.
 */

plflush() { }

/*
 * PLFRAME: advance to next frame.
 */

plframe()
{

    mouse(1000*3600);	/* wait one hour or click mouse */

#if RETAINED
    if (keep_old_segments)
	    close_retained_segment();
    else
    	    delete_retained_segment(segment_name);
#else
    close_temporary_segment();
#endif

    new_frame();			/* clear screen and that stuff */

#if RETAINED
    create_retained_segment(++segment_name);
#else
    create_temporary_segment();
#endif
}

/*
 * PLSTOP: finish up a display.
 */

#define ONEMIN (60 * 1000 * 1000)

plstop()
{
    char msg[80];
    int len, button;

    button = mouse(1000*3600);

    terminate_device(BUTTON,1);
    terminate_device(BUTTON,2);
    terminate_device(BUTTON,3);    
#if RETAINED    
    close_retained_segment();
#else
    close_temporary_segment();
#endif        
/*    await_keyboard(ONEMIN, 1, msg, &len);   */

    

    deselect_view_surface(surface);
    terminate_core();
}


/*
 *  MOUSE: waits a maximum number of miliseconds and return either
 *	   no button touched (0) or the button (1,2,3) touched
 *         BUG: it seems that the documentation is wrong here:
 *	   a 'micro'second seems to be more like 1/4 'mili'second
 */

mouse (msec)
int msec;
{
	int k=0;
        int bell=7;
        
/*	printf ("MOUSE waiting %d\n",msec);  */

/*	set_echo(KEYBOARD,1,1);  	*/
/*	printf("De Bell Volgt:%c\n",bell);  */
	printf ("Enter any button to continue:");
	putchar(bell);
	putchar('\n');		/* send a line feed to flush the buffer */
/*	set_echo(KEYBOARD,1,0);   */
	while ( (k==0) && (msec-- != 0) ) 
		await_any_button (1,&k);   /* here may be a timing problem  */

/*	printf ("MOUSE OUT: k= %d msec=%d\n",k,msec);  */
	
	return(k);
}



#ifdef TOOL1

main()
{
	plinit ("***",0,20,0,20);
	
	plstop();
}

#endif
