/*
 * YAPPCGI: Yet Antoher Plotting Package, using CGI interface to SUN
 * Peter Teuben   Nov 86   IAS Princeton, NJ
 *
 */
 
#include <cgidefs.h>		/*  definitions of ENUM types etc. */
#include <stdio.h>

#define BWScreen 1		/*  1 for BW, 0 for Color */



/*	Global Static constants, used by the plotting package */

Cvwsurf 	device;
Cint		name;

double		ax,bx,ay,by;		/* transformation parameters */
double		xp,yp;			/* current pen position	 */

Ccoor		old,new;
Ccoor		list[2];
Ccoorlist	line;			/* used in plotting lines */

int		gxmin=0,		/* initial default VDC space */
		gymin=0,
		gxmax=32767,
		gymax=32676;

/*==========================================================================*/
/*
 * PLINIT
 */
 
plinit(pltdev, xmin, xmax, ymin, ymax)
char	*pltdev;		/* output device name (ignored) */
double	xmin, xmax;	/* user plotting area */
double	ymin, ymax;
{
	Ccoor	*c1,*c2;
	
	printf ("YAPPCGI \n");
				/******   initialize screen   *****/
	open_cgi();
	if (getenv("WINDOW_ME")) {
#if BWScreen
	NORMAL_VWSURF (device, PIXWINDD);	/* BW */
#else
	NORMAL_VWSURF (device, CGPIXWINDD);	/* color */
#endif
	}
	else {
#if BWScreen
	NORMAL_VWSURF (device, BW2DD);		/* BW */
#else
	NORMAL_VWSURF (device, CG2DD);		/* color */
#endif
	}
	open_vws(&name,&device);

				/***** set scaling properties *****/
				/*  perhaps use a window() later  */
	ax = (gxmax-gxmin)/(xmax-xmin);
	bx = gxmin - ax*xmin;
	ay = (gymax-gymin)/(ymax-ymin);
	by = gymin - ay*ymin;

	
				/**** some other useful things *****/
	text_precision(CHARACTER);
	text_font_index(ROMAN);
				/**** more junk *****/
	line.ptlist = list;	/* used in plotting lines */
	line.n=2;		/* .. just a litte awkward */
}


/*
 * PLFRAME
 */
 
plframe()
{
	clear_control (CLEAR, NULL, NULL, VIEWPORT);
	clear_view_surface (name, OFF, NULL);
}


/*
 * PLFLUSH: output pending graphics
 */
 
plflush()
{
	/*	does nothing yet */
} 


/*
 * PLSTOP
 */
 
plstop()
{
	close_vws(name);
	close_cgi();
}


/*
 * PLLINE, PLMOVE, PLPOINT
 */
 
plline(x,y)
double x,y;
{
	trans (xp,yp,&list[0]);
	trans (x ,y ,&list[1]);
	polyline(&line);
	xp=x;					/* now reset pen pointer */
	yp=y;
}

plmove(x,y)
double x,y;
{
	xp=x;
	yp=y;
}
	
plpoint(x,y)
double x,y;
{
	xp=x;					/*  first reset pen pointer */
	yp=y;
	trans (xp,yp,&list[0]);
	trans (x ,y ,&list[1]);
	polyline(&line);
}


/*
 * PLTEXT: 
 */
 
pltext (msg, x, y, hgt, ang)
char *msg;
double x,y;
double hgt;
double ang;
{
	double c,s,sin(),cos();
	Ccoor  point,box;
		
	c=cos(ang/57.296);
	s=sin(ang/57.296);
	character_orientation(c,s,-s,c);
	character_path(RIGHT);

	trans(0.75*hgt,hgt,&box);
	character_height(box.y);
	
	trans(x,y,&point);
	text (&point,msg);		/* still unjustified etc. */
}

/*
 * PLJUST:
 */
 
static int	textjust =  -1;

pljust(jus)
int jus;
{	
	textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));
}


/*
 * PLLTYPE: does nothing now
 */
 
plltype (lwid, lpat)
int lwid;			/* line width */
int lpat;			/* line pattern */
{
}


/*
 * PLSWAP: does nothing
 */
 
plswap()
{
}

/*
 * PLXSCALE, PLYSCALE:  maybe this should be trans ?? 
 *			at the moment they are identity transformations
 */
 
double plxscale(x,y)
double x,y;
{
	return(x);
}

double plyscale(x,y)
double x,y;
{
	return(y);
}


/*
 * TRANS: transforms from user (world) coordinates to VDC coordinates
 */
 
trans(x,y,g)
double x,y;
Ccoor  *g;
{
	g->x = (int) (ax*x+bx);
	g->y = (int) (ay*y+by);
}


/*
 * WINDOW:  set user coordinates
 */
 
window (xmin, xmax, ymin, ymax)
double xmin, xmax;
double ymin, ymax;
{
	/*  does nothing, see plinit() */
}

	
/*
 *  MOUSE: does nothing but wait a little
 */
 
mouse(msec)
int msec;
{
/*	sleep(msec*1000);  */
	return(3);		/* do as if button 3 was hit */
}


plcircle()
{
}

plbox()
{
}

