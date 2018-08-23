/* 
 *  CCDPLOT: display an image in contour/gray-scale format 
 *
 *	30-jun-87  V1.0 created from CCD as separate module		PJT
 * 	 6-jul-87  V1.2 (does read image(5) files, order keywords optimized for cnt/gray  PJT
 *	 7-jul-87  V1.3 allow matrix size and plot size different
 *	 8-jul-87  V1.4 proper cell position definition (in the middle) PJT
 *	 1-jun-88  V2.0 rename programname look -> ccdplot	PJT
 *			also new filestruct, although code same
 *	31-jan-88  V2.1 axes names from header  PJT
 *			added headline= keyword and 
 *	 1-nov-90  V2.2 drange -> nemoindp, helpvec
 *      12-jan-92  V2.3 greyscale level info in plot header     PJT
 *	23-aug-92  V2.4 new version of contour() with lineto()  PJT
 *	 6-may-95  V2.5 minor fix
 *	22-feb-97  V2.5a support for SINGLEPREC			pjt
 *	 5-dec-99  V2.6 optional compilation with new pl_contour()	PJT
 *	 3-nov-00  V2.7 added blankval=		PJT
 * 	 7-may-01     a cleaned up some superfluous #define's	PJT
 *      19-feb-02  V2.8 added ltype=				pjt
 *      16-mar-05  V3.0 loop over all images if more found      PJT
 *                     a   added blankval to the new pl_matrix routine    PJT
 *      23-aug-18  V3.1  xscale,yscale,xlab,ylab                          PJT
 *	
 */

/* TO DO:
 *	sort input contours before plotting, now they have to be in order
 *	If used on a cube, can only plot first plane
 *	allow color input, use the color tables in $NEMODAT/lut
 *
 * BUGS:
 *	plotting subregions using xrange,yrange does not work
 */
  
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <yapp.h>
#include <axis.h>

string defv[] = {
	"in=???\n	input filename ",
	"contour=\n	contour levels ",
	"gray=\n	have gray scale ? (t/f) ",
	"mmin=\n	plot min for grayscale",
	"mmax=\n	plot max for grayscale [autoscale]",
	"power=1.0\n	power of transfer function ",
	"xrange=\n	range in x to plot -- BUG -- don't use ",
	"yrange=\n	range in y to plot -- BUG -- don't use ",
	"headline=\n	optional header above plot ",
	"tab=\n		table of contour segments",
	"format=%g\n	Format for above table",
	"cmode=0\n	Contour mode (0=orginal 1=pgplot)",
	"blankval=\n	if used, use this as blankval",
	"ltype=\n	width and type for contour",
	"xlab=\n	(override) Label along X-axis",
	"ylab=\n	(override) Label along Y-axis",
	"xscale=1\n     Scale all X values by this",
	"yscale=1\n     Scale all Y values by this",
	"VERSION=3.1\n	23-aug-2018 PJT",
	NULL,
};

string usage = "display an image using contours and/or gray-scales";

string cvsid = "$Id$";

#define UNDEF HUGE

string	infile;				/* input file */
stream  instr;
stream  tabstr=NULL;

real    mmin,mmax;			/* maximum in map */

bool    gray;
real    power;


imageptr iptr=NULL;			/* will be allocated dynamically */
int    nx,ny,nz,nsize;			/* size of matrix as allocated "" */
real   origin[2];			/* from header: in sky coordinates */
real xsize,ysize,size;		        /* frame dimension */
real xmin,ymin,dx,dy,cell;		/* frame dimensions */
real blankval;                          /* blank value */
real xscale, yscale;

#define MCNTVAL	100			/* maximum contour levels */
string cntstr;				/* string of contour levels */
real cntval[MCNTVAL];			/* contour levels */
int ncntval;				/* actual number of contours */

real xplot[2],   yplot[2];		/* ranges in plot variables */

string xlab=NULL, ylab=NULL;
char   xlabel[80], ylabel[80];              
char   plabel[80], clabel[80], glabel[80], tlabel[80];
string headline;			/* extra label */
char   format4[100];
int    cmode;
int    lwidth=0, ltype=0;

local real xtrans(real), ytrans(real);
void lineto(real, real, real, real);

extern int contour(real*, int, int, real*, int, real, real, real, real, proc);


nemo_main()
{
  int n = 0;
  setparams();                    /* set globals */

  instr = stropen (infile, "r");
  plinit ("***", 0.0, 20.0, 0.0, 20.0);       /* init yapp */
  while (read_image(instr,&iptr)) {   /* loop while more images found */
    dprintf(1,"Time= %g MinMax= %g %g\n",Time(iptr),MapMin(iptr),MapMax(iptr));
    nx=Nx(iptr);			
    ny=Ny(iptr);
    nz=Nz(iptr);
    if (nz > 1) error("Cannot handle 3D images [%d,%d,%d]",nx,ny,nz);
    xmin=Xmin(iptr);
    ymin=Ymin(iptr);
    dx=Dx(iptr) * xscale;
    dy=Dy(iptr) * yscale;
    xsize = nx * dx;
    ysize = ny * dy;
    if (n>0) {
      sleep(1);
      plframe();
    }
    plot_map();                                 /* plot the map */
    n++;
  }
  plstop();                                   /* end of yapp */
  strclose(instr);
}

setparams()
{
	string  tmpstr;
	int	tmpint, i;

	infile = getparam ("in");

	cntstr = getparam ("contour");	/* get contour levels */
	if (*cntstr==0) {
/*	    ncntval = 0;		/* this will cause a grayscale plot */
	    gray=TRUE;
	} else {
	    ncntval = nemoinpr(getparam("contour"),cntval,MCNTVAL);
	    tmpstr = getparam("gray");
	    if (*tmpstr==0)
	    	gray=FALSE;
	    else
	    	gray=getbparam("gray");
	}

	if (hasvalue("mmin"))
		mmin = getdparam("mmin");	
	else
		mmin = UNDEF;

	if (hasvalue("mmax"))
		mmax = getdparam("mmax");
	else
		mmax= UNDEF;
		
	power=getdparam("power");

	if (hasvalue("xrange"))
		setrange(xplot,getparam("xrange"));
	else
		xplot[0] = xplot[1] = UNDEF;
		
	if (hasvalue("yrange"))
		setrange(yplot,getparam("yrange"));
	else
		yplot[0] = yplot[1] = UNDEF;
			
        headline = getparam("headline");
        if (hasvalue("tab")) {
            tabstr = stropen(getparam("tab"),"w");
	    tmpstr = getparam("format");
            sprintf(format4,"%s %s %s %s\n",tmpstr,tmpstr,tmpstr,tmpstr);
	}
	cmode = getiparam("cmode");
	if (hasvalue("blankval")) {
	  blankval = getrparam("blankval");
	  contour_setdef(1,blankval);
	} else
	  blankval = -999.9; /* only relevant for greyscale now , use NaN ?? */
			   
	if (hasvalue("ltype")) {
            lwidth = getiparam("ltype");
            /* oops, skip ltype for now */
        }
	xscale = getrparam("xscale");
	yscale = getrparam("yscale");
	if (hasvalue("xlab")) xlab = getparam("xlab");
	if (hasvalue("ylab")) ylab = getparam("ylab");
            
}

plot_map ()
{
    real m_range, brightness, dcm;
    real m_min, m_max;
    int    i, ix, iy;
    int    cnt;		/* counter of pixels */

    nsize = Nx(iptr);           /* old method forced square .. */			
    cell = Dx(iptr) * xscale;   /* and forced so for gray scale due to LW mem-problems */
    size = nsize*cell;
    
    m_min = MapMin(iptr);	/* get min- and max from header */
    m_max = MapMax(iptr);
    dprintf (1,"Min and max in map from map-header are: %f %f\n",m_min,m_max);
 
    if (mmax==UNDEF)	/* reset default autoscales to user supplied if necessary */
	mmax=m_max;
    if (mmin==UNDEF)
    	mmin=m_min;
    m_range = mmax-mmin;
    if (m_range==0) {
    	mmax=mmin+1.0;
	if (gray) warning("%g; Plot-range was zero, mmax increased by 1",mmin);
    }
        
    dprintf (1,"User reset Min and max are: %f %f\n",mmin,mmax);

    
    sprintf (plabel,"File: %s",infile);			/* filename */
    sprintf (clabel,"Contours: %s",cntstr);             /* contour levels */
    sprintf (glabel,"Gray MinMax: %g %g",mmin,mmax);    /* grey scale minmax */
    sprintf (tlabel,"Time: %g",Time(iptr));             /* time of orig snapshot */

	/* set scales and labels along axes */
    if (xplot[0]==UNDEF || xplot[1]==UNDEF) {
    	xplot[0] = Xmin(iptr) - 0.5*Dx(iptr);
    	xplot[1] = xplot[0] + Nx(iptr)*Dx(iptr);
    	xplot[0] *= xscale;
    	xplot[1] *= xscale;
    }
    if (xlab != NULL)
        strncpy(xlabel,xlab,80);
    else if (Namex(iptr))
        strncpy(xlabel,Namex(iptr),80);
    else
        strcpy (xlabel,"");

    if (yplot[0]==UNDEF || yplot[1]==UNDEF) {
    	yplot[0] = Ymin(iptr) - 0.5*Dy(iptr);
    	yplot[1] = yplot[0] + Ny(iptr)*Dy(iptr);
    	yplot[0] *= yscale;
    	yplot[1] *= yscale;
	
    }
    if (ylab != NULL)
        strncpy(ylabel,ylab,80);
    else if (Namey(iptr))
        strncpy(ylabel,Namey(iptr),80);
    else
        strcpy (ylabel,"");
	
    dprintf (1,"Plotting area x=%f:%f y=%f:%f\n",
    	xplot[0], xplot[1], yplot[0], yplot[1]);

    if (gray) {		
				/* gray-scale */
       dcm = Dx(iptr) / (xplot[1]-xplot[0]) * 16.0 * xscale;
       pl_matrix (Frame(iptr), nx, ny, 
		  xtrans(Xmin(iptr)), ytrans(Ymin(iptr)), dcm , mmin, mmax, power, blankval);
           /*  color_bar (100.0,500.0,32);  */
     } 
     		/* OLD ROUTINE, has to call relocate/frame ---> plcontour */
     plltype(lwidth,ltype);
     if (cmode==0) 
        contour (Frame(iptr),nx,ny,cntval,ncntval,          /*  @todo  scaling not correct, see xplot/yplot[] */
		Xmin(iptr), 
		Ymin(iptr), 
		Xmin(iptr)+(Nx(iptr)-1)*Dx(iptr)*xscale,
		Ymin(iptr)+(Ny(iptr)-1)*Dy(iptr)*yscale, lineto);
    
     else if (cmode==1)
         pl_contour (Frame(iptr),nx,ny,ncntval,cntval);
     plltype(1,1);

	/* draw axes and their labels */
    xaxis ( 2.0, 2.0, 16.0, xplot, -7, xtrans,  xlabel);
    xaxis ( 2.0,18.0, 16.0, xplot, -7, xtrans,  NULL);
    yaxis ( 2.0, 2.0, 16.0, yplot, -7, ytrans, ylabel);
    yaxis (18.0, 2.0, 16.0, yplot, -7, ytrans, NULL);

    pltext (plabel,2.0,18.4, 0.32, 0.0);	/* plot header with file name */
    pltext (clabel,2.0,19.0, 0.32, 0.0);        /* plot header with contour levels */
    pltext (glabel,2.0,19.6, 0.32, 0.0);        /* plot header with greyscale info */
    pltext (tlabel,10.0,19.6,0.32, 0.0);        /* time info */

    pljust(1);
    pltext (headline,10.0,18.4, 0.26, 0.0);     /* plot extra user suplied header */
    pljust(-1);
}

#if 0
color_bar(xmin, xmax, ncolors)	/* should go to yapp */
real xmin, xmax;
int ncolors;
{
	int i;
	real graylevel, x;
	char line[200];

	dx = (xmax-xmin)/ncolors;
	x = xmin;
	graylevel = 0.0;
	for (i=0; i<ncolors; i++) {

		vec_paint (x,50.0,dx,graylevel);
		x += dx;
	}
	sprintf (line,"%.1f %.1f mt %.1f %.1f lt s",xmin,45.0,xmax,45.0);
	vec_write (line);
}
#endif

/*********** SOME LOCAL UTILITIES , PLOTTING ANDSO ***********/
/* some functions for yapp-plotting */

local real xtrans(real x)
{
    return 2.0 + 16.0*(x-xplot[0])/(xplot[1]-xplot[0]) ;
}


local real ytrans(real y)	/* upper panel: mass surface density */
{
    return 2.0 + 16.0*(y-yplot[0])/(yplot[1]-yplot[0]) ;
}

void lineto(real x1, real y1, real x2, real y2)
{
    if (tabstr) fprintf(tabstr,format4,x1,y1,x2,y2);

    plmove (xtrans(x1),ytrans(y1));
    plline (xtrans(x2),ytrans(y2));
}
		
/*******  stolen from snapplot.c **********/

setrange(real *rval, string rexp)
{
    char  *cptr;
	
    cptr = strchr(rexp, ':');
    if (cptr != NULL) {
        rval[0] = atof(rexp);
        rval[1] = atof(cptr+1);
    } else {
        rval[0] = 0.0;
        rval[1] = atof(rexp);
    } 
}
					
					
