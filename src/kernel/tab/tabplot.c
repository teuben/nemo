/*
 * TABPLOT: a simple general table plotter
 *	    for ascii data in tabular format which has now grown quite complex
 *
 *	25-nov-88  V1.0 :  created by P.J.Teuben
 *	 5-jul-89      a:  nemoinp to nemoinp(d/i)
 *	13-nov-90  V1.1 :  helpvec, get_atable()
 *      26-mar-92  V1.2 :  added nmax=
 *       8-apr-92  V1.3 :  skip plotting of no points found, usage;
 *			   no more 'nskip'; error= -> errors=
 *	24-apr-92  V1.4 :  repaired xbin=  - PJT
 *	 5-nov-93      a:  get_atable: warn about partial pipe
 *	24-jan-95  V1.5 :  support for reading multiple columns (ycol=)
 *                     a:  added xcoord= and ycoord=
 *      13-feb-96  V1.6 :  added cursor= keyword to get positions back
 *                         into this file (only works with PGPLOT)
 *			   and added layout= keyword
 *      16-aug-96      b : increased max # Y-columns from 16 to 64
 *	14-sep-96      c : bug in reporting ycol
 *	16-feb-97      d : support for SINGLEPREC
 *	20-feb-97  V1.7  : removed never used skip= keyword, MAXYCOL now 256
 *       1-apr-97      a : options pl_cursor compilation
 *	10-oct-97      b : minor debug output change
 *	 7-may-98  V1.8  : added median option for binning
 *	13-jul-98      a : printf->dprintf for fixed binning
 *      25-jul-98      b : allow line= with histogram-type (ltype<0)
 *      31-mar-99  V2.1  : allow autoscaling with a range set in the
 *                         other coordinate
 *	28-jul-99  V2.2  : add color=
 *      21-jul-00  V2.3  : allow autoscaling on half (min OR max) an axis
 *      23-sep-01      a : ->nemo_file_lines
 *       6-oct-01  V2.4  : fixed Y autoscale if all Y's have the same value
 *       2-aug-02  V2.5  : allow multiple X columns, forces pairwise xcol,ycol plotting
 *                     b : fix old semi-autoscaling bug while fixing it for multicol
 *      31-dec-03  V2.6  : option to do layout first
 *      28-apr-04  V2.7  : added xbox, ybox= options as in snapplot
 *      17-sep-05  V2.8  : added pl_readlines()
 *       2-dec-05  V2.9  : implemented many-to-one column plotting mode
 *                 V3.0  : xscale,yscale,dxcol,dycol=
 *      20-dec-05  V3.0a : fix bug in xbin=min:max:step mode
 *                     b : fix bug in bins with no dispersive data, switch to moment.h
 *       9-oct-06      d : Implemented the dxcol= and dycol=
 *       8-apr-11      e : fixed dycol reference bug
 *      22-aug-2018 V3.1 : print commented line if bin= causes an empty bin
 *       8-jan-2020 V4.0 : template python option
 */

/* TODO:
 *     - automatic scaling with dxcol and/or dycol 
 *     - also allow nxcol>1 and nycol=1 ??
 */

/**************** INCLUDE FILES ********************************/ 
 
#include <stdinc.h>	
#include <getparam.h>
#include <yapp.h>
#include <axis.h>
#include <layout.h>
#include <table.h>
#include <pyplot.h>
#include <moment.h>
                    /* undefined values trick !!!  MACHINE DEP  !!! */
#ifdef SINGLEPREC
#define NaN 0x7FFF
#else
#define NaN 0x7FFFFFFF
#endif

#define MAXCOL  256
#define MAXCOORD 16

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n            Input file name (table)",
    "xcol=1\n		 x coordinate column(s)",
    "ycol=2\n		 y coordinate column(s)",
    "xmin=\n		 X-min value if not autoscaled",
    "xmax=\n		 X-max value if not autoscaled",
    "ymin=\n		 Y-min value if not autoscaled",
    "ymax=\n		 Y-max value if not autoscaled",
    "xbin=\n             binning data in X?",
    "point=1,0.1\n       point type and size pairs for each Y column (yapp)",
    "line=0,0\n          line width and style pairs for each Y column (yapp)",
    "color=\n		 colors for the symbols/lines for each Y column",
    "errors=\n           plot wich of X-Y error bars ?",
    "xlab=\n		 Label along X-axis",
    "ylab=\n		 Label along Y-axis",
    "xcoord=\n		 Draw additional vertical coordinate lines along these X values",
    "ycoord=\n		 Draw additional horizontal coordinate lines along these Y values",
    "xbox=2.0:18.0\n	 extent of frame in x direction",
    "ybox=2.0:18.0\n     extent of frame in y direction",
    "nxticks=7\n         Number of ticks on X axis",
    "nyticks=7\n         Number of ticks on Y axis",
    "xscale=1\n          Scale all X values by this",
    "yscale=1\n          Scale all Y values by this",
    "dxcol=\n            Columns representing error bars in X ",
    "dycol=\n            Columns representing error bars in Y ",
    "nmax=100000\n       Hardcoded allocation space, if needed for piped data",
    "median=f\n          Use median in computing binning averages?",
    "headline=\n	 headline in graph on top right",
    "tab=f\n             Table output also if binning is used?",
    "fullscale=f\n       Use full autoscale in one axis if other autoscaled?",
    "cursor=\n           Optional output file to retrieve cursor coordinates",
    "layout=\n           Optional input layout file",
    "first=f\n           Layout first or last?",
    "readline=f\n        Interactively reading commands",
    "pyplot=\n           Template python plotting script",
    "VERSION=4.0\n	 8-jan-2020 PJT",
    NULL
};

string usage = "general table plotter";

string cvsid="$Id$";


/**************** GLOBAL VARIABLES ************************/

local string input;				/* filename */
local stream instr;				/* input file */

local int xcol[MAXCOL], ycol[MAXCOL], nxcol, nycol;	/* column numbers */
local int dxcol[MAXCOL],dycol[MAXCOL],ndxcol,ndycol;	/* column numbers */
local real xrange[2], yrange[2];		/* range of values */
local int npcol;                                /* MAX(nxcol,nycol) */

local real  *x[MAXCOL], *y[MAXCOL]; 			/* data from file */
local real  *dx[MAXCOL], *dy[MAXCOL]; 			/* data from file */
local real *xp, *yp, *xps, *yps;              /* rebinned data (if needed) */
local int    npt;				/* actual number of data points */
local int    np;                              /* actual number of rebinned data */

local real xmin, xmax;			/* actual min and max as computed */
local real ymin, ymax;
local bool   Qautox0, Qautoy0, Qfull;		/* autoscale ? */
local bool   Qautox1, Qautoy1;
local bool   Qsigx, Qsigy;                    /* show error bars on data points ? */
local bool   Qtab;				/* table output ? */
local bool   Qmedian;
local real *xbin;                             /* boundaries for binned data */
local int    nbin;                            /* number of boundaries " */

local int    nmax;				/* lines to allocate */

local string headline;			/* text string above plot */
local string xlab, ylab;			/* text string along axes */
local real xplot[2],yplot[2];		        /* borders of plot */
local int    yapp_line[2*MAXCOL];            /* thickness and pattern of line */
local int    yapp_color[MAXCOL];	      /* color per column */
local real yapp_point[2*MAXCOL];             /* what sort of point to plot */
local string errorbars;                       /* x or y or both (xy) */
local real xcoord[MAXCOORD], ycoord[MAXCOORD];/* coordinate lines */
local int nxcoord, nycoord;
local int nxticks, nyticks;                   /*& number of tickmarks */
local int ncolors;
local real xbox[3], ybox[3];
local real xscale, yscale;

local plcommand *layout;
local bool layout_first;
local bool Qreadlines;

void setparams(), plot_data();

extern real median(int, real *);
extern int nemo_file_lines(string, int);

local real  xtrans(real),  ytrans(real);    /* WC -> cmXY */
local real ixtrans(real), iytrans(real);    /* cmXY -> WC */
local void setrange(real *rval, string rexp);

/****************************** START OF PROGRAM **********************/

void nemo_main(void)
{
    setparams();

    if (hasvalue("pyplot")) {
      stream pstr = pyplot_init(getparam("pyplot"));
      pyplot_plot(pstr, input, xcol, ycol, dycol, xrange, yrange);
      pyplot_close(pstr);
    }
      

    instr = stropen (input,"r");
    read_data();
    plot_data();
}

void setparams()
{
    char *smin,*smax;
    int  i, j;
   
    input = getparam("in");             /* input table file */
    nxcol = nemoinpi(getparam("xcol"),xcol,MAXCOL);
    nycol = nemoinpi(getparam("ycol"),ycol,MAXCOL);
    if (nxcol < 1) error("Error parsing xcol=%s",getparam("xcol"));
    if (nycol < 1) error("Error parsing ycol=%s",getparam("ycol"));
    if (nxcol > 1 && nycol > 1 && nycol != nxcol) 
      error("nxcol=%d nycol=%d cannot be paired properly",nxcol,nycol);
    npcol = MAX(nxcol,nycol);
    ndxcol = nemoinpi(getparam("dxcol"),dxcol,MAXCOL);
    ndycol = nemoinpi(getparam("dycol"),dycol,MAXCOL);
    // TODO: some ndxycol/ndycol checks
    

    nmax = nemo_file_lines(input,getiparam("nmax"));
    dprintf(1,"Allocated %d lines for table\n",nmax);
    for (j=0; j<nxcol; j++)
        x[j] = (real *) allocate(sizeof(real) * (nmax+1));   /* X data */
    for (j=0; j<nycol; j++)
        y[j] = (real *) allocate(sizeof(real) * (nmax+1));   /* Y data */
    if (ndxcol)
      for (j=0; j<ndxcol; j++)
        dx[j] = (real *) allocate(sizeof(real) * (nmax+1));  /* DX data */
    else
      for (j=0; j<ndxcol; j++)
        dx[j] = 0;
    if (ndycol)
      for (j=0; j<ndycol; j++)
        dy[j] = (real *) allocate(sizeof(real) * (nmax+1));  /* DY data */
    else
      for (j=0; j<ndycol; j++)
        dy[j] = 0;


    Qsigx = Qsigy = FALSE;
    if (hasvalue("xbin")) {
      if (nxcol > 1) error("Cannot bin in X, more than 1 column is used");
      smin = getparam("xbin");            /* binning of data */
      Qsigx = Qsigy = TRUE;           /* binning of data */
      xbin = (real *) allocate(nmax*sizeof(real));
      nbin = nemoinpr(smin,xbin,nmax);   /* get binning boundaries */
      if (nbin==1) {              /*  fixed amount of datapoints to bin */
	(void) nemoinpi(smin,&nbin,1);
	dprintf(0,"Binning with fixed number (%d) of points\n",nbin);
	np = nmax / nbin + 1;
	nbin = -nbin;       /* make it <0 to trigger rebin_data */
      } else if (nbin>1) {
	np = nbin - 1;
	if (np==1)
	  warning("plotting one averaged datapoint");
      } else
	error("Error %d in parsing xbin=%s",nbin,smin);

      xp =  (real *) allocate(np*sizeof(real));    /* X data */
      yp =  (real *) allocate(np*sizeof(real));    /* Y data */
      xps = (real *) allocate(np*sizeof(real));   /* sig-X data */
      yps = (real *) allocate(np*sizeof(real));   /* sig-Y data */

      for (i=0; i<np; i++)
	xp[i] = xps[i] = yp[i] = yps[i] = NaN;      /* init empty */
    } else
      np = 0;
    
    if (hasvalue("xmin") && hasvalue("xmax")) {
	xrange[0] = getdparam("xmin");
	xrange[1] = getdparam("xmax");
	Qautox0=Qautox1=FALSE;
	dprintf (2,"fixed plotrange in X: %f : %f\n",xrange[0],xrange[1]);
    } else if (hasvalue("xmin")) {
	xrange[0] = getdparam("xmin");
	Qautox0=FALSE;
	Qautox1=TRUE;
	dprintf (2,"fixed minimum in X: %f\n",xrange[0]);
    } else if (hasvalue("xmax")) {
	xrange[1] = getdparam("xmax");
	Qautox0=TRUE;
	Qautox1=FALSE;
	dprintf (2,"fixed maximum in X: %f\n",xrange[1]);
    } else {
    	xrange[0] = -HUGE;
    	xrange[1] = HUGE;
	Qautox0=Qautox1=TRUE;
	dprintf (2,"auto plotscaling\n");
    }

    if (hasvalue("ymin") && hasvalue("ymax")) {
	yrange[0] = getdparam("ymin");
	yrange[1] = getdparam("ymax");
	Qautoy0=Qautoy1=FALSE;
	dprintf (2,"fixed plotrange in Y: %f : %f\n",yrange[0],yrange[1]);
    } else if (hasvalue("ymin")) {
	yrange[0] = getdparam("ymin");
	Qautoy0=FALSE;
	Qautoy1=TRUE;
	dprintf (2,"fixed minimum in Y: %f\n",yrange[0]);
    } else if (hasvalue("ymax")) {
	yrange[1] = getdparam("ymax");
	Qautoy0=TRUE;
	Qautoy1=FALSE;
	dprintf (2,"fixed maximum in Y: %f\n",yrange[1]);
    } else {
    	yrange[0] = -HUGE;
    	yrange[1] = HUGE;
	Qautoy0=Qautoy1=TRUE;
	dprintf (2,"auto plotscaling\n");
    }

    setrange(xbox, getparam("xbox"));
    setrange(ybox, getparam("ybox"));

    xscale = getrparam("xscale");
    yscale = getrparam("yscale");


    Qfull = getbparam("fullscale");

    parse_pairsi("line", yapp_line, npcol);  
    parse_pairsr("point",yapp_point,npcol);
    ncolors = nemoinpi(getparam("color"),yapp_color,npcol);
    if (ncolors > 0) {
    	for (i=ncolors; i<npcol; i++)
    	   yapp_color[i] = yapp_color[i-1];
    }
       
    smin = getparam("errors");               /* force error bars in graph ? */
    if (*smin) {
        Qsigx = (strchr(smin,'x') != NULL);	/* error bar in X needed ?*/
        Qsigy = (strchr(smin,'y') != NULL);	/* error bar in Y needed ?*/
	if (nxcol>0 || nycol>0) error("Cannot plot errors from binning and normal column errors");
    } 
    Qtab = getbparam("tab");
    Qmedian = getbparam("median");
    xlab=getparam("xlab");
    ylab=getparam("ylab");
    headline = getparam("headline");
    nxcoord = nemoinpr(getparam("xcoord"),xcoord,MAXCOORD);
    nycoord = nemoinpr(getparam("ycoord"),ycoord,MAXCOORD);
    nxticks = getiparam("nxticks");
    nyticks = getiparam("nyticks");
    if (hasvalue("layout"))
        layout = pl_fread(getparam("layout"));
    else
        layout = NULL;
    layout_first = getbparam("first");
    Qreadlines = getbparam("readline");
}

#define MVAL 		 64
#define MLINELEN	512

read_data()
{
    real *coldat[1+MAXCOL];
    int i, j, k, colnr[1+MAXCOL];
		
    dprintf (2,"Reading datafile, xcol,ycol=%d..,%d,...\n",xcol[0],ycol[0]);
    for (j=0, k=0; j<nxcol; j++, k++) {
        colnr[k]  = xcol[j];
        coldat[k] = x[j];
    }
    for (j=0; j<ndxcol; j++, k++) {
        colnr[k]  = dxcol[j];
        coldat[k] = dx[j];
    }
    for (j=0; j<nycol; j++, k++) {
        colnr[k]  = ycol[j];
        coldat[k] = y[j];
    }
    for (j=0; j<ndycol; j++, k++) {
        colnr[k]  = dycol[j];
        coldat[k] = dy[j];
    }
    

    /* could also find out if any columns duplicated, and
       replace them with pointers */

    npt = get_atable(instr,nxcol+nycol+ndxcol+ndycol,colnr,coldat,nmax);    /* get data */
    dprintf(1,"get_atable: %d\n",npt);
    if (npt < 0) {
    	npt = -npt;
    	warning("Could only read first set of %d data",npt);
    }
    if (xscale != 1.0) {
      for (j=0; j<nxcol; j++)
	for (i=0; i<npt; i++)
	  x[j][i] *= xscale;
      for (j=0; j<ndxcol; j++)
	for (i=0; i<npt; i++)
	  dx[j][i] *= xscale;
    }
    if (yscale != 1.0) {
      for (j=0; j<nycol; j++)
	for (i=0; i<npt; i++)
	  y[j][i] *= yscale;
      for (j=0; j<ndycol; j++)
	for (i=0; i<npt; i++)
	  dy[j][i] *= yscale;
    }

    xmin = ymin =  HUGE;
    xmax = ymax = -HUGE;
    for (i=0; i<npt; i++) {     /* find global min and max in X and all Y */
        for (j=0; j<nxcol; j++) {
	  if (ndxcol > 0) {
	    if (ndxcol != nxcol) error("Cannot use X errors if number of columns do not match");
	    xmax=MAX(x[j][i]+dx[j][i],xmax);
	    xmin=MIN(x[j][i]-dx[j][i],xmin);
	  } else {
	    xmax=MAX(x[j][i],xmax);
	    xmin=MIN(x[j][i],xmin);
	  }
        }
        for (j=0; j<nycol; j++) {
	  if (ndycol > 0) {
	    if (ndycol != nycol) error("Cannot use Y errors if number of columns do not match");
	    ymax=MAX(y[j][i]+dy[j][i],ymax);
	    ymin=MIN(y[j][i]-dy[j][i],ymin);
	  } else {
	    ymax=MAX(y[j][i],ymax);
	    ymin=MIN(y[j][i],ymin);
	  }
        }
    }
}


void plot_data()
{
    int    i, j, k;
    real xcur, ycur, edge;
    char c;
    stream cstr;
    bool first, Qin;

    if (npt<1) {
        warning("Nothing to plot, npt=%d",npt);
        return;
    }
    dprintf (0,"read %d points\n",npt);
    dprintf (0,"min and max value in xcolumns %s: [%f : %f]\n",
	     getparam("xcol"),xmin,xmax);
    dprintf (0,"min and max value in ycolumns %s: [%f : %f]\n",
	     getparam("ycol"),ymin,ymax);
    if (Qautox0 && Qautox1) {
	    edge = (xmax-xmin) * 0.05;        /* add 5% to the edges */
	    xmin -= edge;
	    xmax += edge;
    } else if (Qautox0) {
	    xmax = xrange[1];
            edge = (xmax-xmin) * 0.05;        /* add 5% to the edge */
            xmin -= edge;
    } else if (Qautox1) {
	    xmin = xrange[0];
            edge = (xmax-xmin) * 0.05;        /* add 5% to the edge */
            xmax += edge;
    } else {
	    xmin = xrange[0];
	    xmax = xrange[1];
    }
    if (Qautox0 || Qautox1) 
        dprintf(0,"X:min and max value reset to : [%f : %f]\n",xmin,xmax);

    if (Qautoy0 && Qautoy1) {
      edge = (ymax-ymin) * 0.05;        /* add 5% to the edges */
      if (edge == 0.0) edge = 0.05*ABS(ymax);
      if (edge == 0.0) edge = 1;
      ymin -= edge;
      ymax += edge;
    } else if (Qautoy0) {
      ymax = yrange[1];
      edge = (ymax-ymin) * 0.05;        /* add 5% to the edge */
      if (edge == 0.0) edge = 0.05*ABS(ymin);
      if (edge == 0.0) edge = 1;
      ymin -= edge;
    } else if (Qautoy1) {
      ymin = yrange[0];
      edge = (ymax-ymin) * 0.05;        /* add 5% to the edge */
      if (edge == 0.0) edge = 0.05*ABS(ymax);
      if (edge == 0.0) edge = 1;
      ymax += edge;
    } else {
      ymin = yrange[0];
      ymax = yrange[1];
    }
    if (Qautoy0 || Qautoy1)
      dprintf(0,"Y:min and max value reset to : [%f : %f]\n",ymin,ymax);

    /*
     *  if scale in X fully fixed, and some autoscaling in Y is done
     *  we can recompute the Ymin and/or Ymax based on fixed Xmin/Xmax
     *  for those points which are within Xmin/Xmax
     */
    if (!Qfull && (Qautoy0 || Qautoy1) && !Qautox0 && !Qautox1) {
      first = TRUE;
      for (i=0; i<npt; i++) {     /* go through data, find X min and max again */
	for (j=0; j<nxcol; j++) {
	  if (xmin < xmax && (x[j][i]<xmin || x[j][i]>xmax)) break;
	  if (xmax < xmin && (x[j][i]<xmax || x[j][i]>xmin)) break;
	}
	if (j<nxcol) continue;
	
	if (first) {
	  if (Qautoy0) ymin = y[0][i];
	  if (Qautoy1) ymax = y[0][i];
	  first = FALSE;
	}
	for (j=0; j<nycol; j++) {
	  if (Qautoy0) ymin=MIN(y[j][i],ymin);
	  if (Qautoy1) ymax=MAX(y[j][i],ymax);
	}
      }
      dprintf(0,"Y:min and max value re-reset to : [%f : %f]\n",ymin,ymax);
    }

    /*
     * now the same, but reversed axes...
     */

    if (!Qfull && (Qautox0 || Qautox1) && !Qautoy0 && !Qautoy1) {
        first = TRUE;
        for (i=0; i<npt; i++) {     /* go through data, find X min and max again */
	  for (j=0; j<nycol; j++) {
	    if (ymin < ymax && (y[j][i]<ymin || y[j][i]>ymax)) break;
	    if (ymax < ymin && (y[j][i]<ymax || y[j][i]>ymin)) break;
	  }
	  if (j<nycol) continue;
	  
	  if (first) {
	    if (Qautox0) xmin = x[0][i];
	    if (Qautox1) xmax = x[0][i];
	    first = FALSE;
	  }
	  for (j=0; j<nxcol; j++) {
	    if (Qautox0) xmin=MIN(x[j][i],xmin);
	    if (Qautox1) xmax=MAX(x[j][i],xmax);
	  }
        }
        dprintf(0,"X:min and max value re-reset to : [%f : %f]\n",xmin,xmax);
	/* should 5% be added to the edges ? */
    }
    
    plinit("***",0.0,20.0,0.0,20.0);                 /*	PLOTTING */	
    if (layout && layout_first) pl_exec(layout);

    xplot[0] = xmin;        /* set scales for xtrans() */
    xplot[1] = xmax;
    yplot[0] = ymin;        /* set scales for ytrans() */
    yplot[1] = ymax;
    xaxis (xbox[0],ybox[0],xbox[2], xplot, -nxticks, xtrans, xlab);
    xaxis (xbox[0],ybox[1],xbox[2], xplot, -nxticks, xtrans, NULL);
    yaxis (xbox[0],ybox[0],ybox[2], yplot, -nyticks, ytrans, ylab);
    yaxis (xbox[1],ybox[0],ybox[2], yplot, -nyticks, ytrans, NULL);
    for (i=0; i<nxcoord; i++) {
        plmove(xtrans(xcoord[i]),ytrans(yplot[0]));
        plline(xtrans(xcoord[i]),ytrans(yplot[1]));
    }
    for (i=0; i<nycoord; i++) {
        plmove(xtrans(xplot[0]),ytrans(ycoord[i]));
        plline(xtrans(xplot[1]),ytrans(ycoord[i]));
    }

    pljust(-1);     /* set to left just */
    pltext(input,xbox[0],ybox[1]+0.2,0.32,0.0);             /* filename */
    pljust(1);
    pltext(headline,xbox[1],ybox[1]+0.2,0.24,0.0);         /* headline */
    pljust(-1);     /* return to left just */

    for (k=0, i=0, j=0; k<npcol; k++) {
      if (np && nxcol==1) rebin_data (npt,x[i],y[j], nbin, xbin, np, xp, yp,  xps, yps);

      plot_points( npt, x[i], y[j], dx[i], dy[j], xps, yps,
		   yapp_point[2*k], yapp_point[2*k+1],
		   yapp_line[2*k], yapp_line[2*k+1],
		   ncolors > 0 ? yapp_color[k] : -1,
		   0);
      if (nxcol>1) i++;
      if (nycol>1) j++;
    }

    if (hasvalue("cursor")) {
        dprintf(0,"Interactive mode: enter coordinates with the cursor\n");
        cstr = NULL;
        while (pl_cursor(&xcur, &ycur, &c)) {
            dprintf(0,"%c: %g %g cm (%g %g)\n",
                c,xcur,ycur,ixtrans(xcur),iytrans(ycur));
            if (c=='X') break;
            if (c=='A') {
                if (cstr == NULL) cstr = stropen(getparam("cursor"),"w");
                fprintf(cstr,"%g %g\n",ixtrans(xcur),iytrans(ycur));
            }
        }
        if (cstr) strclose(cstr);
    }
    if (layout && !layout_first) pl_exec(layout);
    if (Qreadlines) pl_readlines();

    plstop();
}

double my_sqrt(double x)
{
  if (x <= 0.0) return 0.0;
  return sqrt(x);
}

/*
 *  REBIN_DATA:
 *
 *      nbin > 0:   array xbin[0..nbin-1] countaints bin boundaries for x
 *      nbin < 0:   bin x data in groups of (-nbin) points into (xp,yp)
 *  
 *      The arrays (xp,yp) contain resulting data, (xps,yps) the dispersion,
 *      which can be used for plotting error bars
 *
 *      It is assumed that the X-data are sorted
 */
 
rebin_data (n,x,y, nbin,xbin, np, xp, yp,  xps, yps)
int n, nbin, np;
real *x, *y, *xbin, *xp, *yp, *xps, *yps;
{
    int    i, j, ip, iold=0, zbin=0;
    Moment mx, my;

    dprintf(0,"Rebinning...n=%d nbin=%d np=%d\n",n,nbin,np);
    
    for (i=0; i<np; i++)
        xp[i] = xps[i] = yp[i] = yps[i] = NaN;      /* init (again) */
   
    for (i=1, ip=0; i<n; i++)		/* loop over all points */
      if (x[i] < x[i-1]) {
	dprintf(1,"Unsorted %g < %g at %d\n",x[i],x[i-1],i);
	ip++;              /* count unsorted data in 'ip' */
      }
    if (ip>0)
        warning("%d out of %d datapoints are not sorted\n", ip, n);

    i = 0;                     /* point to original data (x,y) to be rebinned */
    iold = 0;
    if (nbin>0) {
      while(i<n && x[i] < xbin[0]) i++;
      if (i>0) warning("%d points left of first bin",i);
    }
    for (ip=0; ip<nbin-1; ip++) {		/* for each bin, accumulate */
	ini_moment(&mx, 2, 0);
	ini_moment(&my, 2, 0);
	dprintf(1,"ip=%d    %g : %g\n",ip,xbin[ip],xbin[ip+1]);
        while (i<n && xbin[ip] <= x[i] && x[i] < xbin[ip+1]) {    /* in this bin ? */
	  accum_moment(&mx,x[i],1.0);
	  accum_moment(&my,y[i],1.0);
	  i++;
        }
        if (n_moment(&mx) > 0) {
	  xp[ip]  = mean_moment(&mx);
	  xps[ip] = sigma_moment(&mx);
	  yp[ip]  = mean_moment(&my);
	  yps[ip] = sigma_moment(&my);
	  if (Qmedian) {
	    xp[ip] = median(i-iold, &x[iold]);
	    yp[ip] = median(i-iold, &y[iold]);
	  }
        } else
            zbin++;     /* count bins with no data */
	if(Qtab)
	    if (n_moment(&mx))    /* print non-zero bins */
	      printf("%g %g %g %g %d\n",xp[ip],yp[ip],xps[ip],yps[ip],i-iold);
            else
	      printf("# 0\n");    /* comment line for empty bin */
        iold = i;
    } /* for(ip) */
    if(zbin)warning("There were %d bins with no data",zbin);
    if(i<n) warning("%d points right of last bin",i);
    if (nbin>0)
        return 0;      /* done with variable bins */

    nbin = -nbin;       /* make it positive for fixed binning */
    for (i=0, ip=0; i<n;ip++) {       /* loop over all points */
      ini_moment(&mx, 2, 0);
      ini_moment(&my, 2, 0);
      for(j=0; j<nbin && i<n; j++, i++) {   /* accum the new bin */
	accum_moment(&mx,x[i],1.0);
	accum_moment(&my,y[i],1.0);
      }
      xp[ip] = mean_moment(&mx);
      yp[ip] = mean_moment(&my);
      xps[ip] = sigma_moment(&mx);
      yps[ip] = sigma_moment(&my);
      if(Qtab)
	printf("%g %g %g %g %d\n",xp[ip],yp[ip],xps[ip],yps[ip],i-iold);
      iold = i;
    }
}

/* PLOT_POINTS
 *
 *      Plot datapoints, optionally connecting them with lines,
 *      adding errorbars in X as well as Y
 *
 *      X,Y Points with value NaN are skipped, as well as their errorbars
 */
 
plot_points (np, xp, yp, dx, dy, xps, yps, pstyle, psize, lwidth, lstyle, color, errorbars)
int np, lwidth, lstyle, color, errorbars;
real xp[], yp[], dx[], dy[], xps[], yps[], pstyle, psize;
{
    int i, ipstyle;
    real p1, p2;

    if (color >= 0) {
        plcolor(color);
    }
    
    ipstyle = pstyle+0.1;
    if (ipstyle != 0)
        for (i=0; i<np; i++)
            if (xp[i] != NaN && yp[i] != NaN)
                switch (ipstyle) {
                case 1:
                    plpoint (xtrans(xp[i]), ytrans(yp[i]));
                    break;
                case 2:
                    plcircle (xtrans(xp[i]), ytrans(yp[i]),psize);
                    break;
                case 3:
                    plcross (xtrans(xp[i]), ytrans(yp[i]),psize);
                    break;
                case 4:
                    plbox (xtrans(xp[i]), ytrans(yp[i]),psize);
                    break;
                default:
                    error ("Invalid pstyle = %d\n",pstyle);
                }


    if (lwidth > 0) {
        plltype(lwidth,ABS(lstyle));
        i=0;
        while (xp[i] == NaN || yp[i] == NaN)
            i++;                             /* skip first undefined */
	if (lstyle > 0) {       /* connect the dots */
            plmove (xtrans(xp[i]), ytrans(yp[i]));      /* move to point */
            while (++i < np)
                if (xp[i] != NaN)
                    plline (xtrans(xp[i]), ytrans(yp[i]));  /* draw line */

        } else {                /* histogram approach */
            plmove (xtrans(xp[i]), ytrans(yp[i]));      /* move to 1st point */
            while (++i < np) {
                plline(xtrans(0.5*(xp[i-1]+xp[i])),ytrans(yp[i-1]));
                plline(xtrans(0.5*(xp[i-1]+xp[i])),ytrans(yp[i]));
                plline(xtrans(xp[i]),ytrans(yp[i]));
            }            
        }
        plltype(1,1);
    }
    
    if (errorbars && 0x0001) {
        dprintf(0,"Trying X binning errors\n");
        for (i=0; i<np; i++) {
            if (xps[i] == NaN || xp[i] == NaN)
                continue;
            p1 = xp[i] - xps[i];
            p2 = xp[i] + xps[i];
            plmove (xtrans(p1), ytrans(yp[i]));
            plline (xtrans(p2), ytrans(yp[i]));
        }
    }
    if (errorbars && 0x0002) {
        dprintf(0,"Trying Y binning errors\n");
        for (i=0; i<np; i++) {
            if (yps[i] == NaN || yp[i] == NaN)
                continue;
            p1 = yp[i] - yps[i];
            p2 = yp[i] + yps[i];
            plmove (xtrans(xp[i]), ytrans(p1));
            plline (xtrans(xp[i]), ytrans(p2));
        }
    }

    if (dx) {
        dprintf(0,"Trying X column errors\n");
        for (i=0; i<np; i++) {
            p1 = xp[i] - dx[i];
            p2 = xp[i] + dx[i];
            plmove (xtrans(p1), ytrans(yp[i]));
            plline (xtrans(p2), ytrans(yp[i]));
        }
    }
    if (dy) {
        dprintf(0,"Trying Y column errors\n");
        for (i=0; i<np; i++) {
            p1 = yp[i] - dy[i];
            p2 = yp[i] + dy[i];
            plmove (xtrans(xp[i]), ytrans(p1));
            plline (xtrans(xp[i]), ytrans(p2));
        }
    }
}

local real xtrans(real x)
{
  return xbox[0] + xbox[2] * (x - xplot[0]) / (xplot[1]-xplot[0]);
}

local real ytrans(real y)
{
  return ybox[0] + ybox[2] * (y - yplot[0]) /(yplot[1]-yplot[0]); 
}

local real ixtrans(real x)
{
  return xplot[0] + (x-xbox[0])*(xplot[1]-xplot[0])/xbox[2];
}
 
local real iytrans(real y)
{
  return yplot[0] + (y-xbox[0])*(yplot[1]-yplot[0])/ybox[2];
        
}       

local void setrange(real *rval, string rexp)
{
    char *cptr, *tmpstr;
    double dpar;

    cptr = strchr(rexp, ':');
    if (cptr != NULL) {
        tmpstr = allocate(cptr-rexp+1);
        strncpy(tmpstr,rexp,cptr-rexp);
        if (nemoinpd(tmpstr,&dpar,1) != 1)
            error("setrange: parsing error %s",tmpstr);
        free(tmpstr);
        rval[0] = dpar;

        if (nemoinpd(cptr+1,&dpar,1) != 1)
            error("setrange: parsing error %s",cptr+1);
	rval[1] = dpar;
    } else {
        rval[0] = 0.0;
        if (nemoinpd(rexp,&dpar,1) != 1)
            error("setrange: parsing error %s",rexp);
	rval[1] = dpar;
    }
    rval[2] = rval[1] - rval[0];
}



parse_pairsi(key, pairs, nycol)
string key;
int *pairs;
int nycol;
{
    int n, i;

    n = nemoinpi(getparam(key),pairs,2*nycol);
    if (n>=0 && n<2) error("%s= needs at least 2 values",key);
    if (n<0) error("Parsing error %s=%s",key,getparam(key));
    for (i=n; i<2*nycol; i++) pairs[i] = pairs[i-2];

}

parse_pairsr(key, pairs, nycol)
string key;
real *pairs;
int nycol;
{
    int n, i;
    real r1, r2;

    n = nemoinpr(getparam(key),pairs,2*nycol);
    if (n>=0 && n<2) error("%s= needs at least 2 values",key);
    if (n<0) error("Parsing error %s=%s",key,getparam(key));
    for (i=n; i<2*nycol; i++) pairs[i] = pairs[i-2];
}


/*
 *  YAPP routine for PGPLOT
 */
#if 0
    
int pl_cursor(real *x, real *y, char *c)
{
    char inf[8], ans[8];
    int len, inf_len, ans_len;
    permanent float xsave, ysave;
    
    strcpy(inf,"CURSOR");
    inf_len = strlen(inf);
    ans_len = 1;
    pgqinf_(inf, ans, &len, inf_len, ans_len);
    if (ans[0]=='n' || ans[0]=='N') return 0; 
    pgcurs_(&xsave, &ysave, c, 1);
    *x = xsave;
    *y = ysave;
    return 1;
}
#endif  

