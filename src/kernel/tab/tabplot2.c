/*
 * TABPLOT2:    plot two tables, with the aim to compare
 *
 *      15-feb-2024 V0.1:    cloned off tabplot
 */

/**************** INCLUDE FILES ********************************/ 
 
#include <stdinc.h>	
#include <getparam.h>
#include <yapp.h>
#include <axis.h>
#include <layout.h>
#include <table.h>
#include <mdarray.h>
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
    "in1=???\n           Input first table",
    "in2=???\n           Input second table",
    "xcol=1\n		 x coordinate column",
    "ycol=2\n		 y coordinate column",
    "xmin=\n		 X-min value if not autoscaled",
    "xmax=\n		 X-max value if not autoscaled",
    "ymin=\n		 Y-min value if not autoscaled",
    "ymax=\n		 Y-max value if not autoscaled",
    "point=1,0.1\n       point type and size pairs for each table (yapp)",
    "line=0,0\n          line width and style pairs for each table (yapp)",
    "color=\n		 colors for the symbols/lines for each table",
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
    "headline=\n	 headline in graph on top right",
    "fullscale=f\n       Use full autoscale in one axis if other autoscaled?",
    "cursor=\n           Optional output file to retrieve cursor coordinates",
    "backtrack=t\n       Allow backtrack in line= option",
    "layout=\n           Optional input layout file",
    "first=f\n           Layout first or last?",
    "readline=f\n        Interactively reading commands",
    "pyplot=\n           Template python plotting script",
    "VERSION=0.2\n	 15-feb-2024 PJT",
    NULL
};

string usage = "compare two tables";


/**************** GLOBAL VARIABLES ************************/

local string input1, input2;			/* filename */
local stream instr1, instr2;			/* input file */
local table *tptr1, *tptr2;                     /* table */
local mdarray2 d1, d2;                          /* data[col][row] */

local int xcol[MAXCOL], ycol[MAXCOL], nxcol, nycol;	/* column numbers */
local real xrange[2], yrange[2];		/* range of values */
local int npcol;                                /* MAX(nxcol,nycol) */

local real  *x[MAXCOL], *y[MAXCOL]; 			/* data from file */
local int    npt1;				/* actual number of data points */
local int    npt2;				/* actual number of data points */

local real xmin, xmax;			/* actual min and max as computed */
local real ymin, ymax;
local bool   Qautox0, Qautoy0, Qfull;		/* autoscale ? */
local bool   Qautox1, Qautoy1;

local string headline;			/* text string above plot */
local string xlab, ylab;			/* text string along axes */
local real xplot[2],yplot[2];		        /* borders of plot */
local int    yapp_line[2*MAXCOL];            /* thickness and pattern of line */
local int    yapp_color[MAXCOL];	      /* color per column */
local real yapp_point[2*MAXCOL];             /* what sort of point to plot */

local real xcoord[MAXCOORD], ycoord[MAXCOORD];/* coordinate lines */
local int nxcoord, nycoord;
local int nxticks, nyticks;                   /*& number of tickmarks */
local int ncolors;
local real xbox[3], ybox[3];
local real xscale, yscale;

local plcommand *layout;
local bool layout_first;
local bool Qreadlines;
local bool Qbacktrack;

void setparams(void);
void read_data(void);
void plot_data(void);
void plot_points(int np, real xp[], real yp[], real pstyle,
		 real psize, int lwidth, int lstyle, int color);
void parse_pairsi(string key, int *pairs, int nycol);
void parse_pairsr(string key, real *pairs, int nycol);
double my_sqrt(double x);

extern real smedian(int, real *);
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
      pyplot_plot2(pstr, input1, input2, xcol, ycol, xrange, yrange);
      pyplot_close(pstr);
    }
      
    read_data();
    plot_data();
}

void setparams(void)
{
    int  i;
   
    input1 = getparam("in1");             /* input table file */
    input2 = getparam("in2");
    instr1 = stropen (input1,"r");        /* input stream */
    instr2 = stropen (input2,"r");
    tptr1 = table_open(instr1,0);         /* table */
    tptr2 = table_open(instr2,0);
    
    nxcol = nemoinpi(getparam("xcol"),xcol,MAXCOL);  // should be 1 or 2 
    nycol = nemoinpi(getparam("ycol"),ycol,MAXCOL);  // should be 1 or 2 
    if (nxcol < 1) error("Error parsing xcol=%s",getparam("xcol"));
    if (nycol < 1) error("Error parsing ycol=%s",getparam("ycol"));
    if (nxcol==1) {
      xcol[1]=xcol[0];
      nxcol=2;
    }
    if (nycol==1) {
      ycol[1]=ycol[0];
      nycol=2;
    }
    npcol = MAX(nxcol,nycol);
    dprintf(1,"xcol: %d %d    ycol: %d %d\n",xcol[0],xcol[1],ycol[0],ycol[1]);

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
       
    Qbacktrack = getbparam("backtrack");
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

void read_data(void)
{
    int i, j, colnr[1+MAXCOL];
		
    dprintf (2,"Reading datafile, xcol,ycol=%d..,%d,...\n",xcol[0],ycol[0]);

    colnr[0] = xcol[0];
    colnr[1] = ycol[0];
    d1 = table_md2cr(tptr1, 2,colnr,0,0);
    colnr[0] = xcol[1];
    colnr[1] = ycol[1];
    d2 = table_md2cr(tptr2, 2,colnr,0,0);

    x[0] = &d1[0][0];  // x[0] has npt1
    x[1] = &d2[0][0];  // x[1] has npt2

    y[0] = &d1[1][0];  // npt1
    y[1] = &d2[1][0];  // npt2

    npt1 = table_nrows(tptr1);
    npt2 = table_nrows(tptr2);

    // @todo
    if (xscale != 1.0) {
      for (j=0; j<nxcol; j++)
	for (i=0; i<npt1; i++)
	  x[j][i] *= xscale;
    }
    if (yscale != 1.0) {
      for (j=0; j<nycol; j++)
	for (i=0; i<npt1; i++)
	  y[j][i] *= yscale;
    }

    xmin = ymin =  HUGE;
    xmax = ymax = -HUGE;
    for (i=0; i<npt1; i++) {     /* find global min and max in X and all Y */
      xmax=MAX(x[0][i],xmax);
      xmin=MIN(x[0][i],xmin);
      ymax=MAX(y[0][i],ymax);
      ymin=MIN(y[0][i],ymin);
    }
    for (i=0; i<npt2; i++) {     /* find global min and max in X and all Y */
      xmax=MAX(x[1][i],xmax);
      xmin=MIN(x[1][i],xmin);
      ymax=MAX(y[1][i],ymax);
      ymin=MIN(y[1][i],ymin);
    }
}


void plot_data(void)
{
    int    i, j;
    real xcur, ycur, edge;
    char c;
    stream cstr;
    bool first;

    if (npt1<1 || npt2<1) {
        warning("Nothing to plot, npt1=%d npt2=%d",npt1,npt2);
        return;
    }
    dprintf (0,"read %d points from 1\n",npt1);
    dprintf (0,"read %d points from 2\n",npt2);
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
      for (i=0; i<npt1; i++) {     /* go through data, find X min and max again */
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
        for (i=0; i<npt1; i++) {     /* go through data, find X min and max again */
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
    xaxis(xbox[0],ybox[0],xbox[2], xplot, -nxticks, xtrans, xlab);
    xaxis(xbox[0],ybox[1],xbox[2], xplot, -nxticks, xtrans, NULL);
    yaxis(xbox[0],ybox[0],ybox[2], yplot, -nyticks, ytrans, ylab);
    yaxis(xbox[1],ybox[0],ybox[2], yplot, -nyticks, ytrans, NULL);
    for (i=0; i<nxcoord; i++) {
        plmove(xtrans(xcoord[i]),ytrans(yplot[0]));
        plline(xtrans(xcoord[i]),ytrans(yplot[1]));
    }
    for (i=0; i<nycoord; i++) {
        plmove(xtrans(xplot[0]),ytrans(ycoord[i]));
        plline(xtrans(xplot[1]),ytrans(ycoord[i]));
    }

    pljust(-1);     /* set to left just */
    pltext("2 files",xbox[0],ybox[1]+0.2,0.32,0.0);             /* filename */
    pljust(1);
    pltext(headline,xbox[1],ybox[1]+0.2,0.32,0.0);          /* headline */
    pljust(-1);     /* return to left just */

    plot_points( npt1, x[0], y[0], 
		 yapp_point[0], yapp_point[1],
		 yapp_line[0], yapp_line[1],
		 ncolors > 0 ? yapp_color[0] : -1);
    plot_points( npt2, x[1], y[1], 
		 yapp_point[2], yapp_point[3],
		 yapp_line[2], yapp_line[3],
		 ncolors > 0 ? yapp_color[1] : -1);

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

/* PLOT_POINTS
 *
 *      Plot datapoints, optionally connecting them with lines,
 *
 *      X,Y Points with value NaN are skipped
 */
 
void plot_points (int np, real *xp, real *yp, 
		  real pstyle, real psize, int lwidth, int lstyle, int color)
{
    int i, ipstyle;
    dprintf(1,"plot_points: %d  %g %g %d %d %d\n", np, pstyle, psize, lwidth, lstyle, color);
	    

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
            while (++i < np) {
	      if (xp[i] != NaN) {
		if (Qbacktrack || (xtrans(xp[i]) > xtrans(xp[i-1])))
                    plline (xtrans(xp[i]), ytrans(yp[i]));  /* draw line */
	      }
	      plmove (xtrans(xp[i]), ytrans(yp[i]));      /* move to point */	    
	    }

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



void parse_pairsi(string key, int *pairs, int nycol)
{
    int n, i;

    n = nemoinpi(getparam(key),pairs,2*nycol);
    if (n>=0 && n<2) error("%s= needs at least 2 values",key);
    if (n<0) error("Parsing error %s=%s",key,getparam(key));
    for (i=n; i<2*nycol; i++) pairs[i] = pairs[i-2];

}

void parse_pairsr(string key, real *pairs, int nycol)
{
    int n, i;

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

