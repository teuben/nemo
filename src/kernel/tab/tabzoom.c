/*
 *  TABZOOM: zoom around in a table (tabplot and tabview style)
 *           and define actions with points.
 *
 *  cloned from tabview with pl_cursor() calls
 *
 *  Note: yapp_pgplot needs to be compiled with the -DCOLOR option
 *	  for tabview.
 *
 *	8-oct-92  V0.0 .prototype.			PJT
 *      1-dec-92  V0.2 added | and ! commands in interactive mode PJT
 *	4-mar-94   0.3 using Moment (+ansi headers)		  pjt
 *                     and turbo charged the program a bit (may be
 *                     somewhat pgplot dependant ???
 *     18-may-96   0.4 polygonic or nearest point flagging
 *     21-mar-97   0.4a  SINGLEPREC fix
 *     26-jan-00   0.4b  fix code for new pl_cursor in yapp_pgplot
 *
 *      1-jul-03   1.0 cloned from tabview --- or should we stick to tabview --
 *     21-oct-03   1.2 added xlab= ylab=
 */

#include <stdinc.h>
#include <getparam.h>
#include <strlib.h>
#include <yapp.h>
#include <axis.h>
#include <layout.h>
#include <moment.h>

string defv[] = {
    "in=???\n       Input table filename",
    "col=\n         Columns to be named and extracted (def: all)",
    "xvar=%1\n      X-variable to plot (%1=first extracted column)",
    "yvar=%2\n      Y-variable to plot",
    "color=1\n      Color of displayed points (1,2,3,...)",
    "psize=0\n      Size of displayed points (cm)",
    "ptype=0\n      Type of displayed points (not used)",
    "xrange=\n      Range in X [autoscale]",
    "yrange=\n      Range in Y [autoscale]",
    "xlab=\n        Fix the label along X axis to this",
    "ylab=\n        Fix the label along Y axis to this",
    "layout=\n      Filename with layout to be added to plot?",
    "action=\n      (optional) Clicked Point action",
    "out=\n         Filename for output table",
    "options=i,x,y,data\n   Output options",
    "maxstat=1\n    Max amount of statistics printed out",
    "cursor=t\n     Go into cursor mode directly?",
    "headline=\n    Headline for plot",
    "VERSION=1.2b\n 2-feb-05 PJT",
    NULL,
};

string usage="dynamic query table viewer with interactive zoom and action";

typedef struct point {
  real x,y;		/* points in plot coordinates */
  real xd,yd;		/* points in world coordinates */
  int pcolor;		/* color; 0=invisib 1=visib 2,3,4...=color */
  real psize;           /* size: must be in cm (yapp) */
  int ptype;            /* type: -- not used -- */
  bool visib;		/* currently visible ? */
  bool oldvis;          /* previously visible ? */
  real *par;    	/* pointer to original parameters with this point */
  string spar;          /* -same- but now the original textual info */
} point;

typedef struct slider {
    char *name;
    real smin, smax;
    real slo, shi;
    real step;
    int  stepper;
} slider;

#ifndef MAXSLIDERS
# define MAXSLIDERS 32       /* maximum number of columns to be processed */
#endif

#ifndef MAXPOINTS
# define MAXPOINTS MAXLINES  /* maximum number of points (from stdinc.h) */
#endif

#ifndef MAXLINLEN
# define MAXLINLEN  132      /* maximum linelength in a table */
#endif

#ifndef MAXPOL
# define MAXPOL  64          /* maximum # segments in a polygon */
#endif



Moment mom[MAXSLIDERS];
slider sliders[MAXSLIDERS];
point  points[MAXPOINTS]; 
int nsliders, npoints;
string outname;
string outoptions;

int maxstat;

real gxmin, gxmax, gymin, gymax;
real xmin, xmax, ymin, ymax, xplot[2], yplot[2];
string xname, yname, xlab, ylab, headline;
plcommand *layout;
string action;
bool Qcursor;

/* extern */

extern string *burststring(string,string);

/* local functions */
int zoom(void);
real xtrans(real), ytrans(real), ixtrans(real), iytrans(real);
int pl_cursor(real*, real*, char *);
int pl_cursor1(real, real, real*, real*, char *);
void show_lr(int,int);
void show_mr2(int);
real pearson(int,int);
void do_action(int, string);


nemo_main()
{
    int s;

    setparams();
    gettable();
    ini_display();
    for (;;) {
        s = interact();
        if (s < 0) break;
        re_display(0);
    }
    stop_display();
}

setparams()
{
    if (hasvalue("out"))
        outname = getparam("out");
    else
        outname = NULL;
    outoptions = getparam("options");
    if (hasvalue("layout"))
        layout = pl_fread(getparam("layout"));
    else
        layout = NULL;
    maxstat = getiparam("maxstat");
    action = getparam("action");
    Qcursor = getbparam("cursor");
    headline = getparam("headline");
}

gettable()
{
    char line[MAXLINLEN], name[30];
    stream instr;
    string  fname, *colnames, *sp;
    bool first=TRUE;
    int i, j, k, n, nret, one=1, n_invisib=0;
    real errval=0.0, range[2], drange, dmean, sum, sum1, sum2, sum3;
    real pcolor, psize, cmin, cmax, ptype, meanj, sigj, meank, sigk;

    fname = getparam("in");
    instr = stropen(fname,"r");

    colnames = burststring(getparam("col"),", ");
    nsliders = xstrlen(colnames,sizeof(string)) - 1;
    if (nsliders>MAXSLIDERS) error("Too many sliders");
    for (j=0; j<nsliders; j++)
        sliders[j].name = colnames[j];
    if (nsliders==0) warning("Don't use col= yet, doesn't work yet");

    while (fgets(line,MAXLINLEN,instr) != NULL) {     /* read all lines */
        n = strlen(line);
        if (line[0]=='#' || line[0]==';' || line[0]=='!')  
            continue;                               /* skip comment lines */
        if (line[n-1]=='\n') line[n-1]='\0';      /* patch (new)line */
        sp = burststring(line,", \t");
        n = xstrlen(sp,sizeof(string))-1;       /* number of items found */
        if (first) {                            /* perhaps first time around */
            first = FALSE;			/* prevent re-entry */
            if (nsliders==0) {                  /* need some extra inits ? */
                nsliders = n;
                for (j=0; j<nsliders; j++) {    /* create default names */
                    sprintf(name,"%d",j+1);
                    sliders[j].name = scopy(name);
                }
            }
            npoints = 0;
            for (j=0; j<nsliders; j++)
            	ini_moment(&mom[j],4,0);

        }
        if (npoints==MAXPOINTS) {
            warning("Full table with %d values; no more read",npoints);
            break;
        }
	points[npoints].spar = scopy(line);
        points[npoints].par = (real *) allocate(nsliders*sizeof(real));
        for (j=0; j<nsliders; j++) {        /* parse all numbers */
            nret = nemoinpr(sp[j],&points[npoints].par[j], 1);
            if (nret != 1)
		warning("gettable: (%d,%d): %s",npoints+1,j+1,sp[j]);
	    else
	    	accum_moment(&mom[j],points[npoints].par[j],1.0);

        }
        freestrings(sp);
        npoints++;
    }
    dprintf(0,"Reading table %s with %d columns and %d rows\n",
	fname,nsliders,npoints);
    strclose(instr);


    printf("\n### Overview statistics:\n");
    printf("%-10s: %4s %8s %8s %8s %8s %8s %8s\n",
        "Slider","Npts","min","max","mean","sigma","skewness","kurtosis");
    for (j=0; j<nsliders; j++) {  /* set remaining sliders info + output */
        sliders[j].slo = sliders[j].smin = min_moment(&mom[j]);
        sliders[j].shi = sliders[j].smax = max_moment(&mom[j]);
    	printf("%-10s: %4d %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n",
    		sliders[j].name, n_moment(&mom[j]),
    		min_moment(&mom[j]), max_moment(&mom[j]),
    		mean_moment(&mom[j]), sigma_moment(&mom[j]),
    		skewness_moment(&mom[j]), kurtosis_moment(&mom[j]));
    }

    if (maxstat>0) {
        printf("\n### 'r' Pearson correlation matrix:\n");
        printf("   : ");
        for (j=0; j<nsliders; j++) printf(" %5d",j+1);
        printf("\n");
        for (j=0; j<nsliders; j++) {
            printf("%2d : ",j+1);
            for (k=0; k<j; k++)
                printf(" %5.2f",pearson(j,k));
       	    printf("   *  \n");
        }
    }

    if (maxstat>1) {

        printf("\nLinear regression:\n");
        for (j=0; j<nsliders; j++) {
            for (i=0; i<nsliders; i++)
                show_lr(i,j);
            printf("\n");
        }
    }
    if (maxstat>2) {
        for (j=0; j<nsliders; j++)
            show_mr2(j);
    }


    xname = getparam("xvar");    
    xlab = hasvalue("xlab") ? getparam("xlab") : xname;
    inifie(xname);                     /* parse the X expression */
    for (i=0; i<npoints; i++) {         /* set X for all points */
        dofie(points[i].par, &one, &points[i].x, &errval);
        if (i==0) {
            xmin = points[i].x;
            xmax = points[i].x;
        } else {
            xmin = MIN(points[i].x, xmin);
            xmax = MAX(points[i].x, xmax);
        } 
    }
    if (hasvalue("xrange")) {
        if (nemoinpr(getparam("xrange"),range,2) != 2)
            error("Need 2 values for xrange");
        xmin = range[0];
        xmax = range[1];
    } else {
        if (xmin==xmax) error("No range in X data");
        drange = xmax-xmin;
        dmean = 0.5*(xmin+xmax);
        xmin = dmean - 0.5*drange*1.05;
        xmax = dmean + 0.5*drange*1.05;
    }
    dprintf(0,"Xvar %s: displayed min= %g max= %g\n",xname,xmin,xmax);

    yname = getparam("yvar");    
    ylab = hasvalue("ylab") ? getparam("ylab") : yname;
    inifie(yname);                       /* parse the Y expression */
    for (i=0; i<npoints; i++) {           /* set Y for all points */
        dofie(points[i].par, &one, &points[i].y, &errval);
        if (i==0) {
            ymin = points[i].y;
            ymax = points[i].y;
        } else {
            ymin = MIN(points[i].y, ymin);
            ymax = MAX(points[i].y, ymax);
        } 
    }
    if (hasvalue("yrange")) {
        if (nemoinpr(getparam("yrange"),range,2) != 2)
            error("Need 2 values for xrange");
        ymin = range[0];
        ymax = range[1];
    } else {
        if (ymin==ymax) error("No range in Y data");
        drange = ymax-ymin;
        dmean = 0.5*(ymin+ymax);
        ymin = dmean - 0.5*drange*1.05;
        ymax = dmean + 0.5*drange*1.05;
    }
    dprintf(0,"Yvar %s: displayed min= %g max= %g\n",yname,ymin,ymax);
    gxmin = xmin;
    gxmax = xmax;
    gymin = ymin;
    gymax = ymax;

    /* careful: the color stuff doesn't work */
    inifie(getparam("color"));
    for (i=0; i<npoints; i++) {
        dofie(points[i].par, &one, &pcolor, &errval);
        points[i].pcolor = pcolor;
        if (i==0) {
            cmin = pcolor;
            cmax = pcolor;
        } else {
            cmin = MIN(cmin, pcolor);
            cmax = MAX(cmax, pcolor);
        }
        if(points[i].pcolor < 1) {
            points[i].pcolor = 0;
            n_invisib++;
        }
    }
    dprintf(0,"Color : min= %g max= %g\n",cmin,cmax);
    if (n_invisib) warning("%d points will never be visible",n_invisib);

    inifie(getparam("psize"));
    for (i=0; i<npoints; i++) {
        dofie(points[i].par, &one, &psize, &errval);
        points[i].psize = psize;
    }

    inifie(getparam("ptype"));
    for (i=0; i<npoints; i++) {
        dofie(points[i].par, &one, &ptype, &errval);
        points[i].ptype = ptype;
    }
    
    for (i=0; i<npoints; i++) {                 /* DEBUG OUTPUT */
        dprintf(2,"(x,y) %g %g (sliders)",points[i].x, points[i].y);
        for (j=0; j<nsliders; j++) {
            dprintf(2," %g",points[i].par[j]);
        }
        dprintf(2,"\n");
    }

} /* gettable */


int zoom(void)
{
  int mode, posn;
  char c;
  float xref, yref, x, y, tmp;

  /* assume it has cursor mode */
  dprintf(1,"Old box: X: %g %g Y: %g %g\n", xplot[0], xplot[1], yplot[0], yplot[1]);

  mode = 0;
  posn = 0;
  pgband_(&mode, &posn, &xref, &yref, &x, &y, &c, 1);
  if (c=='X') return 0;
  mode = 2;
  xref = x;
  yref = y;
  pgband_(&mode, &posn, &xref, &yref, &x, &y, &c, 1);
  dprintf(0,"Box: %g %g %g %g Char: %c\n", xref, yref, x, y, c);
  if (c=='X') return 0;

  if (x < xref) {  tmp=x; x=xref; xref=tmp; }
  if (y < yref) {  tmp=y; y=yref; yref=tmp; }

  xmin = ixtrans((double)xref);
  xmax = ixtrans((double)x);
  ymin = iytrans((double)yref);
  ymax = iytrans((double)y);
  dprintf(0,"New min/max: X: %g %g Y: %g %g\n", xmin,xmax, ymin,ymax);

  ini_display();
  dprintf(0,"New box: X: %g %g Y: %g %g\n", xplot[0], xplot[1], yplot[0], yplot[1]);
  return 1;
}


/*
 *  Interact with user; return -1 if done with interacting
 *                              0 if all sliders needs updated
 *                             >0 if slider (interact()-1) needs updated
 */

string interact_help = "Menu of commands   (* requires cursor interaction ADX=left/middle/right):\n\n\
    l <num>     modify 'lo' slider\n\
    h <num>     modify 'hi' slider\n\
    l s <step>  set step in 'lo' slider\n\
    h s <step>  set step in 'hi' slider\n\
    b <step>    step both 'lo' and 'hi' slider\n\
    <digits>    change to this slider # to interact with\n\
    c         * cursor mode (NEW)\n\
    z         * zoom with a box (NEW)\n\
    Z           zoom reset to original setting\n\
    f         * cursor flagging (EXPERIMENTAL)\n\
    r           reset lo/hi to min/max for this slider\n\
    s           show min/lo/hi/max for all sliders\n\
    u		update screen new\n\
    o <file>    output table, w/ optional override out= filename\n\
    x           swap lo and hi (invert logic)\n\
    q           quit\n\
    !cmd        execute a shell command 'cmd'\n\
    |cmd        pipe visible data as ascii table to 'cmd'\n\
    <RETURN>    stepper mode\n\
    ?           this help plus status";

interact()
{
  permanent int s=-1;  /* the currently selected slider */
  permanent char cmd_system[MAXLINLEN], cmd_popen[MAXLINLEN]; /* for ! and | */
  permanent bool first = TRUE;
  FILE *fp;
  bool done=FALSE;
  char *word(), *show_slider(), line[MAXLINLEN], cmd[MAXLINLEN];
  int i, itmp, imax;
  real rtmp; 
  real xcur, ycur, dmax, dmaxold, d2;
  real xcm, ycm, xp, yp;
  float xpol[MAXPOL], ypol[MAXPOL];
  char c;

  if (s<0) {  /* first time entry: initialize but return to get display up */
    s=0;            /* set activer slider to 0 */
    strcpy(cmd_system,"");
    strcpy(cmd_popen,"");
    return s;       /* and return code to update all sliders */
  }
  
  do {		/* loop this while 'done==FALSE' */
    printf("%s",show_slider(s));
    
    /* gets() is dangerous and should not be used */
    if (Qcursor && first) {
      strcpy(line,"c");
      first = FALSE;
    } else
      if (gets(line) == NULL) return -1;
    
    switch (line[0]) {
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':  
      itmp = atoi(line)-1; 
      if (itmp>=0 && itmp<nsliders)
	s = itmp;
      else
	warning("Not a legal slider");
      break;
    case '!':  
      if (line[1]) strcpy(cmd_system, &line[1]);
      if (cmd_system[0]==0) break;
      printf("%s\n",cmd_system);
      system(cmd_system);
      break;
    case '|': 
      if (line[1]) strcpy(cmd_popen,&line[1]);
      if (cmd_popen[0]==0) break;
      printf("%s\n",cmd_popen);
      fp = popen(cmd_popen,"w");
      if (fp==NULL) break;
      funny_table_output(fp);
      pclose(fp);
      break;
    case 'u':  
      re_display(0);
      break;
    case 'l':  
      if (*word(line,2) == 's') {
	sliders[s].step = atof(word(line,3)); 
	sliders[s].stepper = 1;
      } else {
	sliders[s].slo = atof(word(line,2)); 
	sliders[s].stepper = 0;
	done = TRUE;
      }
      break;
    case 'h': 
      if (*word(line,2) == 's') {
	sliders[s].step = atof(word(line,3)); 
	sliders[s].stepper = 2;
      } else {
	sliders[s].shi = atof(word(line,2)); 
	sliders[s].stepper = 0;
	done = TRUE;
      }
      break;
    case 'b':  
      sliders[s].step = atof(word(line,2));
      sliders[s].stepper = 3;
      break;
    case 's':  
      printf("\n");
      for (i=0; i<nsliders; i++) 
	printf("%s\n",show_slider(i));
      printf("\n");
      break;
    case 'o':  
      table_output(word(line,2));
      break;
    case 'r':
      sliders[s].shi = sliders[s].smax;
      sliders[s].slo = sliders[s].smin;
      done = TRUE;
      break;
    case '?':  
      printf("%s\n\n",interact_help);
      printf("Xvar=%s Yvar=%s\n",xname,yname);
      printf("!cmd:   %s\n",cmd_system);
      printf("|cmd:   %s\n",cmd_popen);
      printf("action: %s\n",action);
      break;
    case 'q':  
      return -1;
    case 'x':  
      rtmp = sliders[s].slo;
      sliders[s].slo = sliders[s].shi;
      sliders[s].shi = rtmp;
      done = TRUE;
      break;
    case 'f':  
      while(pl_cursor(&xcm, &ycm, &c)) {
	if (c=='X') break;
	if (c=='D') {         /* try and flag point */
	  if (xcm > 2 && xcm < 18 &&
	      ycm > 2 && ycm < 18) {
	    xp = ixtrans(xcm);
	    yp = iytrans(ycm);
	    dprintf(1,"POINT: %g %g\n",xp,yp);
	  }
	} else if (c=='A') {
	  pl_getpoly(xpol, ypol, MAXPOL);
	}
      }
    case 'c':      /* interactive mode */
      for (;;) {                      /* loop forever until done */
	pl_cursor(&xcm, &ycm, &c);
	if (c=='?') {
	  printf("Cursor keys:\n");
	  printf("q   quit\n");
	  printf("A   (left mouse)   action\n");
	  printf("D   (middle mouse) test where the action is\n");
	  printf("X   (right mouse)  zoom\n");
	  continue;
	} else if (c=='q') 
	  break;
#if 0
	if (c == 'X') break;
	if (c == 'z') {
	  while(zoom())
	    re_display(0);
	}
#else
	if (c == 'X') {
	  dprintf(0,"Going into zoom mode\n");
	  zoom();
	  re_display(0);
	  continue;
	}
#endif
	if (xcm > 2 && xcm < 18 &&
	    ycm > 2 && ycm < 18) {
	  xp = ixtrans(xcm);
	  yp = iytrans(ycm);
	  dmax = 999.999;
	  for (i=0; i<npoints; i++) {
	    if (!points[i].visib) continue;
	    d2 = sqr(points[i].x - xcm) + sqr(points[i].y - ycm);
	    if (d2 < dmax) {
	      dmaxold = dmax;
	      dmax = d2;
	      imax = i;
	    }
	  }
	  if (dmax/dmaxold < 0.1) {
	    dprintf(0,"Found point %d\n",imax);
	    do_action((c=='A'), points[imax].spar);
	  } else
	    warning("Not near enuf to point %d",imax);
	}
      }
      break;
    case 'z':  
      while (zoom())
	re_display(0);
      break;
    case 'Z':  
      xmin = gxmin;
      xmax = gxmax;
      ymin = gymin;
      ymax = gymax;
      ini_display();
      re_display(0);
      break;
    case '\0':     /* this is when you hit the RETURN key */
      if (sliders[s].stepper == 1) {
	sliders[s].slo += sliders[s].step;
	done = TRUE;
      } else if (sliders[s].stepper == 2) {
	sliders[s].shi += sliders[s].step;
	done = TRUE;
      } else if (sliders[s].stepper == 3) {
	sliders[s].slo += sliders[s].step;
	sliders[s].shi += sliders[s].step;
	done = TRUE;
      } else
	warning("Not in stepper mode");
      break;
    default:
      warning("%s: Not a valid command",line);
      
    } /* switch(line[0]) */
  } while (!done);
  return s+1;
}

char *word(char *s,int n)
{
    string *sp;
    permanent char wordn[40];

    sp = burststring(s," \t");
    if (xstrlen(sp,sizeof(string))-1 >= n) {
        strcpy(wordn,sp[n-1]);
        freestrings(sp);
    } else
        strcpy(wordn,"");
    return wordn;
}

char *show_slider(int s)
{
    permanent char line[MAXLINLEN];

    sprintf(line, "Slider %s [%g %g %g %g]: ", sliders[s].name, 
        sliders[s].smin, sliders[s].slo, sliders[s].shi, sliders[s].smax);
    return line;
}

table_output(string name)
{
    int i, j, nvis=0;
    stream ostr;
    bool scanopt(), Qi, Qx, Qy, Qd;
    char *show_slider(), *fname;

    if (name!=NULL && *name != 0)
        fname = name;
    else if (outname == NULL) {
        warning("No out= was ever supplied");
        return;
    } else
        fname = outname;
    Qi = scanopt(outoptions,"i");
    Qx = scanopt(outoptions,"x");
    Qy = scanopt(outoptions,"y");
    Qd = scanopt(outoptions,"data");
    if (!Qi && !Qx && !Qy && !Qd) {
        warning("No valid options=: select any of i,x,y,data");
        return;
    }

    ostr = stropen(fname,"w+");        /* force overwrite */
    fprintf(ostr,"# Output created by %s\n",getargv0());
    fprintf(ostr,"# Ncols = %d\n",nsliders);
    fprintf(ostr,"# Nrows = %d\n",npoints);
    for (j=0; j<nsliders; j++)
        fprintf(ostr,"# %s\n",show_slider(j));
    fprintf(ostr,"# Output = ");
    if (Qi) fprintf(ostr," i");
    if (Qx) fprintf(ostr," x");
    if (Qy) fprintf(ostr," y");
    if (Qd) fprintf(ostr," data");
    fprintf(ostr,"\n");
    for (i=0; i<npoints; i++) {
        if (!points[i].visib) continue;
        if (Qi) fprintf(ostr,"%d ",i+1);
        if (Qx) fprintf(ostr,"%g ",points[i].xd);
        if (Qy) fprintf(ostr,"%g ",points[i].yd);
        if (Qd) for (j=0; j<nsliders; j++)
                    fprintf(ostr,"%g ",points[i].par[j]);
        fprintf(ostr,"\n");
        nvis++;
    }
    strclose(ostr);
    dprintf(0,"[%d data written to \"%s\"]\n",nvis,fname);

}

funny_table_output(FILE *fp)
{
    int i, nvis=0;

    for (i=0; i<npoints; i++) {
        if (!points[i].visib) continue;
        fprintf(fp,"%22.16g %22.16g\n",points[i].xd, points[i].yd);
        nvis++;
    }
    dprintf(0,"[%d data piped]\n",nvis);

}



/*
 *  Display engine 
 */

ini_display()
{
    int i;
    static bool first = TRUE;

    plinit("***",0.0,20.0,0.0,20.0);

    xplot[0] = xmin;        /* set scales for xtrans() */
    xplot[1] = xmax;
    yplot[0] = ymin;        /* set scales for ytrans() */
    yplot[1] = ymax;

    if (first) {
      dprintf(1,"Initializing points\n");
      dprintf(1,"0: %g %g %g %g\n",points[0].x,points[0].y,points[0].xd,points[0].yd);
      for (i=0; i<npoints; i++) {     /* rescale all points to cm */
        points[i].xd = points[i].x;             /* xd,yd in physical */
        points[i].yd = points[i].y;
        points[i].x = xtrans(points[i].x);      /* x,y now in cm (plot) */
        points[i].y = ytrans(points[i].y);
        points[i].visib = TRUE;
        points[i].oldvis = TRUE;
      }
      dprintf(0,"1: %g %g %g %g\n",points[0].x,points[0].y,points[0].xd,points[0].yd);
    } else {
      plframe();
      dprintf(0,"Resetting points\n");
      dprintf(0,"1: %g %g %g %g\n",points[0].x,points[0].y,points[0].xd,points[0].yd);
      for (i=0; i<npoints; i++) { 
        points[i].x = xtrans(points[i].xd);    /* put x,y back in cm */
        points[i].y = ytrans(points[i].yd);
      }
      dprintf(0,"1: %g %g %g %g\n",points[0].x,points[0].y,points[0].xd,points[0].yd);
    }
    box_display();
    if (layout) pl_exec(layout);
    first = FALSE;
}

box_display()
{
    plcolor(1);
    /* should have a more intelligent way to figure out #nidigits needed
     */
    xaxvar(3,-1,-1,-1,-1);
    yaxvar(3,-1,-1,-1,-1);

    xaxis ( 2.0,  2.0, 16.0, xplot, -7, xtrans, xlab);
    xaxis ( 2.0, 18.0, 16.0, xplot, -7, xtrans, NULL);
    yaxis ( 2.0,  2.0, 16.0, yplot, -7, ytrans, ylab);
    yaxis (18.0,  2.0, 16.0, yplot, -7, ytrans, NULL);
    if (*headline) {
      pljust(1);
      pltext(headline,18.0,18.2,0.24,0.0); 
      pljust(-1);
    }
}


stop_display()
{
    plstop();
}

re_display(int k)       /* redisplay all, or slider 'k' (1..nsliders) */
{
    permanent int color=-1;	/* remember current paint color */
    int i, j, nplot, jhi, jlo, p_color, nzoom;
    real s, t1, t2, slo, shi;
    bool vis, Qerase;
    slider *sp;
    point *pp;
	
    nplot = 0;  /* counter of points plotted this turn */
    nzoom = 0;  /* counter of points within zoomed box */
    Qerase = TRUE;	/* can this yapp erase points or not ? */
    if (k==0) {
       jlo = 0;   jhi = nsliders;	/* do all */
    } else {
       jlo = k-1; jhi=k;		/* use only one slider k=1..nsliders */
    }

    t1 = 60000.0*cputime();
    for (i=0; i<npoints; i++) {			/* loop over all npoints */
        for (j=jlo; j<jhi; j++) {		/* check these sliders */
            s = points[i].par[j];
            slo = sliders[j].slo;
            shi = sliders[j].shi;
	    if (slo < shi)          /* if LO slider left of HI slider */
	        vis = (slo <= s && s <= shi);
	    else
	        vis = (s >= slo || s <= shi);
	    if(!vis) break;  /* if already invisible, break the sliders loop */
	}
	if (vis) nplot++;
	
        if (points[i].visib == vis && !vis) continue;	/* invisible */
        
        if (Qerase) {
            p_color = vis ? points[i].pcolor : 0;
            if (color != p_color) { 
                color = p_color; 
                plcolor(color); 
            }
            dprintf(1,"%d %d %g %g\n",i+1,color,points[i].xd, points[i].yd);
#if 0
            plpoint(points[i].x, points[i].y);
#else
	    switch (points[i].ptype) {
	    case 0:
	      plpoint(points[i].x, points[i].y);
	      break;
	    case 1:
	      plcircle(points[i].x, points[i].y, points[i].psize);
	      break;
	    case 2:
	      plcross(points[i].x, points[i].y, points[i].psize);
	      break;
	    case 3:
	      plbox(points[i].x, points[i].y, points[i].psize);
	      break;
	    default:
	      error("Bad point type");
	    }
#endif
        }
      	points[i].oldvis = points[i].visib;
	points[i].visib  = vis;
	if (2 < points[i].x && points[i].x < 18 &&
	    2 < points[i].y && points[i].y < 18) nzoom++;
    }

    t2 = 60000.0*cputime();
    dprintf(1,"CPU select = %g ms\n",t2-t1);
    t1=t2;
		
    if (!Qerase) { /* erase screen, draw box/layout/points */
#if 0
        /* erase all points: brute force method */
        plframe();
        box_display();
#else
        /* erase all points: overwriting visible ones */
        plcolor(0);        
        for (i=0; i<npoints; i++)
            plpoint(points[i].x, points[i].y);
        plcolor(color);

#endif        
        for (i=0; i<npoints; i++)
            if (points[i].visib) {
                p_color = points[i].pcolor;
                if (color != p_color) {
                    color = p_color;
                    plcolor(color);
                }
                plpoint(points[i].x, points[i].y);
	    }
    }
    t2 = 60000.0*cputime();
    dprintf(1,"CPU display = %g ms\n",t2-t1);
    t1=t2;
    dprintf(0,"DISPLAY: %d/%d/%d points visible\n",nzoom,nplot,npoints);
    plflush();
}	

real xtrans(real x)
{
    return (2.0 + 16.0*(x-xplot[0])/(xplot[1]-xplot[0]));
}

real ytrans(real y)
{
    return (2.0 + 16.0*(y-yplot[0])/(yplot[1]-yplot[0]));
}

real ixtrans(real x)
{
        return xplot[0] + (x-2.0)*(xplot[1]-xplot[0])/16;
}

real iytrans(real y)
{
        return yplot[0] + (y-2.0)*(yplot[1]-yplot[0])/16;
}

/* 
 *  YAPP routine for PGPLOT
 */

#if 0
/* this is now moved into yapp_pgplot.c */
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

int pl_cursor1(real xa, real ya, real *x, real *y, char *c)
{
    char inf[8], ans[8];
    int mode, posn, len, inf_len, ans_len;
    float xref, yref, xreq, yreq;

    strcpy(inf,"CURSOR");
    inf_len = strlen(inf);
    ans_len = 1;
    pgqinf_(inf, ans, &len, inf_len, ans_len);
    if (ans[0]=='n' || ans[0]=='N') return 0;

    mode = 1;
    posn = 1;
    xref = xreq = xa;
    yref = yreq = ya;
    pgband_(&mode, &posn, &xref, &yref, &xreq, &yreq, c);
    *x = xreq;
    *y = yreq;
    return 1;
}




void show_lr(int i, int j)
{
    int k;
    real x, y, s, sx, sy, sxx, syy, sxy, d;

    s = sx = sy = sxx = syy = sxy = 0;
    
    for (k=0; k<npoints; k++) {
        x = points[k].par[i];
        y = points[k].par[j];

        s += 1.0;     sx += x;     sy += y;
        sxx += x*x;   sxy += x*y;  syy += y*y;
    }
    d = s*sxx - sx*sx;
    
    printf("%d = a + b * %d: a= %g +/- %g  b= %g +/- %g rab= %g\n",
        i+1, j+1, (sxx*sy-sx*sxy)/d, sxx/d,
              (s*sxy-sx*sy)/d, s/d,
              -sx/sqrt(s*sxx));
    

}
/*
 *  Show 2D multiple regression, assuming variable 'j' is the
 *  dependant variable. It selects all other combinations
 *  (l,m) as the independant variables and computes the
 *  regression coeffient R^2
 */
 
void show_mr2(int jdep)
{
    real ry1, ry2, r12;
    int i,j;
    
    printf("\n###'R^2' correlation matrix for variable %d:\n",jdep+1);
    printf("   : ");
    for (j=0; j<nsliders; j++) printf(" %5d",j+1);
    printf("\n");

    for (j=0; j<nsliders; j++) {
    	printf("%2d : ",j+1);
    	for (i=0; i<j; i++) {
    	    if (i==jdep || j==jdep)
    	        printf("   *  ");
    	    else {
                ry1 = pearson(jdep,i);
                ry2 = pearson(jdep,j);
                r12 = pearson(i,j);
      	        printf(" %5.2f",(ry1*ry1+ry2*ry2-2*ry1*ry2*r12)/(1-r12*r12));
      	    }
    	}
        printf("   *  \n");
    }
}

real pearson(int i, int j)
{
    real sum1, sum2, sum3, meani, meanj, sigi, sigj;
    int k;
    
    meani = mean_moment(&mom[i]);
    sigi = sigma_moment(&mom[i]);
    meanj = mean_moment(&mom[j]);
    sigj = sigma_moment(&mom[j]);
    sum1=sum2=sum3=0;
    for (k=0; k<npoints; k++) {
  	sum1 += (points[k].par[i]-meani)*(points[k].par[j]-meanj);
    	sum2 += sqr(points[k].par[i]-meani);
        sum3 += sqr(points[k].par[j]-meanj);
    }
    return sum1/sqrt(sum2*sum3);
}

void do_action(int run, string s)
{
  char cmd[MAXLINLEN];

  if (*action) {
    dprintf(run ? 1 : 0,"ACTION(%s): %s\n",action,s);
    sprintf(cmd,"%s %s",action,s);
    if (run) system(cmd);
  } else
    dprintf(1,"noACTION(): %s\n",s);
}
