/*
 * TABHIST: a general histogram plotter program for ascii data in tabular format
 *          
 *
 *	17-Mar-88  V1.0 :  created by P.J.Teuben
 *	15-apr-88  V1.1 :  added higher order moment, scaling in Y  PJT
 *	 1-jun-88  V2.0 :  just name changed nemohist->tabhist      PJT
 *				mallocation of x[] now done
 *			a: added labels along axes		PJT
 *	25-nov-88  V2.1 :  minmax calulation corrected	PJT
 *	26-sep-89  V2.2 :  tab= keyword added to ignore plot    PJT
 *	13-nov-90  V2.3 :  drange->nemoinpd PJT
 *	20-nov-90  V2.3a:  fixed out-of-range bug		PJT
 *	14-jun-91      b:  format statement needed more acc.	PJT
 *	26-mar-92  V2.4 :  stole crummy ascii histo from miriad PJT
 *	 2-nov-93  V2.6 :  use moments() routines - fix upper edge bug PJT
 *		      a :  use file_lines() and  nmax		pjt
 *	13-nov-93  V2.7 : add optional gaussian overlay plot 	pjt
 *	 4-mar-94     a : ansi header - fixed ylabel bug	pjt
 *	 3-aug-95     b : also report true minmax if min=max= used	pjt
 *	21-aug-95     c : fixed default X and Y labels on plot	pjt
 *	 1-sep-95     d : prototypes
 *	11-jul-96  V2.8 : log units now based 10, was base e	pjt
 *	20-nov-96     a : report median
 *      10-jan-96  V2.9 : optionally not show the residual      pjt
 *	20-feb-97     a : fixed for SINGLEPREC			pjt
 *      25-mar-97  V2.10: allow all columns to be used          pjt
 *                          ** unfinished **
 *	12-apr-97  V3.0:  allow cumulative histogram            pjt
 *	24-apr-97   3.0a: fix median calculation when limits set   pjt
 *	15-feb-99   3.1:  median can be turned off now		pjt
 *	22-dec-99   3.1a: fix reporting bug when 1 point read	pjt
 *       3-jun-01   3.2 : added nsigma=                         pjt
 * 
 * TODO:
 *     option to do dual-pass to subtract the mean before computing
 *     the higher order moments.
 */

/**************** INCLUDE FILES ********************************/ 

#include <stdinc.h>
#include <strlib.h>
#include <getparam.h>
#include <moment.h>
#include <yapp.h>
#include <axis.h>

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n                     Input file name",
    "xcol=1\n			  Column to use (or use: all)",
    "xmin=\n			  In case minimum used (need both minmax)",
    "xmax=\n			  In case maximum used (need both minmax)",
    "bins=16\n			  Number of bins",
    "maxcount=0\n		  Maximum along count-axis",
    "nmax=10000\n		  maximum number of data to be read if pipe",
    "ylog=f\n			  log scaling in Y?",
    "xlab=value\n	          Optional Label along X",
    "ylab=N\n			  Optional Label along Y",
    "headline=\n		  Optional headline in graph",
    "tab=f\n			  Table (t) or Plot( f)",
    "gauss=t\n			  Overlay gaussian plot?",
    "residual=t\n		  Overlay residual data-gauss(fit)",
    "cumul=f\n                    Override and do cumulative histogram instead",
    "median=t\n			  Compute median too (can be time consuming)",
    "nsigma=-1\n                  delete points more than nsigma",
    "VERSION=3.2a\n		  7-jun-01 PJT",
    NULL
};

string usage = "General tabular 1D statistics and histogram plotter";

/**************** SOME GLOBAL VARIABLES ************************/

#ifndef MAXHIST
#define MAXHIST	1024
#endif

local string input;			/* filename */
local stream instr;			/* input file */

local int col;				/* histogram column number */
local real   xrange[2];			/* range of  histogram values */
local int    nsteps;			/* number of divisions */

local real   *x = NULL;			/* pointer to array of nmax points */
local int    *ix = NULL;		/* pointer to index array */
local int    *iq = NULL;                /* boolean masking array */
local int    nmax;			/* lines to use at all */
local int    npt;			/* actual number of points */
local real   xmin, xmax;		/* actual min and max as computed */
local real   nsigma;                    /* outlier rejection attempt */
local bool   Qauto;			/* autoscale ? */
local bool   Qgauss;                    /* gaussian overlay ? */
local bool   Qresid;                    /* gaussian residual overlay ? */
local bool   Qtab;                      /* table output ? */
local bool   Qcumul;                    /* cumulative histogram ? */
local bool   Qmedian;			/* compute median also ? */
local int    maxcount;

local string headline;			/* text string above plot */
local string xlab, ylab;		/* text string along axes */
local bool   ylog;			/* count axis in logarithmic scale? */
local real   xplot[2],yplot[2];		/* borders of plot */

local real xtrans(real), ytrans(real);
local void setparams(void), read_data(void), histogram(void);



/****************************** START OF PROGRAM **********************/

nemo_main()
{
    setparams();			/* read the parameters */
    read_data();
    histogram();
}

local void setparams()
{
    input = getparam("in");
    if (streq(getparam("xcol"),"all"))
        col =  -1;
    else
        col = getiparam("xcol");

    if (col < 0) error("doing all columns is not implemented yet");
    
    nsteps=getiparam("bins");
    if (nsteps > MAXHIST) 
        error("bins=%d too large; MAXHIST=%d",nsteps,MAXHIST);
    if (hasvalue("xmin") && hasvalue("xmax")) {
	Qauto=FALSE;
    	xrange[0] = getdparam("xmin");
    	xrange[1] = getdparam("xmax");
	dprintf (2,"fixed plotrange %g : %g\n",xrange[0],xrange[1]);
	if (xrange[0] >= xrange[1]) error("Need xmin < xmax");
    } else {
	Qauto=TRUE;
	dprintf (2,"auto plotscaling\n");
    }
    maxcount=getiparam("maxcount");
    headline = getparam("headline");
    ylog=getbparam("ylog");
    xlab=getparam("xlab");
    ylab=getparam("ylab");
    Qgauss = getbparam("gauss");
    Qresid = getbparam("residual");
    Qtab = getbparam("tab");
    Qcumul = getbparam("cumul");
    if (Qcumul) {
        Qgauss=Qresid=FALSE;
        ylog=FALSE;
    }
    Qmedian = getbparam("median");
    if (ylog && streq(ylab,"N")) ylab = scopy("log(N)");

    nmax = file_lines(input,getiparam("nmax"));
    if (nmax<1) error("Problem reading from %s",input);
    x = (real *) allocate(sizeof(real)*(nmax+1));

    nsigma = getdparam("nsigma");

    instr = stropen (input,"r");
}



local void read_data()
{
    real *coldat[2];
    int   colnr[2], i;
		
    dprintf (2,"Reading datafile, col=%d\n",col);

    coldat[0] = x;  colnr[0] = col;     /* set to read 1 column */
    npt = get_atable(instr,1,colnr,coldat,nmax);        /* read it */
    if (npt == -nmax) {
    	warning("Could only read %d data",nmax);
    	npt = nmax;
    }
    minmax(npt, x, &xmin, &xmax);
    if (Qauto) {
	xrange[0] = xmin;
	xrange[1] = xmax;
    }
    /*  allocate index arrray , and compute sorted index for median */
    if (Qmedian) {
        ix = (int *) allocate(sizeof(int)*npt);
        sortptr(x,ix,npt);
    } 
}


local void histogram(void)
{
        int i,j,k, l, kmin, kmax, lcount;
	real count[MAXHIST], under, over;
	real xdat,ydat,xplt,yplt,dx,r,sum,sigma2, q, qmax;
	real mean, sigma, skew, kurt, lmin, lmax, median;
	Moment m;

	dprintf (0,"read %d values\n",npt);
	dprintf (0,"min and max value in column %d: [%g : %g]\n",col,xmin,xmax);
	if (!Qauto) {
	    xmin = xrange[0];
	    xmax = xrange[1];
	    dprintf (0,"min and max value reset to : [%g : %g]\n",xmin,xmax);
	    lmin = xmax;
	    lmax = xmin;
	    for (i=0; i<npt; i++) {
	    	if (x[i]>xmin && x[i]<=xmax) {
	    	    lmin = MIN(lmin, x[i]);
   	    	    lmax = MAX(lmax, x[i]);
	    	}
	    }
	    dprintf (0,"min and max value in range : [%g : %g]\n",lmin,lmax);
	} 

	for (k=0; k<nsteps; k++)
		count[k] = 0;		/* init histogram */
	under = over = 0;

	ini_moment(&m,4);
	for (i=0; i<npt; i++) {
		if (xmax != xmin)
		    k = (int) floor((x[i]-xmin)/(xmax-xmin)*nsteps);
		else
		    k = 0;
		dprintf(1,"%d k=%d %g\n",i,k,x[i]);
		if (k==nsteps && x[i]==xmax) k--;     /* include upper edge */
		if (k<0)       { under++; continue; }
		if (k>=nsteps) { over++;  continue; }
		count[k] = count[k] + 1;
		dprintf (4,"%d : %f %d\n",i,x[i],k);
		accum_moment(&m,x[i],1.0);
	}
	mean = mean_moment(&m);
	sigma = sigma_moment(&m);
	skew = skewness_moment(&m);
	kurt = kurtosis_moment(&m);

	if (nsigma > 0) {    /* remove outliers iteratively, starting from the largest */
	  iq = (int *) allocate(npt*sizeof(int));
	  for (i=0; i<npt; i++) {
#if 1
	    iq[i] = x[i] < xmin  || x[i] > xmax;
#else
	    iq[i] = 0;
#endif
	  }
	  lcount = 0;
	  do {               /* loop to remove outliers one by one */
	    qmax = -1.0;
	    for (i=0, l=-1; i<npt; i++) {     /* find largest deviation from current mean */
	      if (iq[i]) continue;            /* but skip previously flagged points */
	      q = (x[i]-mean)/sigma;
	      q = ABS(q);
	      if (q > qmax) {
		qmax = q;
		l = i;
	      }
	    }
	    if (qmax > nsigma) {
	      lcount++;
	      iq[l] = 1;
	      decr_moment(&m,x[l],1.0);
	      mean = mean_moment(&m);
	      sigma = sigma_moment(&m);
	      skew = skewness_moment(&m);
	      kurt = kurtosis_moment(&m);
	      dprintf(0,"%d/%d: removing point %d, m/s=%g %g qmax=%g\n",
		      lcount,npt,l,mean,sigma,qmax);
	      if (sigma <= 0) {
		/* RELATED TO presetting MINMAX */
		warning("BUG");
		accum_moment(&m,x[l],1.0);
		mean = mean_moment(&m);
		sigma = sigma_moment(&m);
		skew = skewness_moment(&m);
		kurt = kurtosis_moment(&m);
		dprintf(0,"%d/%d: LAST removing point %d, m/s=%g %g qmax=%g\n",
		      lcount,npt,l,mean,sigma,qmax);
		break;
		
	      }

	    } else
	      dprintf(0,"%d/%d: keeping point %d, m/s=%g %g qmax=%g\n",
		      lcount,npt,l,mean,sigma,qmax);

	    /* if (lcount > npt/2) break; */
	  } while (qmax > nsigma);
	  free(iq);
	}

	dprintf (0,"Number of points     : %d\n",n_moment(&m));
	dprintf (0,"Mean and dispersion  : %g %g\n",mean,sigma);
	dprintf (0,"Skewness and kurtosis: %g %g\n",skew,kurt);
	if (Qmedian) {
#if 0
        if (npt % 2) 
            median = x[ix[(npt-1)/2]];
        else
            median = 0.5 * (x[ix[npt/2]] + x[ix[npt/2-1]]);
        dprintf (0,"Median               : %g\n",median);
#else
	if (npt != n_moment(&m)) {
            kmin = kmax = -1;
            for (k=0; k<npt; k++) {
                if (x[ix[k]] >= xmin) {
                    kmin = k;
                    break;
                }   
            }
            for (k=npt-1; k>=0; k--) {
                if (x[ix[k]] <= xmax) {
                    kmax = k;
                    break;
		}   
            }
            if (kmin == -1)
                error("New median search: kmin=%d, no data in range?",kmin);
            if (kmax == -1)
                error("New median search: kmax=%d, no data in range?",kmax);
        } else {
            kmin = 0;
            kmax = npt-1;
        }
        if ((kmax-kmin+1)%2) {
                median = x[ix[kmin+(kmax-kmin)/2]];
        } else {
                median = 0.5 * (x[ix[kmin+(kmax-kmin+1)/2]] +
                                x[ix[kmin+(kmax-kmin+1)/2-1]]);

        }
        dprintf (0,"Median               : %g\n",median);
        dprintf (0,"Range for median     : %d %d \n",kmin,kmax);
#endif
        }

	dprintf (3,"Histogram values : \n");
        dx=(xmax-xmin)/nsteps;
	kmax=0;
	sum=0.0;
	for (k=0; k<nsteps; k++) {
	        sum = sum + dx*count[k];
		if (ylog) {
		    if (count[k]>0.0)
                        count[k] = log10(count[k]);
                    else
                        count[k] = -1.0;
		}
		if (count[k]>kmax)
			kmax=count[k];
		dprintf (3,"%f ",count[k]);
                if (Qcumul) {
                    if (k==0)
                        count[k] += under;
                    else
                        count[k] += count[k-1];
                }
	}
	dprintf (3,"\n");
        sigma2 = 2.0 * sigma * sigma;	/* gaussian */
        sum /= sigma * sqrt(2*PI);	/* scaling factor for equal area gauss */

	if (ylog && over>0)  over =  log10(over);
	if (ylog && under>0) under = log10(under);

	kmax *= 1.1;		/* add 10% */
        if (Qcumul) kmax = npt;
	if (maxcount>0)		/* force scaling by user ? */
		kmax=maxcount;	

	if (Qtab) {
            maxcount = 0;
            for (k=0; k<nsteps; k++)
                maxcount = MAX(maxcount,count[k]);
            if (maxcount>0)
                r = 29.0/maxcount;
            else
                r = 1.0;
            printf("  Bin    Value          Number\n");
            printf("       Underflow   %g\n",under);
            for (k=0; k<nsteps; k++) {
                j = (int) (r*count[k]) + 1;
                if (ylog) printf("%3d %13.6g %13.6g ", 
				k+1, xmin+(k+0.5)*dx, count[k]);
		else printf("%3d %13.6g %8d ", 
				k+1, xmin+(k+0.5)*dx, (int)count[k]);
                while (j-- > 0) printf("*");
                printf("\n");
            }
            printf("       Overflow    %g\n",over);
	    stop(0);
	}

#ifdef YAPP
	/*	PLOTTING */	
	plinit("***",0.0,20.0,0.0,20.0);

	xplot[0] = xmin;
	xplot[1] = xmax;
	yplot[0] = 0.0;
	yplot[1] = (real) kmax;
	xaxis (2.0, 2.0, 16.0, xplot, -7, xtrans, xlab);
	xaxis (2.0, 18.0,16.0, xplot, -7, xtrans, NULL);
	yaxis (2.0, 2.0, 16.0, yplot, -7, ytrans, ylab);
	yaxis (18.0, 2.0, 16.0, yplot, -7, ytrans, NULL);

        pljust(-1);     /* set to left just */
        pltext(input,2.0,18.2,0.32,0.0);             /* filename */
        pljust(1);
        pltext(headline,18.0,18.2,0.24,0.0);         /* headline */
        pljust(-1);     /* return to left just */

	xdat=xmin;
	dx=(xmax-xmin)/nsteps;
	plmove(xtrans(xmin),ytrans(0.0));
	for (k=0; k<nsteps; k++) {	/* nsteps= */
		xplt = xtrans(xdat);
		yplt = ytrans((real)count[k]);
		plline (xplt,yplt);
		xdat += dx;
		xplt = xtrans(xdat);
		plline (xplt,yplt);	
	}
	plline(xplt,ytrans(0.0));
	
	if (Qgauss) {                   /* plot model and residuals */
            if (ylog)
		plmove(xtrans(xmin),ytrans(-1.0));
            else
		plmove(xtrans(xmin),ytrans(0.0));
	    for (k=0; k<100; k++) {
	        xdat = xmin + (k+0.5)*(xmax-xmin)/100.0;
	        ydat = sum * exp( -sqr(xdat-mean)/sigma2);
	        if (ylog) ydat = log10(ydat);
	        plline(xtrans(xdat), ytrans(ydat));
	    }
        }

        if (Qresid) {

	    plltype(0,2);   /* dotted residuals */
	    xdat = xmin+0.5*dx;
            dprintf(1,"# residuals from gauss\n");
	    for (k=0; k<nsteps; k++, xdat +=dx) {
	        ydat = sum * exp( -sqr(xdat-mean)/sigma2);
	        dprintf(1,"%g %g %g\n",xdat,count[k],ydat);
	        if (ylog) ydat = log10(ydat);
                ydat = count[k] - ydat;
                if (k==0)
                    plmove(xtrans(xdat),ytrans(ydat));
                else
                    plline(xtrans(xdat),ytrans(ydat));
	    }
	    plltype(0,1);   /* back to normal line type */

	}
	plstop();
#endif
}

local real xtrans(real x)
{
	return (2.0 + 16.0*(x-xplot[0])/(xplot[1]-xplot[0]));
}

local real ytrans(real y)
{
	return (2.0 + 16.0*(y-yplot[0])/(yplot[1]-yplot[0]));
}
