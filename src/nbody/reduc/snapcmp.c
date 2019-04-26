/*
 * SNAPCMP.C: compare two N-body snapshots.
 *
 *	22-jan-90   V1.2					JEB
 *	25-nov-90   V1.3 helpvec, allocate from library 	PJT
 *	23-mar-94   V1.3a
 *	29-mar-94   V1.4 process all snapshots - added time to output  pjt
 *	27-jan-00   V1.4b  more usefule error messages                 pjt
 *       9-apr-01       c  added header in output   pjt
 *      14-apr-04   V1.5  added headline= to plots                     pjt
 *      15-apr-05         add RMS to output
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>
#include <vectmath.h>
#include <moment.h>
#include <yapp.h>
#include <axis.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {			/* DEFAULT INPUT PARAMETERS         */
    "in1=???\n				  first input file name",
    "in2=???\n				  second input file name",
    "time1=all\n			  time selected from 1st file",
    "time2=all\n			  time selected from 2nd file",
    "obs=r\n				  observable used in comparison",
    "diffpart=true\n			  difference particles first",
    "relative=false\n			  scale difference by ref value",
    "log=f\n                              use logarithm ",
#if HISTOGRAM
    "nbins=8\n				  number of bins in histogram",
    "xrange=0.0:1.0\n			  range of values to histogram",
    "xlabel=\n				  label for X axis; default is obs ",
    "ymax=\n				  if given, use for max N value",
    "ylabel=N\n				  label for Y axis",
#endif
#if SCATTERPLOT
    "xvar=\n				  indp var to plot; default is obs",
    "xrange=0.0:2.0\n			  range of values to display",
    "xlabel=\n				  X axis label; default is xvar",
    "yrange=0.0:1.0\n			  range of delta values to display",
    "ylabel=delta\n			  label for Y axis",
#endif
#if HISTOGRAM || SCATTERPLOT
    "nxticks=3\n			  number of tickmarks along X axis",
    "nyticks=3\n			  number of tickmarks along Y axis",
    "xbox=2.5:17.5\n			  extent of graph in x direction",
    "ybox=5.0:15.0\n			  extent of graph in y direction",
    "formal=false\n			  if true, make publication plot",
    "headline=\n                          header",
#endif
    "VERSION=1.6a\n			  25-apr-2019 PJT",
    NULL,
};

string usage = "compare two N-body snapshots";

string cvsid = "$Id$";

rproc obsfunc;
bool diffpart;
bool relative;

#if HISTOGRAM || SCATTERPLOT

real xrange[2], yrange[2];
real xbox[2], ybox[2];
real xtrans(real), ytrans(real);
string headline;

#endif

local real small_dt = 1.0e-14;

local bool Qlog;



local real *snapcmp(Body*, Body*, int, real);
local int cmpreal(real*, real*);
local void printquart(real*, int, real);

extern rproc btrtrans(string);

void nemo_main()
{
    stream instr1, instr2;
    string time1, time2, file1, file2;
    Body *btab1 = NULL, *btab2 = NULL;
    int nbody1, bits1, nbody2, bits2;
    real dt, tsnap1 = 0.0, tsnap2 = 0.0, *result;

    file1 = getparam("in1");
    file2 = getparam("in2");
    instr1 = stropen(file1, "r");
    get_history(instr1);
    instr2 = stropen(file2, "r");
    get_history(instr2);
    time1 = getparam("time1");
    time2 = getparam("time2");
    obsfunc = btrtrans(getparam("obs"));
    diffpart = getbparam("diffpart");
    relative = getbparam("relative");
    Qlog = getbparam("log");
#if HISTOGRAM || SCATTERPLOT
    headline = getparam("headline");
#endif

    for(;;) {
      do {
	if (!get_snap_by_t(instr1, &btab1, &nbody1, &tsnap1, &bits1, time1))
		break;
        dprintf(1,"snap1: reading time=%g\n",tsnap1);
      } while (bits1 == TimeBit);
      do {
	if (!get_snap_by_t(instr2, &btab2, &nbody2, &tsnap2, &bits2, time2))
		break;
        dprintf(1,"snap2: reading time=%g\n",tsnap2);
      } while (bits2 == TimeBit);
      if (bits1 == 0 && bits2 == 0)
        break;
      else if (bits1 == 0)
    	error("no more snapshots found in %s at t=%g",file1,tsnap2);
      else if (bits2 == 0)
    	error("no more snapshots found in %s at t=%g",file2,tsnap1);
      if (nbody1 != nbody2) 
	error("%s = %d, %d different", NobjTag, nbody1, nbody2);
      if (tsnap1 != tsnap2) {
	dt = ABS(tsnap2-tsnap1);
	if (dt > small_dt)
	  warning("times = %f, %f are different (%g)", tsnap1, tsnap2, dt);
      }
      if (bits1 != bits2)
	warning("bits = 0x%x, 0x%x are different", bits1, bits2);
      result = snapcmp(btab1, btab2, nbody1, tsnap1);
#if STATISTICS
      printquart(result, nbody1, tsnap1);
#endif
#if HISTOGRAM
      histogramplot(result, nbody1);
#endif
#if SCATTERPLOT
      snapcmpplot(result, btab1, nbody1, tsnap1);
#endif
    }
}

local real *snapcmp(Body *btab1, Body *btab2, int nbody, real tsnap)
{
    real *result, oref, odiff;
    Body *bp1, *bp2, dbody;
    int i;

    result = (real *) allocate(nbody * sizeof(real));
    for (i = 0; i < nbody; i++) {
	bp1 = btab1 + i;
	bp2 = btab2 + i;
	oref = (*obsfunc)(bp1, tsnap, i);
	if (diffpart) {
	    Mass(&dbody) = Mass(bp1) - Mass(bp2);
	    SUBV(Pos(&dbody), Pos(bp1), Pos(bp2));
	    SUBV(Vel(&dbody), Vel(bp1), Vel(bp2));
	    Phi(&dbody) = Phi(bp1) - Phi(bp2);
	    SUBV(Acc(&dbody), Acc(bp1), Acc(bp2));
	    Aux(&dbody) = Aux(bp1) - Aux(bp2);
	    Key(&dbody) = Key(bp1) - Key(bp2);
	    odiff = (*obsfunc)(&dbody, tsnap, i);
	} else
	    odiff = oref - (*obsfunc)(bp2, tsnap, i);
	if (relative)
	    odiff = odiff / oref;
	result[i] = odiff;
    }
    return result;
}

real show(real x) {
  if (!Qlog) return x;
  if (x==0.0) return -99.99;
  return log10(x);
}
#if STATISTICS

local void printquart(real result[], int nbody, real tsnap)
{
  static int Qheader = 1;
  int i;
  Moment m;

  if (Qheader) {
    printf("# time  Min  Qlow Median Qhigh  Max   Mean Sigma\n");
    printf("# obs=%s\n",getparam("obs"));
    Qheader = 0;
  }
  qsort(result, nbody, sizeof(real), cmpreal);
  ini_moment(&m,2,0);
  for (i=0; i<nbody; i++)
    accum_moment(&m,result[i],1.0);
  printf("%g   %g %g %g %g %g  %g %g\n",
	 tsnap,
	 show(result[0]),
	 show(result[nbody/4]),
	 show(result[nbody/2-1]),
	 show(result[3*nbody/4-1]),
	 show(result[nbody-1]),
	 show(mean_moment(&m)),
	 show(sigma_moment(&m)));
	 
}

/* should have prototype :: int (*compar)(const void *, const void *)) */

local int cmpreal(real *ap, real *bp)
{
    return (*ap < *bp ? -1 : *ap > *bp ? 1 : 0);
}





#endif

#if HISTOGRAM

histogramplot(result, nbody)
real result[];
int nbody;
{
    int nbins, *hgram, *histogram(), i;
    string xlabel, ylabel;
    real x, y;

    nbins = getiparam("nbins");
    setrange(xrange, getparam("xrange"));
    hgram = histogram(result, nbody, nbins, xrange);
    fitrange(yrange, getiparam("nyticks"), hgram, nbins);
    plinit("", 0.0, 20.0, 0.0, 20.0);
    xlabel = getparam(streq(getparam("xlabel"), "") ? "obs" : "xlabel");
    ylabel = getparam("ylabel");
    axisplot(xlabel, ylabel);
    pltext(headline,18.0,18.2,0.24,0.0);
    for (i = 1; i <= nbins; i++) {
	x = xrange[0] + (i - 1.0) * (xrange[1] - xrange[0]) / nbins;
	y = hgram[i];
	if (i == 1)
	    plmove(xtrans(x), ytrans(y));
	else
	    plline(xtrans(x), ytrans(y));
	x = xrange[0] + i * (xrange[1] - xrange[0]) / nbins;
	plline(xtrans(x), ytrans(y));
    }
    plstop();
}

int *histogram(values, nvals, nbins, range)
real values[];
int nvals;
int nbins;
real range[];
{
    int *hgram, i, j;
    real x;

    hgram = (int *) allocate((nbins + 2) * sizeof(int));
    for (j = 0; j <= nbins+1; j++)
	hgram[j] = 0;
    for (i = 0; i < nvals; i++) {
	x = (values[i] - range[0]) / (range[1] - range[0]);
	if (0.0 <= x && x < 1.0)
	    j = 1 + floor(nbins * x);
	else
	    j = (x < 0.0 ? 0 : nbins + 1);
	hgram[j]++;
    }
    return (hgram);
}

fitrange(range, nyticks, hgram, nbins)
real range[];
int nyticks;
int hgram[];
int nbins;
{
    int hmax, i, hstep;

    if (streq(getparam("ymax"), "")) {		/* derive ymax from hgram[] */
	hmax = 0;
	for (i = 1; i <= nbins; i++)
	    hmax = MAX(hmax, hgram[i]);
	hstep = 5;
	while ((nyticks + 0.5) * hstep < hmax)
	    hstep = hstep + 5;
	range[0] = 0.0;
	range[1] = (nyticks + 1) * hstep;
    } else {					/* use value supplied       */
	range[0] = 0.0;
	range[1] = getdparam("ymax");
    }
}

#endif

#if SCATTERPLOT

snapcmpplot(result, btab, nbody, tsnap)
real result[];
Body *btab;
int nbody;
real tsnap;
{
    rproc indfunc;
    string xlabel, ylabel;
    int nxticks, nyticks, i;
    real x, y;
    Body *bp;

    if (! streq(getparam("xvar"), "")) {
	indfunc = btrtrans(getparam("xvar"));
	xlabel = getparam(streq(getparam("xlabel"), "") ? "xvar" : "xlabel");
    } else {
	indfunc = obsfunc;
	xlabel = getparam(streq(getparam("xlabel"), "") ? "obs" : "xlabel");
    }
    setrange(xrange, getparam("xrange"));
    setrange(yrange, getparam("yrange"));
    ylabel = getparam("ylabel");
    plinit("", 0.0, 20.0, 0.0, 20.0);
    axisplot(xlabel, ylabel);
    pltext(headline,18.0,18.2,0.24,0.0);
    for (i = 0, bp = btab; i < nbody; i++, bp++) {
	x = xtrans((*indfunc)(bp, tsnap, i));
	y = ytrans(result[i]);
	if (xbox[0] <= x && x <= xbox[1] &&
	      ybox[0] <= y && y <= ybox[1])
	    plpoint(x, y);
    }
    plstop();
}

#endif

#if HISTOGRAM || SCATTERPLOT

axisplot(xlabel, ylabel)
string xlabel, ylabel;
{
    int nxticks, nyticks;

    setrange(xbox, getparam("xbox"));
    setrange(ybox, getparam("ybox"));
    nxticks = getiparam("nxticks");
    nyticks = getiparam("nyticks");
    if (getbparam("formal")) {
	formalaxis = TRUE;
	xaxisvar.labdn = 0.44;
	xaxisvar.szlab = 0.40;
	yaxisvar.numdn = 0.24;
	yaxisvar.labdn = 0.24;
	yaxisvar.szlab = 0.40;
    }
    xaxis(xbox[0], ybox[0], xbox[1] - xbox[0],
	  xrange, - nxticks, xtrans, xlabel);
    xaxis(xbox[0], ybox[1], xbox[1] - xbox[0],
	  xrange, - nxticks, xtrans, NULL);
    yaxis(xbox[0], ybox[0], ybox[1] - ybox[0],
	  yrange, - nyticks, ytrans, ylabel);
    yaxis(xbox[1], ybox[0], ybox[1] - ybox[0],
	  yrange, - nyticks, ytrans, NULL);
}

setrange(rval, rexp)
real rval[];
string rexp;
{
    char *cptr;

    cptr = strchr(rexp, ':');
    if (cptr != NULL) {
        rval[0] = atof(rexp);
	rval[1] = atof(cptr+1);
    } else {
        rval[0] = 0.0;
	rval[1] = atof(rexp);
    }
}

real xtrans(real x)
{
    return (xbox[0] + (xbox[1] - xbox[0]) *
	      (x - xrange[0]) / (xrange[1] - xrange[0]));
}

real ytrans(real y)
{
    return (ybox[0] + (ybox[1] - ybox[0]) *
	      (y - yrange[0]) / (yrange[1] - yrange[0]));
}

#endif
