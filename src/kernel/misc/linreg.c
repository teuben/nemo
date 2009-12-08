
/*
 * NEMO history:
 *
 *	29-oct-90 1.0 	written				PJT
 *	21-feb-92 1.0a  usage, skip comment lines	PJT
 *                      parsing via nemoinpd() now
 *	24-feb-94 1.1   ansi header - moved stdinc.h up PJT
 *      22-sep-95 1.2a  allow xcol=, ycol= and fixed counting bug at maxpt PJT
 *      17-apr-04 1.2b  TAB's in the data were confused with numbers   PJT
 *                      blank lines were not handled properly either
 *       7-aug-09 1.2d  output of 6th method more friendly for parsers  PJT
 *
 */

 
/*
 *
 *                     SIX LINEAR REGRESSIONS 
 *     WRITTEN BY T. ISOBE, G. J. BABU AND E. D. FEIGELSON
 *               CENTER FOR SPACE RESEARCH, M.I.T.
 *                             AND
 *
 *              THE PENNSYLVANIA STATE UNIVERSITY
 *
 *                   REV. 1.0,   SEPTEMBER 1990
 *
 *       THIS SUBROUTINE PROVIDES LINEAR REGRESSION COEFFICIENTS 
 *    COMPUTED BY SIX DIFFERENT METHODS DESCRIBED IN ISOBE,
 *    FEIGELSON, AKRITAS, AND BABU 1990, ASTROPHYSICAL JOURNAL
 *    AND BABU AND FEIGELSON 1990, SUBM. TO TECHNOMETRICS.
 *    THE METHODS ARE OLS(Y/X), OLS(X/Y), OLS BISECTOR, ORTHOGONAL,
 *    REDUCED MAJOR AXIS, AND MEAN-OLS REGRESSIONS.
 *
 *    INPUT 
 *         X(I) : INDEPENDENT VARIABLE
 *         Y(I) : DEPENDENT VARIABLE
 *            N : NUMBER OF DATA POINTS 
 *
 *    OUTPUT
 *         A(J) : INTERCEPT COEFFICIENTS
 *         B(J) : SLOPE COEFFICIENTS
 *      SIGA(J) : STANDARD DEVIATIONS OF INTERCEPTS
 *      SIGB(J) : STANDARD DEVIATIONS OF SLOPES
 *     WHERE J = 1, 6. 
 *
 *    ERROR RETURNS
 *         CALCULATION IS STOPPED WHEN DIVISION BY ZERO WILL
 *         OCCUR (I.E. WHEN SXX, SXY, OR EITHER OF THE OLS
 *         SLOPES IS ZERO).
 */


#include <stdinc.h>

int linreg (int n,                                /* input:  number of data */
            real *x, real *y,                     /* input:  data points */
	    real *a, real *va, real *b, real *vb) /* output: fits */
{

/*
 *  local variables
 */
    double sxx, sxy, syy, sum1, sum2, sum3, rn;
    double xavg, yavg, sign, gm1, gm2, cov;
    int    i, j;

/*
 *  INITIALIZATIONS
 */
    xavg = yavg = 0.0;
    sxx  = syy  = sxy  = 0.0;
    sum1 = sum2 = sum3 = 0.0;
    rn = n;
    for (j=0; j<6; j++)
        va[j] = 0.0;

/*
 *  AVERAGE OF X AND Y arrays
 */
    for (i=0; i<n; i++) {
         xavg += x[i];
         yavg += y[i];
    }
    xavg /= rn;
    yavg /= rn;    

/*
 *  COMPUTE THE SLOPE COEFFICIENTS
 */
    for (i=0; i<n; i++) {       
        x[i] -= xavg;           /* reset data for stability */
	y[i] -= yavg;           /* reset data for stability */
        sxx += x[i]*x[i];
        syy += y[i]*y[i];
	sxy += x[i]*y[i];
    }
    if (sxx==0.0 || syy==0.0) {
        warning("SXX or SYY is zero - linreg has no solution");
        return(0);
    }

/*
 *  CHECK THE SIGN OF SXY
 */
    sign = (sxy < 0.0) ? -1.0 : 1.0;
 
    b[0] = sxy/sxx;
    b[1] = syy/sxy;
    b[2] = (b[0]*b[1] - 1.0 
              + sqrt((1.0 + sqr(b[0]))*(1.0 + sqr(b[1]))))/(b[0]+b[1]);
    b[3] = 0.5*(b[1] - 1.0/b[0]
               + sign*sqrt(4.0 + sqr(b[1] - 1.0/b[0])));
    b[4] = sign*sqrt(b[0]*b[1]);
    b[5] = 0.5*(b[0]+b[1]);

/*
 * COMPUTE INTERCEPT COEFFICIENTS
 */
    for (j=0; j<6; j++)
        a[j] = yavg - b[j]*xavg;

/*
 *     PREPARATIONS FOR COMPUTATIONS OF VARIANCES
 */
    if (b[0]==0.0 || b[2]==0.0) {
        warning("b[0] or b[1] zero; linreg has no solution");
        return(0);
    }
    gm1 = b[2]/((b[0] + b[1])
              * sqrt((1.0 + sqr(b[0]))*(1.0 + sqr(b[1]))));
    gm2 = b[3]/(sqrt(4.0*sqr(b[0]) + sqr(b[0]*b[1] - 1.0)));

    for (i=0; i<n; i++) {
        sum1 += sqr(x[i]*(y[i] - b[0]*x[i]));
        sum2 += sqr(y[i]*(y[i] - b[1]*x[i]));
        sum3 += x[i]*y[i]*(y[i] - b[0]*x[i])*(y[i] - b[1]*x[i]);
    }

/*
 *       COMPUTATION OF COVARIANCE OF B[0] AND B[1]
 */
    cov = sum3/(b[0]*sqr(sxx));

/*
 *    COMPUTATION OF THE STANDARD DEVIATIONS OF SLOPE COEFFICIENTS
 *       FIRST COMPUTE THE VARIANCES
 */
    vb[0] = sum1/sqr(sxx);
    vb[1] = sum2/sqr(sxy);
    vb[2] = sqr(gm1)*((sqr(1.0 + sqr(b[1])))*vb[0]
                  + 2.0*(1.0 + sqr(b[0]))*(1.0 + sqr(b[1]))*cov
                  + (sqr(1.0 +sqr(b[0])))*vb[1]);
    vb[3] = sqr(gm2)*(vb[0]/sqr(b[0]) + 2.0*cov
                 + sqr(b[0])*vb[1]);
    vb[4] = 0.25*(b[1]*vb[0]/b[0]
                    + 2.0*cov + b[0]*vb[1]/b[1]);
    vb[5] = 0.25*(vb[0] + 2.0*cov + vb[1]);

/*
 *   COMPUTATION OF THE STANDARD DEVIATIONS OF A(K)'S
 */
    for (i=0; i<n; i++) {
        va[0] += sqr((y[i] - b[0]*x[i])
                        *(1.0 - rn*xavg*x[i]/sxx));
        va[1] += sqr((y[i] - b[1]*x[i])
                        *(1.0 - rn*xavg*y[i]/sxy));
        va[2] += sqr((x[i]*(y[i] - b[0]*x[i])*(1.0 + sqr(b[1]))/sxx
                + y[i]*(y[i] - b[1]*x[i])*(1.0 + sqr(b[0]))/sxy)
                 *gm1*xavg*rn - y[i] + b[2]*x[i]);
        va[3] += sqr((x[i]*(y[i] - b[0]*x[i])/sxx
                 + y[i]*(y[i] - b[1]*x[i])*(sqr(b[0]))/sxy)*gm2
                   *xavg*rn/sqrt(sqr(b[0])) - y[i] + b[3]*x[i]);
        va[4] += sqr((x[i]*(y[i] - b[0]*x[i])*sqrt(b[1]/b[0])/sxx
                + y[i]*(y[i] - b[1]*x[i])*sqrt(b[0]/b[1])/sxy)
                  *0.5*rn*xavg - y[i] + b[4]*x[i]);
        va[5] += sqr((x[i]*(y[i] - b[0]*x[i])/sxx 
                + y[i]*(y[i] - b[1]*x[i])/sxy)
                  *0.5*rn*xavg - y[i] + b[5]*x[i]);
    }

/*
 *    CONVERT VARIANCES TO STANDARD DEVIATIONS
 */
    for (i=0; i<6; i++) {
        vb[i] = sqrt(vb[i]);
        va[i] = sqrt(va[i])/rn;
    }
/*
 *    RETURN DATA ARRAYS TO THEIR ORIGINAL STATE
 */
    for (i=0; i<n; i++) {
        x[i] += xavg;
	y[i] += yavg;
    }

    return 1;          /* return success code for C */
}
	

#if defined(TESTBED)
/*
 *                   SIX LINEAR REGRESSIONS 
 *                     WRITTEN BY T. ISOBE
 *                THE PENNSYLVANIA STATE UNIVERSITY
 *                            AND
 *    CENTER FOR SPACE RESEARCH, MASS. INSTITUTE OF TECHNOLOGY 
 *                    REV. 0.1   JULY 1990
 *			   1.0   Sept 1990
 *
 *       THIS PROGRAM PROVIDES LINEAR REGRESSION COEFFICIENTS 
 *    COMPUTED BY SIX DIFFERENT METHODS DESCRIBED IN ISOBE,
 *    FEIGELSON, AKRITAS, AND BABU 1990, ASTROPHYSICAL JOURNAL
 *    AND BABU AND FEIGELSON 1990, SUBM. TO TECHNOMETRICS.
 *    THE METHODS ARE OLS(Y/X), OLS(X/Y), OLS BISECTOR, ORTHOGONAL,
 *    REDUCED MAJOR AXIS, AND MEAN-OLS REGRESSIONS.
 *
 *    INPUT (FROM A FILE : FILE NAME FORMAT=<A9>
 *         X(I) : INDEPENDENT VARIABLE
 *         Y(I) : DEPENDENT VARIABLE
 *    WHERE I CAN BE UP TO 2000 AND THE FORMAT IS 2F10.3.
 *
 *    OUTPUT
 *         A(J) : INTERCEPT COEFFICIENTS
 *         B(J) : SLOPE COEFFICIENTS
 *        VA(J) : STANDARD DEVIATIONS OF INTERCEPTS
 *        VB(J) : STANDARD DEVIATIONS OF SLOPES
 *     WHERE J = 1, 6. 
 *
 *    OTHER
 *      FILE    :  INPUT FILE NAME (A9)
 *      NTOT    :  TOTAL NUMBER OF DATA POINTS
 *
 *    SUBROUTINE
 *    LINREG(NTOT,X,Y,A,VA,B,VB) : SIX LINEAR REGRESSION 
 *                                              SUBROUTINE.
 */

#include <getparam.h>

string defv[] = {
    "in=???\n       Input ascii table",
    "xcol=1\n       Column for X coordinates",
    "ycol=2\n       Column for Y coordinates",
    "mode=0\n       Output mode (0=all, 1..6=various methods",
    "maxline=10000\n Maximum size of columns",
    "VERSION=1.2e\n  27-aug-09 PJT",
    NULL,
};

string usage = "six linear regressions";

string cvsid="$Id$";

#define MAXCOL  256


patch_line(char *cp)      /* replace TAB and CR/NL with spaces */
{
  while (*cp) {
    if (*cp == '\t' || *cp == '\n' || *cp == '\r') *cp = ' ';
    cp++;
  }
}

nemo_main()
{
    stream fp;
    char   line[256];
    double *x, *y, a[6],va[6],b[6],vb[6], rn, data[MAXCOL];
    int    maxp, i, j, ntot, mode, xcol, ycol;

    mode = getiparam("mode");
    if (mode<0 || mode>6) {
        warning("Illegal mode, 0 assumed");
        mode = 0;
    }
    xcol = getiparam("xcol");
    ycol = getiparam("ycol");
    dprintf(1,"%s: Reading from columns %d and %d\n",getparam("in"),xcol,ycol);
/*
 *   OPEN FILE, READ THE DATA and CLOSE THE FILE
 */

    fp = stropen(getparam("in"),"r");
    maxp = getiparam("maxline");
    x = (double *) allocate(maxp * sizeof(double));
    y = (double *) allocate(maxp * sizeof(double));
    i=0;
    while (fgets(line,255,fp)) {
	if(line[0] == '#') continue;
	patch_line(line);
        if (i >= maxp) {
            printf("Maximum number of points (%d) reached\n",maxp);
            break;
        }
	if ((j=nemoinpd(line,data,MAXCOL)) < 0) {
	    warning("nemoinp error %d, skipping line: %s\n",j,line);
            continue;
        }
	if (j==0) continue;
        x[i] = data[xcol-1];
        y[i] = data[ycol-1];
        dprintf(2,"Data %d: %g %g\n",i,x[i],y[i]);
        i++;
    }
    ntot=i;
    fclose(fp);

/*
 *  CALL function 'linreg' IN WHICH SIX REGRESSIONS ARE PERFORMED.
 *      and PRINT OUT RESULTS when function performed well
 */

    if (linreg (ntot,x,y,a,va,b,vb)) {
        dprintf(1," Linear regression\n file: %s # of data: %d\n",
            getparam("in"),ntot);
        if (mode==0 || mode==1)    
            printf("OLS(Y/X): %g %g %g %g\n", a[0],va[0],b[0],vb[0]);
        if (mode==0 || mode==2)    
            printf("OLS(X/Y): %g %g %g %g\n", a[1],va[1],b[1],vb[1]);
        if (mode==0 || mode==3)    
            printf("BISECTOR: %g %g %g %g\n", a[2],va[2],b[2],vb[2]);
        if (mode==0 || mode==4)    
            printf("ORTHOG:   %g %g %g %g\n", a[3],va[3],b[3],vb[3]);
        if (mode==0 || mode==5)    
            printf("RED.AX.:  %g %g %g %g\n", a[4],va[4],b[4],vb[4]);
        if (mode==0 || mode==6)    
            printf("MEAN.OLS: %g %g %g %g\n", a[5],va[5],b[5],vb[5]);
    }
}


/*
 Sample DATA:
  1   11.0
  2   13.2
  3   14.5
  5   15.0
  6   16.5
  7   18.1
  8   19.2
  9   21.3
 RESULTS for Sample Data
  OLS(Y/X): 10.2851 .466391 1.13461 0.0765053
  OLS(X/Y): 10.0390 .462436 1.18263 0.0690194
  BISECTOR: 10.1635 .460563 1.15834 0.07211
  ORTHOG:   10.1452 .466985 1.16192 0.0730374
  RED.AX.:  10.1633 .460602 1.15837 0.0721114
  MEAN OLS: 10.1621 .460403 1.15862 0.0720172
 */
#endif /* TESTBED */

