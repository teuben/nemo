/*  ft_open, ft_intpol, ft_close
 *
 * Function Table Processor: open an ascii table with the purpose
 *      of designating two columns as a function lookup table, with
 *      either a linear or spline interpolation function evaluator
 *      An otherwise irrelevant structure is passed around as control
 *      block; at closing all allocated space is freed again.
 *
 * Example Usage:
 *      ftp = ft_open(fname, mode, xcol, ycol);
 *      value = ft_intpol(ftp, x0);
 *      errors = ft_close(ftp);
 *
 *          4-oct-93    Created             peter teuben
 *          3-nov-93    Inverse function? (TOOLBOX only)  pjt
 *	    6-nov-93    fixed incomplete pipe data problem  pjt
 *          5-mar-94    use MAXLINES from header
 *     	    4-sep-97    1.1  added linear, made it default (spline has bug) PJT
 *         20-jun-01	gcc 3
 *	   23-sep-01    nemo_file_line
 *	    8-dec-01    MAX_LINES now in maxsizes.h
 */



#include <stdinc.h>
#include <strlib.h>
#include <funtab.h>
#include <spline.h>

extern int nemo_file_lines(string, int);
extern int get_atable(stream, int, int*, real **, int);

FunctionTable *ft_open(string fname, int mode, int xcol, int ycol)
{
    FunctionTable *ftp;
    stream instr;
    real *coldat[2];   /* 2, for X and Y columns */
    int n, colnr[2];

    ftp = (FunctionTable *) allocate(sizeof(FunctionTable));
    ftp->name = scopy(fname);
    ftp->mode  = mode;
    
    n = nemo_file_lines(fname,MAXLINES);  /* initial guess for number of points */
    if (n < 1) {
        warning("ft_open: No data can be read from %s",fname);
        return NULL;
    }
    /* allocate a bit conservatively */
    coldat[0] = ftp->x = (real *) allocate(n * sizeof(real));
    coldat[1] = ftp->y = (real *) allocate(n * sizeof(real));
    colnr[0] = xcol;
    colnr[1] = ycol;

    instr = stropen(fname,"r");
    ftp->n = get_atable(instr, 2, colnr, coldat, n);
    strclose(instr);
    
    if (ftp->n < 0) {
    	ftp->n = -ftp->n;
    	warning("Could only read %d data from file",ftp->n);
    }

    if (ftp->n < 1) {
        warning("ft_open: no data has been read from %s",fname);
        return NULL;
    }
    if (mode == FUNTAB_SPLINE) {
        ftp->coeff = (real *) allocate(3*ftp->n*sizeof(real));
        spline(ftp->coeff, ftp->x, ftp->y, ftp->n);        
    } else if (mode == FUNTAB_LINEAR) {
        ftp->coeff = NULL;
    } else {
        error("ft_open: wrong mode (%d)",mode);
    }

    ftp->errors = 0;

    return ftp;
}
        
    
real ft_spline(FunctionTable *ftp, real x0)
{
    if (ftp==NULL) {
        warning("ft_spline: not initialized (ftp=NULL)");
        return 0.0;
    } 

    if (ftp->n == 1)
        return ftp->y[0];
    else if (ftp->n > 1)
        return seval(x0, ftp->x, ftp->y, ftp->coeff, ftp->n);
    else
        error("ft_spline: bad call, n=%d",ftp->n);

    return 0.0;
}

real ft_linear(FunctionTable *ftp, real x0)
{
    int i, n = ftp->n;
    real val;

    if (ftp==NULL) {
        warning("ft_linear: not initialized (ftp=NULL)");
        return 0.0;
    } 

    if (n == 1)
        return ftp->y[0];
    else if (n > 1) {
        if (x0 < ftp->x[0]) {
            warning("ft_linear: %g too small",x0);
            return ftp->y[0];
        } if (x0 > ftp->x[n-1]) {
            warning("ft_linear: %g too big",x0);
            return ftp->y[n-1];
        } else {
            for (i=0; i<n-1; i++) {       /* assumes sorted !!! */
                if (ftp->x[i] <= x0 && x0 < ftp->x[i+1]) {
                    val = ftp->y[i] + (x0 - ftp->x[i]) *
                          (ftp->y[i+1]-ftp->y[i])/(ftp->x[i+1]-ftp->x[i]);
                    return val;
                }
            }
            return ftp->y[n-1];
        }
    } else
        error("ft_linear: bad call, n=%d",ftp->n);

    return 0.0;
}

int ft_close(FunctionTable *ftp)
{
    int errors;
    
    if (ftp==NULL) {
        warning("ft_close: not initialized (ftp=NULL)");
        return 0.0;
    }
    errors = ftp->errors; 
    if (errors != 0)
        warning("FunctionTable[%s]: %d errors during processing",
                    ftp->name, ftp->errors);

    if (ftp->x)     free(ftp->x);
    if (ftp->y)     free(ftp->y);
    if (ftp->coeff) free(ftp->coeff);
    free(ftp->name);
    free(ftp);

    return errors;
}

#ifdef TOOLBOX

#ifndef MAX_LINES  
#define MAX_LINES   10000
#endif


#include <getparam.h>

string defv[] = {
    "in=???\n       Input table with Y=f(X)",
    "xcol=1\n       Column for X",
    "ycol=2\n       Column for Y",
    "x=\n           X values to lookup (Y to return)",
    "y=\n           Y values to lookup (X to return)",
    "format=%g\n    Output format",
    "mode=linear\n  Lookup mode (spline, linear, near)",
    "VERSION=1.2\n  8-dec-01 PJT",
    NULL,
};

string usage="Function Table lookup";

nemo_main()
{
    FunctionTable *ftp;
    string fmt1, fmode = getparam("mode");
    char fmt[100];
    real x[MAX_LINES], y;
    bool Qinv;
    int i, n, mode=FUNTAB_SPLINE;

    if (hasvalue("x") && hasvalue("y")) error("Can only handle x= or y=");

    Qinv = hasvalue("y");
    if (streq(fmode,"linear")) mode = FUNTAB_LINEAR;

    if (Qinv) {
        n = nemoinpd(getparam("y"),x,MAX_LINES);
        if (n<0) error("Error # %d parsing y=; or too many?",n);
    } else {
        n = nemoinpd(getparam("x"),x,MAX_LINES);
        if (n<0) error("Error # %d parsing x=; or too many?",n);
    }
    fmt1 = getparam("format");
    sprintf(fmt,"%s %s\n",fmt1,fmt1);
    dprintf(1,"Using format=\"%s\"\n",fmt);

    if (Qinv)
        ftp = ft_open(getparam("in"), mode,
                  getiparam("ycol"), getiparam("xcol"));
    else
        ftp = ft_open(getparam("in"), mode,
                  getiparam("xcol"), getiparam("ycol"));


    for (i=0; i<n; i++) {
        if (mode == FUNTAB_LINEAR)
            y = ft_linear(ftp, x[i]);
        else if (mode == FUNTAB_SPLINE)
            y = ft_spline(ftp, x[i]);
        else
            error("Bad mode");        
        printf(fmt,x[i], y);
    }

    ft_close(ftp);                  
}

#endif
