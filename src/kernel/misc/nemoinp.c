/*   nemoinp (archaic), nemoinpd, nemoinpf, nemoinpi, nemoinpl, nemoinpb
 *
 * NEMOINP:	derivatives of herinp(), more c-user friendly
 *		interface of the same
 *  returns a negative number (the herinp() error code) when a 
 *  parsing error occurred, and a positive (or 0) number of the
 *  the number of elements returned in 'a' when no error detected.
 *
 *  18-May-88   Peter Teuben
 *   6-Dec-88   added 'elen' parameter because of new herinp PJT
 *   6-Mar-89   split into nemoinpX functions
 *   5-jul-89   can't make up my mind To have an _ or not to have.. PJT
 *  14-oct-90   allow option with burststring and using getXparam.. PJT
 *  25-feb-92	happy gcc2.0 - kludged up NOHERINP		    PJT
 *   7-jan-93   bit more help for ignorant user in the TOOLBOX version PJT
 *   4-mar-94   V1.3 ansi - added nemoinpf for float's              pjt
 *   2-jul-95   V1.4 added tab= to have optional table output file  pjt
 *  18-sep-96   V1.5 changed default of newline from 'f' to 't'     pjt
 *  16-jun-97   V1.6 use integer  routines if %d is used	    pjt
 *  13-mar-00   V1.7 added seed= parameter			    pjt
 *                   (but it doesn't appear to do anything useful)
 */

#include <stdinc.h>
#include <getparam.h>

#if defined(NOHERINP)

int nemoinpi(expr,a,na)         /* parse integers */
char *expr;
int  *a;
int   na;
{
    int nret, atoi();
    string *vals, *burststring();

    vals = burststring(expr,", ");
    for (nret=0; vals[nret] != NULL; nret++) {
        if (nret>=na) {
            return(-23);
	}
        a[nret] = atoi(vals[nret]);
    }
    return(nret);    
}

int nemoinpl(expr,a,na)         /* parse longs */
char *expr;
long *a;
int   na;
{
    int nret, atoi();
    string *vals, *burststring();

    vals = burststring(expr,", ");
    for (nret=0; vals[nret] != NULL; nret++) {
        if (nret>=na) {
            return(-23);
        }
        a[nret] = atoi(vals[nret]);
    }
    return(nret);
}


int nemoinpb(expr,a,na)         /* parse booleans */
char *expr;
bool *a;
int   na;
{
    int nret;    
    string *vals, *burststring();

    vals = burststring(expr,", ");
    for (nret=0; vals[nret] != NULL; nret++) {
        if (nret>=na) {
            return(-23);
	}
        a[nret] = (*vals[nret]=='t' || *vals[nret]=='T') ? TRUE : FALSE;
    }
    return(nret);    
}
    

int nemoinpd(expr,a,na)         /* parse doubles */
char   *expr;
double *a;
int     na;
{
    int nret;
    string *vals, *burststring();
    double atof();

    vals = burststring(expr,", ");
    for (nret=0; vals[nret] != NULL; nret++) {
        if (nret>=na) {
            return(-23);
        }
        a[nret] = atof(vals[nret]);
    }
    return(nret);
}
    

int nemoinpf(expr,a,na)         /* parse floats */
char   *expr;
float  *a;
int     na;
{
    int nret;
    string *vals, *burststring();
    double atof();

    vals = burststring(expr,", ");
    for (nret=0; vals[nret] != NULL; nret++) {
        if (nret>=na) {
            return(-23);
        }
        a[nret] = atof(vals[nret]);
    }
    return(nret);
}

#else
extern void herinp();

int nemoinpi(expr,a,na)         /* parse integers */
char *expr;
int  *a;
int   na;
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(int);
    type = 'I';     /* integer */
    herinp (expr, &elen, &type, &tlen, a, &na, &nret, &ierr);
    if (ierr < 0)
        return(ierr);
    else
        return(nret);
}

int nemoinpl(expr,a,na)         /* parse longs */
char *expr;
long *a;
int   na;
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(long);
    type = 'I'; /* integer */
    herinp (expr, &elen, &type, &tlen, a, &na, &nret, &ierr);
    if (ierr < 0)
        return(ierr);
    else
        return(nret);
}


int nemoinpb(expr,a,na)         /* parse booleans */
char *expr;
bool *a;
int   na;
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(bool);
    type = 'L';     /* logical */
    herinp (expr, &elen, &type, &tlen, a, &na, &nret, &ierr);
    if (ierr < 0)
        return(ierr);
    else
        return(nret);
}
    

int nemoinpd(expr,a,na)         /* parse doubles */
char   *expr;
double *a;
int     na;
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(double);
    type = 'F';
    herinp (expr, &elen, &type, &tlen, a, &na, &nret, &ierr);
    if (ierr < 0)
        return(ierr);
    else
        return(nret);
}

int nemoinpf(expr,a,na)         /* parse floats */
char   *expr;
float  *a;
int     na;
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(float);
    type = 'F';
    herinp (expr, &elen, &type, &tlen, a, &na, &nret, &ierr);
    if (ierr < 0)
        return(ierr);
    else
        return(nret);
}
#endif

#if defined(TOOLBOX)

string defv[] = {
    "expression=\n	Expression to parse [help]",
    "separ=\n		Separator between numbers on output",
    "format=%g\n	Format used to print numbers",
    "newline=t\n	Use newline between numbers?",
    "nmax=32768\n	Size of output buffer",
    "tab=\n             Optional table output",
    "seed=0\n		Seed for xrandom",
    "VERSION=1.7\n	13-mar-00 PJT",
    NULL,
};


string usage = "expression parser and evaluator; also does lists";


nemo_main()
{
    char   fmt1[20], fmt2[20], *cp;
    double *x;
    int    *ix;
    int    i,nret, nx, seed;
    bool   Qnl, Qint;
    stream outstr;

    if (hasvalue("tab"))
        outstr = stropen(getparam("tab"),"w");
    else
        outstr = stdout;                /* is that really ok? */
    cp = getparam("format");            /* get format string to print with */
    if (strchr(cp,'%')==NULL)
      warning("%s badly formed printf-conversion specification ???",cp);
    if (strchr(cp,'d')) {
      dprintf(1,"Using integer math\n");
      Qint = TRUE;
    } else {
      dprintf(1,"Using floating point math\n");
      Qint = FALSE;
    }
    seed = set_xrandom(getiparam("seed"));
    dprintf(1,"set_xrandom: seed=%d\n",seed);

    strcpy (fmt1,cp);                   /* store it in 'fmt' */
    strcpy (fmt2,cp);                   /* and here */
    cp = getparam("separ");             
    if (hasvalue("separ"))              /* separator between numbers ? */
        strcat(fmt2,getparam("separ"));
    else
        strcat (fmt2," ");              /* else use blank as separator */

    Qnl = getbparam("newline");         /* use newlines also ? */

    if (!hasvalue("expression")) {      /* extra help if nothing given */
    	/* this somewhat unusual exit is because aliens often use this program */
	/* what really should have been done is make "expression=???" default */
        dprintf(0,"Usage: %s <expression>\n",getargv0());
        dprintf(0,"\n<expr> can be of form:  start[:end][:incr][,start::repeat]...\n");
        dprintf(0,"Also try keyword 'help= or help=h'\n");
        stop(0);
    }
    nx = getiparam("nmax");

    if (Qint) {
        ix = (int *) allocate(nx * sizeof(int));
   
        nret = nemoinpi(getparam("expression"),ix,nx);
        if (nret == -23)
            error("Too many items in list, use a bigger nmax=%d",nx);
        else if (nret < 0)
            error("nemoinp (%s) parsing error (%d)",cp,nret);

        for (i=0; i< nret-1; i++) {
            fprintf(outstr,fmt2,ix[i]);
            if (Qnl)
                fprintf(outstr,"\n");
        }
        fprintf (outstr,fmt1,ix[nret-1]);
        if (Qnl)
            fprintf (outstr,"\n");
    } else {
        x = (double *) allocate(nx * sizeof(double));
   
        nret = nemoinpd(getparam("expression"),x,nx);
        if (nret == -23)
            error("Too many items in list, use a bigger nmax=%d",nx);
        else if (nret < 0)
            error("nemoinp (%s) parsing error (%d)",cp,nret);

        for (i=0; i< nret-1; i++) {
            fprintf(outstr,fmt2,x[i]);
            if (Qnl)
                fprintf(outstr,"\n");
        }
        fprintf (outstr,fmt1,x[nret-1]);
        if (Qnl)
            fprintf (outstr,"\n");
    }
}
#endif
