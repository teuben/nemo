/*   nemoinp (archaic), nemoinpd, nemoinpf, nemoinpi, nemoinpl, nemoinpb
 *                      nemorinpd, nemorinpf, nemorinpi, nemorinpl, nemorinpb
 *                      natof, natoi
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
 *  31-may-01   V1.8 added natof, natoi                             pjt
 *   7-jun-01       a:  added casting to shut up compiler           pjt
 *   8-sep-01       b:  init_xrandom
 *   4-mar-03   V1.9: added nemoinpx for sexa decoding into degrees pjt
 *  28-jun-03       a: fixed prototype for darwin :-)               pjt
 *  24-nov-03       b: fixed another prototype for gcc3             pjt
 *  28-jan-04       c: recognize nan or NaN and return same         pjt
 *   6-jan-05   V1.10:  added the repeat nemorinpX routines         pjt
 *   5-may-05       a: add newline                                  pjt
 *   5-feb-06       b: also recognize -nan or -NaN ....             pjt
 *
 * TODO:  how to handle largest and smallest number (old style MIN/MAXLOG 38)
 */

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>

extern string *burststring(string,string);
extern void freestrings(string *);

#if defined(NOHERINP)

int nemoinpi(
	     char *expr,
	     int  *a,
	     int   na)
{
  int nret;
  string *vals;

    vals = burststring(expr,", ");
    for (nret=0; vals[nret] != NULL; nret++) {
        if (nret>=na) {
            return(-23);
	}
        a[nret] = atoi(vals[nret]);
    }
    return(nret);    
}

int nemoinpl(
	     char *expr,
	     long *a,
	     int   na)
{
  int nret;
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


int nemoinpb(
	     char *expr,
	     bool *a,
	     int   na)
{
  int nret;    
  string *vals;

    vals = burststring(expr,", ");
    for (nret=0; vals[nret] != NULL; nret++) {
        if (nret>=na) {
            return(-23);
	}
        a[nret] = (*vals[nret]=='t' || *vals[nret]=='T') ? TRUE : FALSE;
    }
    return(nret);    
}
    

int nemoinpd(
	     char   *expr,
	     double *a,
	     int     na)
{
  int nret;
  string *vals;

    vals = burststring(expr,", ");
    for (nret=0; vals[nret] != NULL; nret++) {
        if (nret>=na) {
            return(-23);
        }
        a[nret] = atof(vals[nret]);
    }
    return(nret);
}
    

int nemoinpf(
	     char   *expr,
	     float  *a,
	     int     na)
{
  int nret;
  string *vals;

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

extern void herinp(char *expr, int *nchr, char *type, int *length,
		   char *outv, int *nout, int *nret, int *ierd);


int nemoinpi(
	     char *expr,
	     int  *a,
	     int   na)
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(int);
    type = 'I';     /* integer */
    herinp (expr, &elen, &type, &tlen, (char *)a, &na, &nret, &ierr);
    if (ierr < 0)
        return(ierr);
    else
        return(nret);
}

int nemoinpl(
	     char *expr,
	     long *a,
	     int   na)
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(long);
    type = 'I'; /* integer */
    herinp (expr, &elen, &type, &tlen, (char *)a, &na, &nret, &ierr);
    if (ierr < 0)
        return(ierr);
    else
        return(nret);
}


int nemoinpb(
	     char *expr,
	     bool *a,
	     int   na)
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(bool);
    type = 'L';     /* logical */
    herinp (expr, &elen, &type, &tlen, (char *)a, &na, &nret, &ierr);
    if (ierr < 0)
        return(ierr);
    else
        return(nret);
}
    

int nemoinpd(
	     char   *expr,
	     double *a,
	     int     na)
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(double);
    type = 'F';
    herinp (expr, &elen, &type, &tlen, (char *)a, &na, &nret, &ierr);
    if (ierr < 0)
        return(ierr);
    else
        return(nret);
}

int nemoinpf(
	     char   *expr,
	     float  *a,
	     int     na)
{
    int ierr, nret, elen, tlen;
    char type;

    elen = strlen(expr);
    tlen = sizeof(float);
    type = 'F';
    herinp (expr, &elen, &type, &tlen, (char *)a, &na, &nret, &ierr);
    if (ierr < 0)
        return ierr;
    else
        return nret;
}
#endif


int nemoinpx(
	     char *expr,
	     double *a,
	     int     na)
{
  int nret, ncomp;
  string *vals, *comp;

  vals = burststring(expr,",");
  for (nret=0; vals[nret] != NULL; nret++) {
    if (nret>=na)
      return -23;
    comp = burststring(vals[nret],":");
    ncomp = xstrlen(comp,sizeof(string))-1;
    if (ncomp < 1 || ncomp > 3)
      return -13;
    a[nret] = atof(comp[0]);
    if (ncomp==1) continue;
    a[nret] += atof(comp[1])/60.0;
    if (ncomp==2) continue;
    a[nret] += atof(comp[2])/3600.0;
    freestrings(comp);
  }
  freestrings(vals);
  return nret;
}


int nemorinpi(string var, int *xvar,  int n, int defvar, bool repeated)
{
  int i, nret;
  if (var == 0 || *var == 0) {
    for (i=0; i<n; i++) 
      xvar[i] = defvar;
    return 0;
  } else {
    nret = nemoinpi(var,xvar,n);
    if (nret <= 0) error("nemorinpi: parsing error %d in %s",nret,var);
    for (i=nret; i<n; i++)
      xvar[i] = repeated ? xvar[i-1] : defvar;
    return nret;
  }
  return -1;
}


int nemorinpd(string var,  double *xvar, int n,double defvar, bool repeated)
{
  int i, nret;
  if (var == 0 || *var == 0) {
    for (i=0; i<n; i++) 
      xvar[i] = defvar;
    return 0;
  } else {
    nret = nemoinpd(var,xvar,n);
    if (nret <= 0) error("nemorinpd: parsing error %d in %s",nret,var);
    for (i=nret; i<n; i++)
      xvar[i] = repeated ? xvar[i-1] : defvar;
    return nret;
  }
  return -1;
}

int nemorinpf(string var, float *xvar,  int n, float defvar, bool repeated)
{
  int i, nret;
  if (var == 0 || *var == 0) {
    for (i=0; i<n; i++) 
      xvar[i] = defvar;
    return 0;
  } else {
    nret = nemoinpf(var,xvar,n);
    if (nret <= 0) error("nemorinpf: parsing error %d in %s",nret,var);
    for (i=nret; i<n; i++)
      xvar[i] = repeated ? xvar[i-1] : defvar;
    return nret;
  }
  return -1;
}

int nemorinpl(string var, long *xvar,  int n, long defvar, bool repeated)
{
  int i, nret;
  if (var == 0 || *var == 0) {
    for (i=0; i<n; i++) 
      xvar[i] = defvar;
    return 0;
  } else {
    nret = nemoinpl(var,xvar,n);
    if (nret <= 0) error("nemorinpl: parsing error %d in %s",nret,var);
    for (i=nret; i<n; i++)
      xvar[i] = repeated ? xvar[i-1] : defvar;
    return nret;
  }
  return -1;
}

int nemorinpb(string var, bool *xvar,  int n, bool defvar, bool repeated)
{
  int i, nret;
  if (var == 0 || *var == 0) {
    for (i=0; i<n; i++) 
      xvar[i] = defvar;
    return 0;
  } else {
    nret = nemoinpb(var,xvar,n);
    if (nret <= 0) error("nemorinpb: parsing error %d in %s",nret,var);
    for (i=nret; i<n; i++)
      xvar[i] = repeated ? xvar[i-1] : defvar;
    return nret;
  }
  return -1;
}





double natof(char *expr)
{
  double x;
  int n;
  if (streq(expr, "nan") || streq(expr, "NaN") || streq(expr,"-nan") || streq(expr,"-NaN"))
    return atof("nan");

  n = nemoinpd(expr,&x,1);
  return x;
}

int natoi(char *expr)
{
  int x, n;
  n = nemoinpi(expr,&x,1);
  return x;
}


#if defined(TOOLBOX)

string defv[] = {
    "expression=\n	Expression to parse [help]",
    "separ=\n		Separator between numbers on output",
    "format=%g\n	Format used to print numbers",
    "newline=t\n	Use newline between numbers?",
    "nmax=32768\n	Size of output buffer",
    "tab=\n             Optional table output",
    "seed=0\n		Seed for xrandom",
    "atof=\n            test (n)atof single value expression",
    "dms=f\n            Use D:M:S.SS parsing instead of regular",
    "VERSION=1.10c\n	23-feb-2019 PJT",
    NULL,
};


string usage = "expression parser and evaluator; also does lists";

string cvsid="$Id$";

void nemo_main(void)
{
    char   fmt1[20], fmt2[20], *cp;
    double dms[32];
    double *x;
    int    *ix;
    int    i,nret, nx, seed;
    bool   Qnl, Qint,Qdms;
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
    seed = init_xrandom(getparam("seed"));
    dprintf(1,"init_xrandom: seed=%d\n",seed);

    strcpy (fmt1,cp);                   /* store it in 'fmt' */
    strcpy (fmt2,cp);                   /* and here */
    
    Qdms = getbparam("dms");
    if (Qdms) {
      strcat(fmt1,"\n");
      nret = nemoinpx(getparam("expression"),dms,32);
      if (nret < 0) error("Parsing dms expression");
      for (i=0; i<nret; i++)
	printf(fmt1,dms[i]);
      return;  /* for now */
    }

    cp = getparam("separ");             
    if (hasvalue("separ"))              /* separator between numbers ? */
        strcat(fmt2,getparam("separ"));
    else
        strcat (fmt2," ");              /* else use blank as separator */

    Qnl = getbparam("newline");         /* use newlines also ? */

    if (hasvalue("atof")) {
      printf(fmt1,natof(getparam("atof")));
      printf("\n");
      return;
    }

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
    if (!Qnl) fprintf(outstr,"\n");  /* always make sure to add 1 newline */

}
#endif
