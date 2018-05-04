/*
 * XRANDOM: better (portable) floating-point random number routine.
 *
 *	-DNUMREC: a portable one from NumRec can be used (currently ran3)
 *      -DRAND48: the standard 'rand48' pseudo-random are used 
 *                (lrand48/srand48 instead of random/srandom)
 *
 *  xx-xxx-86   created                                     JEB
 *  xx-xxx-90   introduced the NUMREC portable option       PJT
 *  15-apr-91   set_xrandom() now returns it's seed!        PJT
 *   7-jan-92	fixed bug which made seed's in set_xrandom() useless	PJT
 *  18-may-92   documentation 
 *  29-jul-92   convex declaration fix - no more math       PJT
 *  18-feb-94   changed default to NUMREC (Solaris)  - ansi PJT
 *  24-mar-94   made -DRAND48 a full option		    PJT
 *   7-sep-95   added security to keep xrandom() in bounds  pjt (for old ran3!!)
 *  21-nov-96   added set_xrandom(-1) for centisecond control	    pjt
 *		TESTBED is now a TOOLBOX, default is n=0
 *   8-dec-96   rearranged order in TOOLBOX, also report output for unifom
 *  20-jan-99   add tab= to optionally output all random numbers  PJT
 *  24-feb-01   added comments to grandrom() and return both #'s  PJT
 *   7-apr-01   fixed grandom() bug, introduced 24-feb
 *   8-sep-01   (V2.0) added GSL 
 *  24-nov-03   V2.0b  some prototypes for -Wall added            PJT
 *  13-may-05   
 *  23-oct-07   added random0() and srandom0(), from Koda, based on ran1 from NumRec     PJT
 *  19-jun-08   enable init_xrandom(0) (previously caused Segmentation fault)  WD
 *  27-Sep-10   MINGW32/WINDOWS support  JCL
 *   9-oct-12   test the miriad method of drawing a gaussian                  PJT
 *  24-jan-18   faster version of grandom() with better caching               PJT
 *              1e7 grandom:   V2.2 -> 1.40"    V2.3 -> 0.85
 *
 *  See also: getrandom(2LINUX)
 */

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>

#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#ifndef __MINGW32__
#include <sys/times.h>
#endif

extern string *burststring(string,string);

#if defined(mc68k)
# define RAND48
#endif


#if defined(RAND48)
 extern double drand48(void);
 extern void   srand48(long);
#endif


#if defined(NUMREC)
# define portable_ran  ran3
real portable_ran(int *);     /* ieck, this is a long */
#endif

#if defined(RANDOM0)
extern void srandom0(long seed0);
extern double random0(void);
#endif

local int idum;            /* local variable to store used seed */

#ifdef HAVE_GSL
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
 static gsl_rng *my_r = NULL;
 static const gsl_rng_type *my_T;

 static string env_type = "GSL_RNG_TYPE";    /* environment variables    */
 static string env_seed = "GSL_RNG_SEED";    /* used by GSL_RNG routines */

#endif


int init_xrandom(string init)
{
#ifdef HAVE_GSL
    char my_gsl_type[64], my_gsl_seed[64], *cp;     /* need two strings, bug in putenv? */
    string *is;
    int nis, iseed;

    is = burststring(init,", ");                      /* parse init as "[seed[,name]]"  */
    nis = xstrlen(is,sizeof(string))-1;
    if (nis > 0) {                                          /* seed is first, but optional */
        iseed = natoi(is[0]);
	if (iseed > 0 || streq(is[0],"+0")) {
	  sprintf(my_gsl_seed,"%s=%s",env_seed,is[0]);
	} else {
	  iseed = set_xrandom(iseed);
	  sprintf(my_gsl_seed,"%s=%u",env_seed,iseed);
	}
	putenv(my_gsl_seed);
	dprintf(1,"putenv: %s\n",my_gsl_seed);
        if (nis > 1) {                                      /* name is second, also optional */
            sprintf(my_gsl_type,"%s=%s",env_type,is[1]);
            putenv(my_gsl_type);
	    dprintf(1,"putenv: %s\n",my_gsl_type);
        }
    }

    gsl_rng_env_setup();                          /* initialize the rng (name/seed) setup */
    my_T = gsl_rng_default;
    my_r = gsl_rng_alloc(my_T);


    dprintf(1,"GSL generator type: %s\n",gsl_rng_name(my_r));
    dprintf(1,"GSL seed = %u\n",gsl_rng_default_seed);
    dprintf(1,"GSL first value = %u\n",gsl_rng_get(my_r));

    return (int) gsl_rng_default_seed;
#else
    return set_xrandom(init? natoi(init) : 0);     /*  18/06/2008: allow for init=0 WD */
#endif
}


int set_xrandom(int dum)
{
    int retval;
#ifndef __MINGW32__
    struct tms buffer;
#endif
    if (dum <= 0) {
	if (dum == -1)
#ifndef __MINGW32__
            retval = idum = (int) times(&buffer);   /* clock cycles */
#else
	    ;
#endif
        else if (dum == -2)
            retval = idum = (int) getpid();         /* process id */
        else            /* normally if dum==0 */
	    retval = idum = (int) time(0);          /* seconds since 1970 */
    } else
    	retval = idum = dum;	           /* use supplied seed in argument */

#if !defined(HAVE_GSL)
#if defined(NUMREC)
    dprintf(2,"set_xrandom(NUMREC portable) seed=%d\n",idum);
    if (idum > 0) idum = -idum; /* ran needs a negative idum as seed */
    portable_ran(&idum);        /* use portable seed: will only work once */
    idum = 0;                   /* reset to avoid reentry ... */
#elif defined(RAND48)
    dprintf(2,"set_xrandom(UNIX srand48) seed=%d\n",retval);
    srand48(retval);
#elif defined(RANDOM0)
    dprintf(2,"set_xrandom(random0) seed=%d\n",retval);
    srandom0(retval);
#else    
    dprintf(2,"set_xrandom(UNIX srandom) seed=%d\n",retval);
    srandom(retval);
#endif
#endif
    return retval;

}

double xrandom(double xl, double xh)
{
    double retval;

    for(;;) {
#if defined(HAVE_GSL)
        if (my_r == NULL) error("GSL init_xrandom was never called");
        retval = gsl_rng_uniform(my_r);
#elif defined(NUMREC)
        retval = ((double) portable_ran(&idum));
#elif defined(RAND48)
        retval = drand48();
#elif defined(RANDOM0)
        retval = random0();
#else    
        retval = (double) random() / 2147483647.0;
#endif
        if (retval<0.0 || retval>1.0) {
            warning("xrandom: spinning again, out of bounds [%g]",retval);
            continue;
        } 
        break;
    }
    return xl + retval*(xh-xl);
}

/*
 * GRANDOM: normally distributed random number (polar method)
 *	    referred to as the Box-Mueller method, although
 *	    Gauss also attributes it to Laplace.
 * See also: Knuth, vol. 2, p. 104.
 */

#if 0
	/* better check */
do {
     v1=2.0*ran2()-1.0;
     v2=2.0*ran2()-1.0;
     s=v1*v1+v2*v2;
} while(s>1.0);
ar=log(s)
r1=v1*sqrt(-(ar+ar)/s)*sig-mu;
r2=v2*sqrt(-(ar+ar)/s)*sig-mu;

#endif

double grandom(double mean, double sdev)
{
    static double v1, v2, s;
    static int gcount = 0;

    if (gcount) {
        gcount = 0;
        return mean + sdev * v2 * s;
    }

    do {				/* loop until */
	v1 = xrandom(-1.0, 1.0);	/* two points within */
	v2 = xrandom(-1.0, 1.0);	/* the unit square */
	s = v1*v1 + v2*v2;		/* have radius such that */
    } while (s >= 1.0);			/* point fall inside unit circle */
    s = sqrt(-2.0 * log(s) / s);
    
    gcount = 1;
    return mean + sdev * v1 * s;

}




#if defined(TOOLBOX)

#include <moment.h>

string defv[] = {
#ifdef HAVE_GSL
    "seed=0,mt19937\n Seed [0=seconds_1970, -1=centisec_boot -2=pid], and optional GSL name",
#else
    "seed=0\n       Seed [0=seconds_1970, -1=centisec_boot -2=pid]",
#endif
    "n=0\n          Number of random numbers to draw",
    "m=1\n          Number of times to repeat experiment (enforces tab=f)",
    "gauss=f\n      gaussian or uniform noise?",
    "report=f\n     Report mean/dispersian/skewness/kurtosis?",
    "tab=t\n        Tabulate all random numbers?",
    "offset=0.0\n   Offset of distribution from 0",
#ifdef HAVE_GSL
    "gsl=\n         If given, GSL distribution name",
    "pars=\n        Parameters for GSL distribution",
#endif
    "VERSION=2.3\n  24-jan-2018 PJT",
    NULL,
};

string usage="Return seed for random numbers and optionally random numbers";


/* argsssssssssss , this is for seed=1 */

string defaults[] = {
    "srandom:        0.968071  0.0667306 0.478281  0.909534",
    "srand48:        0.0416303 0.454492  0.834817  0.335986",
    "numrec ran3:    0.715119  0.0330211 0.874394  0.534194",
    NULL,
};

#define MAXPARS 5

static bool check(string, string, string, int, int);

void nemo_main()
{
  int i, j, k, seed, n, m, npars;
  double sum[5], mean, sigma, skew, kurt, x, y, s;
  double p[MAXPARS];
  unsigned int up[MAXPARS];
  bool   Qgauss = getbparam("gauss");
  bool   Qreport = getbparam("report");
  bool   Qtab = getbparam("tab");
  bool   Qbench;
  string *sp, ran_name;
  Moment mom;
  real   offset = getdparam("offset");

  
  n = getiparam("n");
  m = getiparam("m");
  seed = init_xrandom(getparam("seed"));
  Qbench = (n==4 && seed==1);
  
  if (m>1) {
    Qtab= Qreport = FALSE;   /* override */
    ini_moment(&mom,4,0);
  } else {
    ini_moment(&mom,0,0);
  }

  printf("Seed used = %d\n",seed);
  
  for (k=0; k<m; k++) {
    sum[0] = sum[1] = sum[2] = sum[3] = sum[4] = 0.0;
    if (Qgauss) {
      for (j=0; j<n;j++) {
	for (i=0, s=1.0, y = grandom(offset,1.0); i<5; i++, s *= y)
	  sum[i] += s;
	if (Qtab) printf("%g\n",y);
      }
    } else {
      for (j=0; j<n;j++) {
	for (i=0, s=1.0, y = xrandom(offset,1.0); i<5; i++, s *= y)
	  sum[i] += s;
	if (Qtab) printf("%g\n",y);	
      }
    }

    if (m>1) {
      accum_moment(&mom, sum[1]/sum[0], 1.0);
    }
    
    if (Qreport && n > 1) {
      mean  = sum[1]/sum[0];
      sigma = sqrt(sum[2]/sum[0]-mean*mean);
      skew  = ( (sum[3]-3*sum[2]*mean)/sum[0] + 2*mean*mean*mean) /
	(sigma*sigma*sigma);
      kurt  = ( (sum[4]-4*sum[3]*mean+6*sum[2]*mean*mean)/sum[0]
		- 3*mean*mean*mean*mean) /
	(sigma*sigma*sigma*sigma)  - 3.0;
      
      printf ("Mean and dispersion  : %f %f\n",mean,sigma); 
      printf ("Skewness and kurtosis: %f %f\n",skew,kurt);
    }
  } /* k < m */

  if (m>1) {
    printf ("N                    : %d\n",n_moment(&mom));
    printf ("Mean and dispersion  : %f %f\n",mean_moment(&mom),sigma_moment(&mom));
    printf ("Skewness and kurtosis: %f %f\n",skewness_moment(&mom),kurtosis_moment(&mom));
  }



  if (!Qgauss && Qbench) {
    dprintf(0,"Known n=4 seed=1 cases are:\n");
    for (sp=defaults; *sp; sp++) dprintf(0,"%s\n",*sp);
  }


#ifdef HAVE_GSL
  if (hasvalue("gsl")) {
      ran_name = getparam("gsl");
      npars = nemoinpd(getparam("pars"),p,MAXPARS);
      if (npars < 0) error("Error parsing pars=%s",getparam("pars"));
      for (i=0; i<npars; i++) {
	up[i] = (unsigned int) p[i];
	dprintf(1,"par(%d) = %g\n",i+1,p[i]);
      }

      if (       check(ran_name,"gaussian",              "sigma",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_gaussian(my_r, p[0]));
      } else if (check(ran_name,"gaussian_ratio_method", "sigma",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_gaussian_ratio_method(my_r, p[0]));
      } else if (check(ran_name,"gaussian_tail",         "a sigma",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_gaussian_tail(my_r, p[0],p[1]));
      } else if (check(ran_name,"exponential",           "mu",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_exponential(my_r, p[0]));
      } else if (check(ran_name,"laplace",               "a",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_laplace(my_r, p[0]));
      } else if (check(ran_name,"exppow",                "a b",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_exppow(my_r, p[0],p[1]));
      } else if (check(ran_name,"cauchy",                "a",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_cauchy(my_r, p[0]));
      } else if (check(ran_name,"rayleigh",              "sigma",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_rayleigh(my_r, p[0]));
      } else if (check(ran_name,"rayleigh_tail",         "a sigma",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_rayleigh_tail(my_r, p[0],p[1]));
      } else if (check(ran_name,"landau",                "-none-",0,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_landau(my_r));
      } else if (check(ran_name,"levy",                  "c alpha",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_levy(my_r, p[0],p[1]));
      } else if (check(ran_name,"levy_skey",             "c alpha beta",3,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_levy_skew(my_r, p[0],p[1],p[2]));
      } else if (check(ran_name,"gamma",                 "a b",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_gamma(my_r, p[0],p[1]));
      } else if (check(ran_name,"flat",                  "a b",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_flat(my_r, p[0],p[1]));
      } else if (check(ran_name,"lognormal",             "zeta sigma",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_lognormal(my_r, p[0],p[1]));
      } else if (check(ran_name,"chisq",                 "nu",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_chisq(my_r, p[0]));
      } else if (check(ran_name,"fdist",                 "nu1 nu2",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_fdist(my_r, p[0],p[1]));
      } else if (check(ran_name,"tdist",                 "nu",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_tdist(my_r, p[0]));
      } else if (check(ran_name,"beta",                  "a b",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_beta(my_r, p[0],p[1]));
      } else if (check(ran_name,"logistic",              "a",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_logistic(my_r, p[0]));
      } else if (check(ran_name,"pareto",                "a b",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_pareto(my_r, p[0],p[1]));
      } else if (check(ran_name,"weibull",               "a b",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_weibull(my_r, p[0],p[1]));
      } else if (check(ran_name,"gumbel1",               "a b",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_gumbel1(my_r, p[0],p[1]));
      } else if (check(ran_name,"gumbel2",               "a b",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", gsl_ran_gumbel2(my_r, p[0],p[1]));

	/* the next ones are actually discrete (integer) functions  */

      } else if (check(ran_name,"poisson",               "mu",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", (double)gsl_ran_poisson(my_r, p[0]));
      } else if (check(ran_name,"bernoulli",             "p",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", (double) gsl_ran_bernoulli(my_r, p[0]));
      } else if (check(ran_name,"binomial",              "p n",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", (double) gsl_ran_binomial(my_r, p[0],up[1]));
      } else if (check(ran_name,"negative_binomial",     "p n",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", (double) gsl_ran_negative_binomial(my_r, p[0],p[1]));
      } else if (check(ran_name,"pascal",                "p k",2,npars)) {
	for (i=0; i<n; i++) printf("%g\n", (double) gsl_ran_pascal(my_r, p[0],up[1]));
      } else if (check(ran_name,"geometric",             "p",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", (double) gsl_ran_geometric(my_r, p[0]));
      } else if (check(ran_name,"hypergeometric",        "n1 n2 nt",3,npars)) {
	for (i=0; i<n; i++) printf("%g\n", (double) gsl_ran_hypergeometric(my_r, up[0],up[1],up[2]));
      } else if (check(ran_name,"logarithmic",           "sigma",1,npars)) {
	for (i=0; i<n; i++) printf("%g\n", (double) gsl_ran_logarithmic(my_r, p[0]));
      } else
	error("Bad name");
  }
#endif
}

static bool check(string a, string b, string msg, int k, int l)
{
  if (!streq(a,b)) return FALSE;
  if (k!=l) error("%s needs %d parameters (%s)",a,k,msg);
  return TRUE;
}

#endif
