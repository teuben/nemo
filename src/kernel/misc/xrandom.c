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
 *   7-sep-01   added GSL testing for ran3 :-)
 */

#include <stdinc.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>

#if defined(mc68k)
# define RAND48
#endif


#if defined(RAND48)
extern double drand48(void);
extern void   srand48(long);
#endif


#if defined(NUMREC)
# define portable_ran  ran3
real portable_ran(int *);
#endif

local int idum;            /* local variable to store used seed */

#ifdef HAVE_GSL
#include <gsl/gsl_rng.h>
#endif


int set_xrandom(int dum)
{
    int retval;
    struct tms buffer;

    if (dum <= 0) {
	if (dum == -1)
            retval = idum = (int) times(&buffer);   /* clock cycles */
        else if (dum == -2)
            retval = idum = (int) getpid();         /* process id */
        else            /* normally if dum==0 */
	    retval = idum = (int) time(0);          /* seconds since 1970 */
    } else
    	retval = idum = dum;	           /* use supplied seed in argument */

#if defined(NUMREC)
    dprintf(2,"set_xrandom(NUMREC portable) seed=%d\n",idum);
    if (idum > 0) idum = -idum; /* ran needs a negative idum as seed */
    portable_ran(&idum);        /* use portable seed: will only work once */
    idum = 0;                   /* reset to avoid reentry ... */
#elif defined(RAND48)
    dprintf(2,"set_xrandom(UNIX srand48) seed=%d\n",retval);
    srand48(retval);
#else    
    dprintf(2,"set_xrandom(UNIX srandom) seed=%d\n",retval);
    srandom(retval);
#endif    
    return retval;
}

double xrandom(double xl, double xh)
{
    double retval;

    for(;;) {
#if defined(NUMREC)
        retval = ((double) portable_ran(&idum));
#elif defined(RAND48)
        retval = drand48();
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
        return mean + sdev * v2 * sqrt(-2.0 * log(s) / s);
    }

    do {				/* loop until */
	v1 = xrandom(-1.0, 1.0);	/* two points within */
	v2 = xrandom(-1.0, 1.0);	/* the unit square */
	s = v1*v1 + v2*v2;		/* have radius such that */
    } while (s >= 1.0);			/* point fall inside unit circle */

    gcount = 1;
    return mean + sdev * v1 * sqrt(-2.0 * log(s) / s);

}



#if defined(TOOLBOX)

#include <getparam.h>

string defv[] = {
    "seed=0\n       Seed [0=seconds_1970, -1=centisec_boot -2=pid",
    "n=0\n          Number of random numbers to draw",
    "gauss=f\n      gaussian or uniform noise?",
    "report=f\n     Report mean/dispersian/skewness/kurtosis?",
    "tab=t\n        Tabulate all random numbers?",
    "VERSION=1.2\n  7-20-jan-99 PJT",
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

nemo_main()
{
    int i, seed, n;
    double sum[5], mean, sigma, skew, kurt, y, s;
    bool   Qgauss = getbparam("gauss");
    bool   Qreport = getbparam("report");
    bool   Qtab = getbparam("tab");
    bool   Qbench;
    string *sp;

    n = getiparam("n");
    seed = getiparam("seed");
    Qbench = (n==4 && seed==1);

    seed = set_xrandom(seed);
    printf("Seed used = %d\n",seed);

    sum[0] = sum[1] = sum[2] = sum[3] = sum[4] = 0.0;
    if (Qgauss) {
        while (n-- > 0) {
            for (i=0, s=1.0, y = grandom(0.0,1.0); i<5; i++, s *= y)
                sum[i] += s;
            if (Qtab) printf("%g\n",y);
        }
    } else {
        while (n-- > 0) {
            for (i=0, s=1.0, y = xrandom(0.0,1.0); i<5; i++, s *= y)
                sum[i] += s;
            if (Qtab) printf("%g\n",y);	
        }
    }

    n = getiparam("n");
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

        
    if (!Qgauss && Qbench) {
            dprintf(1,"Known n=4 seed=1 cases are:\n");
            for (sp=defaults; *sp; sp++) dprintf(1,"%s\n",*sp);
    }
}

#endif
