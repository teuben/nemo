/*
 *  SNAPOPT: Ostriker-Peebles little-t calculator
 *
 *	28-oct-93	V1.0 created      		PJT
 *			with hardcoded bodytrans() routines	PJT
 *      12-nov-93       V1.1 free() to force alloc()    PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n		   Input file (snapshot)",
    "dr=0.1\n              Ringsize",
    "rmax=\n               Force rmax? Otherwise autoscale",
    "w=0\n                 Additive external self-potential",
    "times=all\n	   Times to select snapshot",
    "tab=\n		   Standard output or table file?",
    "VERSION=1.1a\n	   25-mar-94 PJT",
    NULL,
};

string usage="Special Ostriker-Peebles 't' calculator";

real btr_opt();

nemo_main()
{
    stream instr;
    real   tsnap, *x=NULL, *y, *w, rmax, dr, rmax_all;
    real *sum0, *sum1, *sum2, sum, mean, sigma, sumt, sumw, t, wext;
    bool Qpot;
    string headline=NULL, times;
    Body *btab = NULL, *bp;
    int i, j, nbody, bits, bins, xmom, ymom;

    instr = stropen(getparam("in"), "r");	/* open input file */

    times = getparam("times");
    wext = getdparam("w");
    dr = getdparam("dr");
    rmax_all = hasvalue("rmax") ? getdparam("rmax") : -1;
    
    get_history(instr);                 /* read history */

    for(;;) {                /* repeating until first or all times are read */
	get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
            break;                       /* done */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ( (bits & PhaseSpaceBit) == 0)
            continue;                   /* skip work, only diagnostics here */
        dprintf(1,"Processing snapshot @ 0x%x time=%g\n",btab,tsnap);
        Qpot = bits & PotentialBit;

        x = (real *) allocate(nbody*sizeof(real));
        y = (real *) allocate(nbody*sizeof(real));
        w = (real *) allocate(nbody*sizeof(real));
        rmax = 0.0;
        sumw = wext;      /* start with external self-potential DUBIOUS */
        for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
            x[i] = sqrt(sqr(Pos(bp)[0]) + sqr(Pos(bp)[1]));
            y[i] = btr_opt(bp,tsnap,i);
            w[i] = Mass(bp);
            rmax = MAX(rmax, x[i]);
            if (Qpot) sumw += Mass(bp)*Phi(bp);
        }
        sumw /= 2.0;			/* Potential was double counted before */
        dprintf(1,"RMax = %g, RMax_all\n",rmax,rmax_all);
        if (rmax_all > 0) rmax = rmax_all;
        bins = (int) floor(rmax/dr + 1.0);

        sum0 = (real *) allocate(bins*sizeof(real));
        sum1 = (real *) allocate(bins*sizeof(real));
        sum2 = (real *) allocate(bins*sizeof(real));

        for (j=0; j<bins; j++) 
            sum0[j] = sum1[j] = sum2[j] = 0.0;

        for (i=0; i<nbody; i++) {
            j = (int) floor(x[i]/dr);
            if (j<0 || j>=bins) continue;
            sum0[j] += w[i];
            sum1[j] += w[i]*y[i];
            sum2[j] += w[i]*y[i]*y[i];
        }

        sumt = 0.0;
        for (j=0; j<bins; j++) {
            sum   = sum0[j];					/* n\alpha */
            mean  = (sum0[j] > 0 ? sum1[j]/sum0[j] : 0.0);	/* v\alpha */
            sigma = (sum0[j] > 1 ? sum2[j]/sum0[j]-mean*mean : 0.0);
            if (sigma>0) sigma=sqrt(sigma);
            dprintf(1,"%g %g %g %g\n",(j+0.5)*dr, sum,mean,sigma);
            sumt += 0.5*sum*mean*mean;
        }
	t = (sumw==0 ? 0.0 : ABS(sumt/sumw));
        printf("%g %g %g %g\n",tsnap,sumt,sumw,t);
	/* free it all again */
        free(btab);  				btab=NULL;
	free(sum0); free(sum1); free(sum2); 	sum0=sum1=sum2=NULL;
	free(x); free(y); free(w);   		x=y=w=NULL;
    }
    strclose(instr);
}

#include <bodytrans.h>
real btr_opt(b,t,i)
Body *b;
real t;
int  i;
{
  real sn,vttot;
#if 0
  vttot=((vx*vx + vy*vy + vz*vz) -
		   sqr(x*vx + y*vy + z*vz) / (x*x + y*y + z*z));
#else  
  vttot=((vx*vx + vy*vy) -
		   sqr(x*vx + y*vy) / (x*x + y*y));
#endif
  sn = (x*vy-y*vx > 0 ? 1 : -1);

  return -sn*sqrt(vttot);

}



