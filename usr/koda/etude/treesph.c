/****************************************************************************/
/* TREESPH.C: routines to compute hydro. force.                             */
/* Public routines: sphcalc(), stephknl(), checkhknl()                      */
/* Copyright (c) 2000 by Jin Koda, Tokyo, JAPAN.                            */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "treecode.h"

/*
 * Local routines and variables to perform neighbor search.
 */

local void sphdensity(bodyptr, int);     /* calculate SPH density of bodies  */
local void sphforce(bodyptr, int);       /* calculate SPH force on bodies    */
local real calcdtcfl(bodyptr, real);     /* calculate CFL time criterion     */
local real kernel(bodyptr, bodyptr, real);
local real gradkernel(bodyptr, bodyptr, real);
local void calckernel(bodyptr, bodyptr, real, real *, real *);
local void calckernel2(bodyptr, bodyptr, real, real *, real *);

/*
 * SPHCALC: common part of sph calculation.
 */

void sphcalc(bodyptr btab, int nbody)
{
    ngbrlen = nbody * nmax;                     /* length of ngbr list      */
    ngbrlist = (bodyptr *) allocate(ngbrlen * sizeof(bodyptr));
                                                /* allocate mem. for list   */
    makengbrlist(btab, nbody);                  /* construct neigbor list   */
    sphdensity(btab, nbody);                    /* compute densities        */
    sphforce(btab, nbody);                      /* compute hydro. forces    */
    free(ngbrlist);                             /* release memory           */
}

/*
 * SPHDENSITY: compute density of SPH particles.
 */

local void sphdensity(bodyptr btab, int nbody)
{
    double cpustart;
    bodyptr p, *q;
    real dr2, dv2, dotrv, adiv, arot, knl, gknl;
    vector dr, dv;
#if defined(THREEDIM)
    vector rotvr;
#else
    real rotvr;
#endif

    cpustart = cputime();
    for (p = btab; p < btab+nbody; p++) {      
	Rho(p) = 0.0;                             /* initialize density      */
	Divv(p) = 0.0;                            /* divergence of vel       */
#if defined(THREEDIM)
	CLRV(Rotv(p));                            /* rotation of vel         */
#else
	Rotv(p) = 0.0;                            /* rotation of vel         */
#endif	
    }
    for (p = btab; p < btab+nbody; p++) {         /* loop over particles     */
	for (q = Ngbrs(p); q < Ngbrf(p); q++) {   /* loop over ngbr list     */
	    if (Hknl(p) > Hknl(*q) || (Hknl(p) == Hknl(*q) && p > *q) ) {
                                                  /* if p larger than *q     */
		DOTPSUBV(dr2, dr, Pos(p), Pos(*q)); /* compute separation    */
		DOTPSUBV(dv2, dv,Velm(p),Velm(*q)); /* diff of velocity      */
		calckernel(p, *q, dr2, &knl, &gknl);/* and kernels           */
		Rho(p)  += Mass(*q) * knl;        /* increment density for p */
		Rho(*q) += Mass( p) * knl;        /* and for *q              */
		DOTVP(dotrv, dr, dv);             /* compute dot product     */
		dotrv *= gknl;                    /* weighted by grd knl     */
		Divv(p) -= Mass(*q) * dotrv;      /* increment div vel for p */
		Divv(*q) -= Mass(p) * dotrv;      /* and for *q              */
		CROSSVP(rotvr, dv, dr);           /* compute cross product   */
#if defined(THREEDIM)
		MULVS(rotvr, rotvr, gknl);        /* weighted by grd knl     */
		ADDMULVS(Rotv(p),rotvr,Mass(*q)); /* increment rot vel for p */
		ADDMULVS(Rotv(*q),rotvr,Mass(p)); /* and for *q              */
#else
		rotvr *= gknl;                    /* weighted by grd knl     */
		Rotv(p) += Mass(*q) * rotvr;      /* increment rot vel for p */
		Rotv(*q) += Mass(p) * rotvr;      /* and for *q              */
#endif
		
	    }
	}
    }
    for (p = btab; p < btab+nbody; p++) {
	Rho(p) += Mass(p) * kernel(p, p, 0.0);    /* density due to itself   */
	Divv(p) /= Rho(p);                        /* get divV and rotV of p  */
	adiv = ABS(Divv(p));                      /* get absolute of div     */
#if defined(THREEDIM)
	DIVVS(Rotv(p), Rotv(p), Rho(p));          /* by deviding by density  */
	ABSV(arot, Rotv(p));                      /* and rot of velocity     */
#else
	Rotv(p) /= Rho(p);                        /* by deviding by density  */
	arot = ABS(Rotv(p));                      /* and rot of velocity     */
#endif
	Vsprs(p) = adiv / (adiv + arot + 0.0001 * Csnd(p) / Hknl(p));
                                                  /* visc suppress factor    */
    }
    cpudens = cputime() - cpustart;
}

/*
 * SPHFORCE: compute hydrodynamical force of SPH particles.
 */

local void sphforce(bodyptr btab, int nbody)
{
    double cpustart;
    bodyptr p, *q;
    real pgp, pgq, cpq, hpq, rpq, mpq, vpq;
    real gknl, iknl, eps2, dr2i, phig;
    real dr2, dv2, dotrv, visc, pgrd, apq, dt0;
    vector dr, dv;

    cpustart = cputime();
    dtcfl = 1.0e9;
    eps2 = rpow(eps, 2.0);
    for (p = btab; p < btab+nbody; p++) {         /* loop over particles     */
	pgp = rsqr(Csnd(p)) / Rho(p);             /* comp grad of prs term   */
	for (q = Ngbrs(p); q < Ngbrf(p); q++) {   /* loop over ngbr list     */
	    if (Hknl(p) > Hknl(*q) || (Hknl(p) == Hknl(*q) && p > *q)) {
                                                  /* if p larger than q      */
		DOTPSUBV(dr2, dr, Pos(p), Pos(*q)); /* compute separation    */
		DOTPSUBV(dv2, dv,Velm(p),Velm(*q)); /* diff of velocity and  */
		DOTVP(dotrv, dr, dv);             /* dot product of p and *q */
		if (dotrv >= 0.0) {               /* if p and *q separating  */
		    mpq  = 0.0;
		    visc = 0.0;                   /* then no viscosity       */
		} else {                          /* else comp mean of       */
		    hpq = 0.5 * (Hknl(p)  + Hknl(*q)); /* smoothing radius   */
		    rpq = 0.5 * (Rho(p)   +  Rho(*q)); /* density            */
		    cpq = 0.5 * (Csnd(p)  + Csnd(*q)); /* sound velocity     */
		    vpq = 0.5 * (Vsprs(p) +Vsprs(*q)); /* visc spprss fact.  */
		    mpq = hpq * dotrv / (dr2 + 0.01 * hpq * hpq);
		    visc = (-alpha * cpq + beta * mpq) * mpq / rpq * vpq;
                                                  /* get viscousity          */
		}
		pgq = rsqr(Csnd(*q)) / Rho(*q);   /* grad of prs term of *q  */
		pgrd = pgp + pgq;                 /* sum of grad prs term    */
		if (selfgrav) {                   /* if calculate self-grav  */
                                                  /* correct grav and potent */
		    calckernel2(p, *q, dr2, &gknl, &iknl);
		    dr2i = 1.0 / (dr2 + eps2);
		    phig = (iknl - 1.0) * rsqrt(dr2i);
		    apq = - (pgrd + visc) * gknl - phig * dr2i;
                                                  /* sum up prg and grv crct */
		    ADDMULVS(Acc(p), dr, Mass(*q) * apq);
		    ADDMULVS(Acc(*q), dr, -Mass(p)* apq);
		    Phi(p) -= Mass(*q) * phig;    /* add corcted acc and pot */
		    Phi(*q) -= Mass(p) * phig;    /* between p and *q        */
		} else {                          /* else calc accel of hydr.*/
		    apq = - (pgrd + visc) * gradkernel(p, *q, dr2);
		    ADDMULVS(Acc(p), dr, Mass(*q) * apq);
		    ADDMULVS(Acc(*q), dr, -Mass(p)* apq);
		}
		dt0 = calcdtcfl(p, mpq);          /* comp CFL timestep for p */
		dtcfl = MIN(dtcfl, dt0);
		dt0 = calcdtcfl(*q, mpq);         /* comp CFL timestep for *q*/
		dtcfl = MIN(dtcfl, dt0);
	    }
	}
    }
    cpusphf = cputime() - cpustart;
}

/*
 * CALCDTCFL: calculate the CFL time criterion.
 */

local real calcdtcfl(bodyptr p, real mpq)
{
    real dt;

    dt = Hknl(p) / (Csnd(p) + 0.6 * (alpha * Csnd(p) +beta * ABS(mpq)));
    return(dt);
}
    
/*
 * STEPHKNL: update SPH smoothing length.
 * reference to eq.(2.17) of Katz & Hernquist (ApJS, 70, 419, 1989).
 */

void stephknl(bodyptr btab, int nbody, real dt)
{
    bodyptr p;
    int numngbr;

    for (p=btab; p<btab+nbody; p++) {             /* loop over particles     */
	numngbr = (int) (Ngbrf(p)-Ngbrs(p));      /* comp number of ngbrs    */
	if (numngbr < 1)
	    error("advancehknl: number of ngbr less than 1.\n");
#if defined(THREEDIM)
	Hknl(p) = Hknl(p) * 0.5 *                 /* advance smooth length   */
	    (1.0 + rpow( (real) nnbr / (real) numngbr, ONETRD));
#else
	Hknl(p) = Hknl(p) * 0.5 *                 /* advance smooth length   */
	    (1.0 + rpow( (real) nnbr / (real) numngbr, 0.5));
#endif
    }
}

/*
 * CHECKHKNL: check and correct SPH smooting radius if not in accepted range.
 */

bodyptr *checkhknl(nodeptr p, nodeptr *nptr, nodeptr *mptr)
{
    nodeptr *np;
    int numngbr;

    numngbr = (int) (Ngbrf(p)-Ngbrs(p));        /* comp number of neighbors */
    if ( numngbr > nmax ) {                     /* if it is too large       */
	Hknl(p) = hcorrect((bodyptr) p, Ngbrs(p), Ngbrf(p), nmax);
                                                /* set to nmax closest dist.*/
	Ngbrf(p) = Ngbrs(p) + nmax;             /* end of ngbr list for p   */
    } else if ( numngbr < nmin ) {              /* else too small           */
	np = nptr;                              /* use residual of act list */
	*np++ = (nodeptr) root;                 /* and descend from root    */
	while ((Ngbrf(p)-Ngbrs(p)) < nmin) {    /* loop to find enuf ngbrs  */
	    Hknl(p) *= 1.2;                     /* set larger searching rad.*/
	    scanngbr(nptr, np, mptr, mptr, p, Ngbrs(p));
                                                /* and scan ngbrs in it     */
	}
	Hknl(p) = hcorrect((bodyptr) p, Ngbrs(p), Ngbrf(p), nmin);
                                                /* set to nmin closest dist.*/
	Ngbrf(p) = Ngbrs(p) + nmin;             /* end of ngbr list for p   */
    }
    return(Ngbrf(p));                           /* return final pointer     */
}

/*
 * HCORRECT: correct SPH smoothing radius if it sticks out from accepted range.
 * select mth smallest distance and return it. neighbor list is also corrected.
 */

#define SWAP(a, b) {int temp; temp = (a); (a) = (b); (b) = temp;} 
#define DR2(x) dr2[indx[(x)]]

real hcorrect(bodyptr p, bodyptr *start, bodyptr *finish, int mth)
{
    bodyptr *q, *nlst;
    int *indx, numngbr, mth1, l, r, m, i, j;
    real *dr2, newh;
    vector dr;

    numngbr = (int) ( finish - start );           /* compute number of ngbrs */
    indx = (int *) allocate(numngbr * sizeof(int));
    nlst = (bodyptr *) allocate(numngbr * sizeof(bodyptr));
    dr2 = (real *) allocate(numngbr * sizeof(real));
    for (i = 0; i < numngbr; i++) {               /* loop over all ngbrs     */
	q = start + i;                            /* set "i"th ngbr pointer, */
	indx[i] = i;                              /* index of ngbrs          */
	nlst[i] = *q;                             /* and their body number   */
	DOTPSUBV(DR2(i), dr, Pos(p), Pos(*q));    /* comp distance to ngbrs  */
    }
    l = 0;                                        /* left and right end of   */
    r = numngbr - 1;                              /* first active partition  */
    mth1 = mth - 1;                               /* find mth nearest dist.  */
    while (r-l > 1) {                             /* start iteration         */
	m = (l + r) / 2;                          /* gete median element     */
	SWAP(indx[m], indx[l+1]);
	if (DR2(l+1) > DR2(r)) SWAP(indx[l+1], indx[r]);
	if (DR2(l  ) > DR2(r)) SWAP(indx[l  ], indx[r]);
	if (DR2(l+1) > DR2(l)) SWAP(indx[l+1], indx[l]);
                                    /* arrange as such  a(l+1) < a(l) < a(r) */
	i = l + 1;
	j = r;                                    /* initialize pointers     */
	while (j >= i) {                          /* scan to find elements,  */
	    do {i++;} while (DR2(i) < DR2(l));    /* larger than DR2(l)      */
	    do {j--;} while (DR2(j) > DR2(l));    /* smaller than DR2(l)     */
	    if (j >= i) SWAP(indx[i], indx[j]);   /* exchange them           */
	}
	SWAP(indx[l], indx[j]);                   /* put partition in active */
	if (j >= mth1) r = j - 1;                 /* set new partition which */
	if (j <= mth1) l = i;                     /* contains mth element    */
    }
    if (r-l == 1)                                 /* just two element remain */
	if (DR2(r) < DR2(l))                      /* if a(r) less than a(l)  */
	    SWAP(indx[r], indx[l]);               /* then exchange them      */
    for (i = 0; i < numngbr; i++)
	*q = nlst[indx[i]];                       /* copy rearranged list    */
    newh = rsqrt(DR2(mth1));                      /* get corrected h         */
    free(dr2);
    free(indx);
    free(nlst);
    return(newh);
}

#undef SWAP
#undef DR2

/*
 * KERNEL: arithmatic mean of kernels between two SPH particles.
 */

local real kernel(bodyptr p, bodyptr q, real dr2)
{
    real hinv, hinv2, coef;
    real ratio, ratio2, knl;
    bodyptr bptr[2];
    int i;

    bptr[0] = p;
    bptr[1] = q;                                 /* set pointers to p and q  */
    knl = 0.0;                                   /* initilize kernel value   */
    for (i=0; i<2; i++) {                        /* loop over two particles  */
	hinv = 1.0 / Hknl(bptr[i]);              /* compute inverse of h     */
	hinv2 = hinv * hinv;                     /* square of inverse of h   */
	ratio2 = dr2 * hinv2;                    /* and square of ratio of   */
                                                 /* separation and h         */
	if (ratio2 < 4.0) {                      /* if ratio is less than 2  */
#if defined(THREEDIM)
	    coef = INV_PI * hinv * hinv2;        /* compute coefficient (3-D)*/
#else
	    coef = TENSVN_IPI * hinv2;           /* compute coefficient (2-D)*/
#endif
	    ratio = rsqrt(ratio2);               /* ratio of separation & h  */
	    if (ratio2 <= 1.0 )                  /* if ratio is less than 1  */
		knl += coef * (1.0 - 0.75 * ratio2 * (2.0 - ratio));
	    else                                 /* else if ratio is in 1-2  */
		knl += coef * (0.25 * rpow((2.0 - ratio), 3.0));
	}
    }
    knl = 0.5 * knl;
    return(knl);
}

/*
 * GRADKERNEL: arithmatic mean of grad. of kernels between two SPH particles.
 */

local real gradkernel(bodyptr p, bodyptr q, real dr2)
{

    real hinv, hinv2, coef;
    real ratio, ratio2, gknl;
    bodyptr bptr[2];
    int i;

    bptr[0] = p;
    bptr[1] = q;                                 /* set pointers to p and q  */
    gknl = 0.0;                                  /* init grad kernel value   */
    for (i=0; i<2; i++) {                        /* loop over two particles  */
	hinv = 1.0 / Hknl(bptr[i]);              /* compute inverse of h     */
	hinv2 = hinv * hinv;                     /* square of inverse of h   */
	ratio2 = dr2 * hinv2;                    /* and square of ratio of   */
                                                 /* separation and h         */
	if (ratio2 < 4.0) {                      /* if ratio is less than 2  */
#if defined(THREEDIM)
	    coef = INV_PI * hinv * hinv2 * hinv2;/* compute coefficient (3-D)*/
#else
	    coef = TENSVN_IPI * hinv2 * hinv2;   /* compute coefficient (2-D)*/
#endif
	    ratio = rsqrt(ratio2);               /* ratio of separation & h  */
	    if (ratio <= TWOTRD)                 /* if ratio less than 2/3   */
		gknl += - coef / ratio;         
	    else if (ratio2 <= 1.0)              /* if ratio is less than 1  */
		gknl += coef * (-3.0 + 2.25 * ratio);
	    else                                 /* else if ratio is in 1-2  */
		gknl += coef * (-0.75) * rpow(2.0-ratio, 2.0) /ratio;
	}
    }
    gknl = 0.5 * gknl;
    return(gknl);
}

/*
 * CALCKERNEL: calculate mean values of karnel and its gradient.
 */

local void calckernel(bodyptr p, bodyptr q, real dr2, real *knl, real *gknl)
{

    real hinv, hinv2, coefk, coefg;
    real ratio, ratio2;
    bodyptr bptr[2];
    int i;

    bptr[0] = p;
    bptr[1] = q;                                 /* set pointers to p and q  */
    *knl = 0.0;                                  /* initialize kernel        */
    *gknl = 0.0;                                 /* and grad. of kernel      */
    for (i=0; i<2; i++) {                        /* loop over two particles  */
	hinv = 1.0 / Hknl(bptr[i]);              /* compute inverse of h     */
	hinv2 = hinv * hinv;                     /* square of inverse of h   */
	ratio2 = dr2 * hinv2;                    /* and square of ratio of   */
                                                 /* separation and h         */
	if (ratio2 < 4.0) {                      /* if ratio is less than 2  */
#if defined(THREEDIM)                            /* compute coefficient (3-D)*/
            coefk = INV_PI * hinv * hinv2;
	    coefg = coefk * hinv2;
#else                                            /* compute coefficient (2-D)*/
            coefk = TENSVN_IPI * hinv2;
	    coefg = coefk * hinv2;
#endif
	    ratio = rsqrt(ratio2);               /* ratio of separation & h  */
	    if (ratio <= TWOTRD) {               /* if ratio less than 2/3   */
		*knl += coefk * (1.0 - 0.75 * ratio2 * (2.0 - ratio));
		*gknl += - coefg / ratio;
	    } else if (ratio2 <= 1.0) {          /* if ratio is less than 1  */
		*knl += coefk * (1.0 - 0.75 * ratio2 * (2.0 - ratio));
 		*gknl += coefg * (-3.0 + 2.25 * ratio);
	    } else {                             /* else if ratio is in 1-2  */
		*knl += coefk * (0.25 * rpow(2.0 - ratio, 3.0));
 		*gknl += coefg * (-0.75) * rpow(2.0-ratio, 2.0) /ratio;
	    }
	}
    }
    *knl = 0.5 * *knl;
    *gknl = 0.5 * *gknl;
}

/*
 * CALCKERNEL2: calculate mean values of karnel gradient and its integral.
 */

local void calckernel2(bodyptr p, bodyptr q, real dr2, real *gknl, real *iknl)
{

    real hinv, hinv2, coefg, coefi;
    real ratio, ratio2, r22;
    bodyptr bptr[2];
    int i;

    bptr[0] = p;
    bptr[1] = q;                                 /* set pointers to p and q  */
    *gknl = 0.0;                                 /* and grad. of kernel      */
    *iknl = 0.0;
    for (i=0; i<2; i++) {                        /* loop over two particles  */

	hinv = 1.0 / Hknl(bptr[i]);              /* compute inverse of h     */
	hinv2 = hinv * hinv;                     /* square of inverse of h   */
	ratio2 = dr2 * hinv2;                    /* and square of ratio of   */
                                                 /* separation and h         */
	if (ratio2 < 4.0) {                      /* if ratio is less than 2  */
#if defined(THREEDIM)                            /* compute coefficient (3-D)*/
	    coefg = INV_PI * hinv * hinv2 * hinv2;
#else                                            /* compute coefficient (2-D)*/
	    coefg = TENSVN_IPI * hinv2 * hinv2;
#endif
	    ratio = rsqrt(ratio2);               /* ratio of separation & h  */
	    if (ratio <= 1.0) {
		if (ratio <= TWOTRD)             /* if ratio less than 2/3   */
		    *gknl += - coefg / ratio;
		else                             /* if ratio more than 2/3   */
		    *gknl += coefg * (-3.0 + 2.25 * ratio);
#if defined(THREEDIM)
		*iknl += ratio2 * ratio *
		    (FORTRD - 0.5 * ratio2 * (2.4 - ratio));
#else
		*iknl += ONESVN * ratio2 *
		    (10.0 - 3.0 * ratio2 * (2.5 - ratio));
#endif
	    } else {                             /* else if ratio is in 1-2  */
		r22 = rpow(2.0-ratio, 2.0);
 		*gknl += coefg * (-0.75) * r22 /ratio;
#if defined(THREEDIM)
		*iknl += 1.0 - ONESIX * r22 * r22 *
		    (r22 - 4.8 * (2.0-ratio) + 6.0) 
#else
		*iknl += ONESVN * (7.0 - (0.5 + ratio) * r22 * r22);
#endif
	    }
	} else {
	  *iknl += 1.0;
	}
    }
    *gknl = 0.5 * *gknl;
    *iknl = 0.5 * *iknl;
}
