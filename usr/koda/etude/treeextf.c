/****************************************************************************/
/* TREEEXTF.C: routines to compute external force.                          */
/* Public routines: initextf(), extforce(), massext()                       */
/* Copyright (c) 2000 by Jin Koda, Tokyo, JAPAN.                            */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "treecode.h"

//#define SK00                  /* Sakamoto, Baker, Scoville (2000) potential */

/* Local routines and variables to calculate external force. */

local real potent1(real);
local real dpot0dr_r(real);
local real dpot1dr_r(real);

static real cfpot;                              /* coef. of external pot.   */
static real acore;                              /* core radius in sys unit  */
static real acorebh;                            /* core radius of BH        */
static real acoresq;                            /* square of core radius    */
static real acorebhsq;                          /* square of BH core radius */

/*
 * INITEXTF: initialize external force parameters.
 */

void initextf(void)
{

    real vmsys;

    acore = rcore / rscale;                      /* core radius in sys unit  */
    acorebh = rcorebh / rscale;                  /* BH radius in sys unit    */
    acoresq = acore * acore;                     /* sqr. of core radius      */
    acorebhsq = acorebh * acorebh;               /* sqr. of BH core radius   */
    vmsys = vmax /rscale *tscale *1.022712169e-9;/* vmax in system unit      */

#if defined(SK00)
    cfpot = 0.5 * vmsys * vmsys;
#else
    cfpot = rsqrt(27.0/4.0) * acore * vmsys * vmsys;
#endif
    
}

/*
 * EXTFORCE: compute external force.
 */

void extforce(bodyptr btab, int nbody, real time)
{
    bodyptr p;
    real abh, abhsq, a, asq, r, rsq, omgnow, accaxi, accbar;
    real c1, c2, dot, crs, sinmt, cosmt, sinp, cosp, ctptch;
    vector acc0, nlon;

    a = acore;                                  /* core radius in sys unit  */
    asq = a * a;                                /* square of it in sys.unit */
    abh = acorebh;                              /* core radius in sys unit  */
    abhsq = abh * abh;                          /* and square               */

    if (pitch < 90.0)                           /* if spiral potential      */
      ctptch = mode / rtan(pitch / 180.0 * PI); /* ctan of pitch angle      */
    else                                        /* if bar potential is used */
      ctptch = 0.0;                             /* switch off spiral term   */

    omgnow = (1.02273e-9 *omgb) *(time *tscale);/* compute angle, and       */
    nlon[0] = rcos(omgnow);                     /* direction of the line-of-*/
    nlon[1] = rsin(omgnow);                     /* node                     */

    for (p = btab; p < btab+nbody; p++) {       /* loop over particles      */

	CLRV(acc0);                             /* clear acceleration       */
	ABSV(r, Pos(p));                        /* compute radus, and       */
	rsq = r * r;                            /* square of radius         */
	accaxi = - cfpot *                      /* acceleration due to      */
	    (dpot0dr_r(rsq) +                   /* axisymmetric potential   */
	     fbh / rpow(rsq + abhsq, 1.5));     /* and central black hole   */
	MULVS(acc0, Pos(p), accaxi);            

	DOTVP(dot, Pos(p), nlon);               /* dot and cross products   */
	crs = nlon[0] * Pos(p)[1] - nlon[1] * Pos(p)[0];
                                                /* of position with bar vec.*/

	cosmt = 2.0 * rpow(dot / r, 2.0) - 1.0; /* cos & sin(m theta)       */
	sinmt = 2.0 * crs * dot / rsq;          /* the case for m=2         */

	cosp = cosmt * cos(ctptch * rlog(r)) - sinmt * sin(ctptch * rlog(r));
	sinp = sinmt * cos(ctptch * rlog(r)) + cosmt * sin(ctptch * rlog(r));
	c1 = dpot1dr_r(rsq) * cosp - ctptch * potent1(rsq) / rsq * sinp;
	c2 = mode * potent1(rsq) / rsq * sinp;

	accbar = - fbar * cfpot;
	acc0[0] += accbar * (c1 * Pos(p)[0] + c2 * Pos(p)[1]);
	acc0[1] += accbar * (c1 * Pos(p)[1] - c2 * Pos(p)[0]);
	ADDV(Acc(p), Acc(p), acc0);

    }
}

/*
 * MASSEXT: mass of external potential within a radius [in unit of Msun].
 * input radius should be in the unit of kpc.
 */

global real massext(real rad)
{
    real ms;
    
#if defined(SK00)
    ms =  2.325138e5 * rcore * vmax * vmax *
	rpow(rad, 3.0) / (rcore * (rad*rad + rcore*rcore));
#else
    ms =  2.325138e5 * rsqrt(27.0/4.0) * rcore * vmax * vmax *
	rpow(rad, 3.0) / rpow(rad*rad + rcore*rcore, 1.5);
#endif

    return(ms);
}

/*
 * DPOT0DR_R: first derivative of axisymmetric potential, devided by r.
 */

local real dpot0dr_r(real rsq)
{
    real dpr0;

#if defined(SK00)
    dpr0 = 2.0 / (rsq + acoresq);
#else
    dpr0 = 1.0 / rpow(rsq + acoresq, 1.5);
#endif

    return(dpr0);
}  

/*
 * DPOT1DR_R: first derivative of bisymmetric potential, devided by r.
 */

local real dpot1dr_r(real rsq)
{
    real dpr1;

#if defined(SK00)
    dpr1 = -2.0 / (rsq + acoresq);
#else
    dpr1 = -2.0 * acore * (acoresq - rsq) / rpow(rsq + acoresq, 3.0);
#endif

    return(dpr1);
}

/*
 * POTENT1: bisymmetric potential.
 */

local real potent1(real rsq)
{
    real bipot;

#if defined(SK00)
    bipot = - rlog(1.0 + (rsq/acoresq));
#else
    bipot = - acore * rsq / rpow(rsq + acoresq, 2.0);
#endif
    return(bipot);
}

