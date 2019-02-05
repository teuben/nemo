/* mkplummer.c - mkplummer, main<TOOL> */

/*
 *  mkplummer.c: constructs a Plummer model, with a spatial or mass cut-off
 *		(c) 1987  Piet Hut  Princeton, NJ, USA
 *			
 *	10-Mar-89	V2.0 modified for new filestruct - Peter Teuben
 *	10-Jul-89	V2.1 headline also output		PJT
 *	17-sep-90	V2.2 experiment with mass spectrum	PJT
 *      15-nov-90       V2.3 default is no mass spectrum - fixed bug PJT
 *	15-apr-91	V2.4 used set_xrandom's return value	PJT
 *	 7-mar-92	V2.5 happy gcc2.0			pjt
 *	19-mar-92	fixed wrong allocate decl. 
 *      20-may-92       V2.5a SGI
 *	17-feb-94       V2.5b fixed for -DSINGLEPREC		pjt
 *	11-jan-95	V2.5c ??? (or was it march 1994)
 *	 6-jun-96       V2.6d report total mass before rescaling  pjt
 *       8-sep-01       e   init_xrandom
 *      19-oct-01       V3.0  add a scale= parameter (defaults to virial units) pjt
 *      22-mar-04       V2.7  merged version with a hole      ncm+pjt
 *      31-mar-05       V2.8  added nmodel=                       pjt
 *      30-may-07       V2.8b allocate() with size_t
 */


#include  <stdinc.h>
#include  <getparam.h>
#include  <vectmath.h>
#include  <filestruct.h>
#include  <history.h>
#include  <loadobj.h>

#include <snapshot/snapshot.h>  
#include <snapshot/body.h>
#include <snapshot/put_snap.c>
#include <bodytransc.h>

#include <moment.h>
#include <grid.h>

extern rproc  getrfunc();
Body    *mkplummer();

local string headline;		/* random text message */

#define MAXNGR2 100


string  defv[] = {                        /* DEFAULT INPUT PARAMETERS */
    "out=???\n		      Output file name",
    "nbody=???\n	      Number of particles",
    "mlow=0\n                 Low mass fraction cutoff of Plummer dist",
    "mfrac=0.999\n            Mass fraction used of Plummer distribution",
    "rfrac=22.8042468\n       Radius fraction used of Plummer distribution\n\
                              NOTE: the above two values are chosen so\n\
                                    that m( rfrac ) = mfrac                ",
    "seed=0\n                 Seed for the random number generator",
    "time=0.0\n               Time at which snapshot is taken",
    "zerocm=t\n               Centrate snapshot (t/f)?",
    "scale=-1\n               Model scale factor (-1=virial 1=natural)",
    "quiet=0\n                0=noisy 1=somewhat quiet 2=more quiet",
    "massname=\n              If used Mass-function name (e.g. n(m))",
    "massexpr=pow(m,p)\n      Mass function expression (e.g. pow(m,p))",
    "masspars=p,0.0\n         Mass function parameters (e.g. p,0.0)",
    "massrange=1,1\n          Range for mass-spectrum (e.g. 1,2)",
    "headline=\n	      Verbiage for output",
    "nmodel=1\n               number of models to produce",
    "mode=1\n                 0=no data,  1=data, no analysis 2=data, analysis",
    "single=true\n            Single output file?",
    "VERSION=3.0\n            7-nov-2018 PJT",
    NULL,
};

string usage="construct a Plummer model";

string cvsid="$Id$";

/*---------------------------------------------------------------------------
 *  main  --  a tool to make a Plummer model, by invoking mkplummer()
 *            note:  only for 3 dimensions and Carthesian coordinates
 *---------------------------------------------------------------------------
 */
void nemo_main(void)
{
    bool    zerocm, single;
    int     nbody, seed, bits, quiet, i, j, n, nmodel, mode;
    real    snap_time, rfrac, mfrac, mlow, mrange[2], scale, rsum;
    Body    **btab, *bp;
    stream  outstr;
    string  massname;
    char    hisline[80], sseed[80], fname[80];
    rproc   mfunc;
    rproc_body   r2func, v2func, vzfunc;
    Grid    gr2;
    int     ngr2, ir2;
    Moment  mgv2[MAXNGR2],  mgvz[MAXNGR2];
    Moment  mmgv2[MAXNGR2], mmgvz[MAXNGR2];
    real    v2, vz, r2, r2min, r2max;
    real    t0, t1, t2, t3, t4;
      

    t0 = cputime();                    // start
    nbody = getiparam("nbody");
    mfrac = getdparam("mfrac");
    mlow  = getdparam("mlow");
    rfrac = getdparam("rfrac");
    seed = init_xrandom(getparam("seed"));
    snap_time = getdparam("time");
    zerocm = getbparam("zerocm");
    quiet = getiparam("quiet");
    scale = getdparam("scale");
    headline = getparam("headline");
    massname = getparam("massname");
    nmodel = getiparam("nmodel");
    mode = getiparam("mode");
    single = getbparam("single");
    if (*massname) {
        mysymbols(getargv0());
        n=1;
        mfunc = getrfunc(massname,getparam("massexpr"),getparam("masspars"),&n);
        if (nemoinpr(getparam("massrange"),mrange,2)!=2)
            error("Need two numbers in massrange=");
    	dprintf(1,"Massrange from %f : %f\n",mrange[0],mrange[1]);
    } else
        mfunc = NULL;

    r2func = btrtrans("r2");
    v2func = btrtrans("v2");
    vzfunc = btrtrans("vz");

    if (mode>0 && single) outstr = stropen(getparam("out"), "w");

    t1 = cputime();   // init
    
    btab = (Body **) allocate(sizeof(Body **) * nmodel);
    for (i=0; i<nmodel; i++) {
      if (i>0) {
	seed++;
	sprintf(sseed,"%d",seed);
	init_xrandom(sseed);
      }
      btab[i] = mkplummer(nbody, mlow, mfrac, rfrac, seed, snap_time, zerocm, scale,
		       quiet,mrange,mfunc);
    }
    t2 = cputime(); // models created now
    if (mode > 0) {
      if (single) {
	bits = (MassBit | PhaseSpaceBit | TimeBit);
	sprintf(hisline,"init_xrandom: seed used %d",seed);
	put_string(outstr, HeadlineTag, hisline);
	put_history (outstr);           /* update history */
	if (*headline)
	  put_string (outstr, HeadlineTag, headline);
	for (i=0; i<nmodel; i++)
	  put_snap (outstr, &btab[i], &nbody, &snap_time, &bits);
	strclose(outstr);
      } else {
	for (i=0; i<nmodel; i++) {
	  sprintf(fname,"%s.%d",getparam("out"),i);
	  outstr = stropen(fname, "w");
	  put_snap (outstr, &btab[i], &nbody, &snap_time, &bits);
	  strclose(outstr);
	}
      }
    } 

    t3 = cputime(); // i/o done

    if (mode != 1) {
      dprintf(0,"analysis now following\n");
      /* match the default grid in cluster_stats.py */
      r2min = 0.0;
      r2max = 3.875;
      ngr2  = 31;

      inil_grid(&gr2,ngr2,r2min,r2max);
      for (j=0; j<ngr2; j++) {
	ini_moment(&mgv2[j],  1, nbody);
	ini_moment(&mgvz[j],  1, nbody);
	ini_moment(&mmgv2[j], 2, nbody);
	ini_moment(&mmgvz[j], 2, nbody);
	// printf("%d %g\n", j, value_grid(&gr2,j));
      }

      for (i=0; i<nmodel; i++) {
	rsum = 0.0;
	for (j=0; j<ngr2; j++) {
	  reset_moment(&mgv2[j]);
	  reset_moment(&mgvz[j]);
	}

	// btab[i] is now the root of each model
	for (j = 0, bp=btab[i]; j < nbody; j++, bp++) {
	  r2  = r2func(bp,0.0,j);
	  ir2 = index_grid(&gr2, r2);
	  
	  v2 = v2func(bp,0.0,j);
	  vz = Vel(bp)[2];	
	  if (ir2 >= 0 && ir2 < ngr2) {
	    accum_moment(&mgv2[ir2], v2, 1.0);
	    accum_moment(&mgvz[ir2], vz, 1.0);
	  }
	  
	}
	for (j=0; j<ngr2; j++) {
	  if (n_moment(&mgv2[j]) > 0) {
	    accum_moment(&mmgv2[j], mean_moment(&mgv2[j]), 1.0);
	    accum_moment(&mmgvz[j], mean_moment(&mgvz[j]), 1.0);
	  }
	} //j<ngr2

      }//i<nmodel

      // final report for the bin statistics
    
      for (j=0; j<ngr2; j++) {
	dprintf(1,"j=%d %g  %g %g   %g %g n=%d %d\n",
	     j,
	     value_grid(&gr2,j+0.5),
	     mean_moment(&mmgv2[j]),
	     sigma_moment(&mmgv2[j]),
	     mean_moment(&mmgvz[j]),
	     sigma_moment(&mmgvz[j]),
	     n_moment(&mmgv2[j]),
	     n_moment(&mmgvz[j]));
      }
    }
    t4 = cputime(); // analysis done
    dprintf(0,"%g %g %g %g\n",(t1-t0)*60,(t2-t1)*60,(t3-t2)*60,(t4-t3)*60);
}


/*-----------------------------------------------------------------------------
 *  mkplummer  --  builds a nbody system according to a Plummer model,
 *                 in VIRIAL units (M=G=-4E=1, with E the total energy),
 *                 and finite spatial extent which can be regulated by
 *                 specifying mfrac or rfrac or using their default values.
 *		   note: the distribution function is spherical and isotropic,
 *		         and is a polytrope of index n = 5.
 *		   litt: S.J. Aarseth, M. Henon and R. Wielen (1974),
 *		         Astron. and Astrophys. 37, p. 183.
 *                 accepts: nbody: the number of particles.
 *                          mlow:  lower mass fraction cutoff
 *                          mfrac: mass fraction of the (infinitely extended)
 *                                 Plummer model; see  rfrac immediately below.
 *                          rfrac: radius fraction of the (infinitely extended)
 *                                 Plummer model. If mfrac = rfrac = 1.0 then 
 *                                 particles will be sprinkled in all over
 *                                 space. If mfrac < 1.0 or rfrac < 1.0 then 
 *                                 each particle is constrained to lie within
 *                                 both the radial and (cumulative) mass bound.
 *                                 For example, if rfrac( mfrac ) > rfrac then
 *                                 rfrac is the limiting factor, but if
 *                                 rfrac( mfrac ) < rfrac then mfrac limits 
 *                                 the extent of the Plummer realization.
 *                                 Note: specifying either value may have no
 *                                       effect if the default value of the
 *                                       other parameter is still the limiting
 *                                       factor; Beware!
 *			    seed: the seed for the random number generator;
 *                          snap_time: the time at which the snapshot applies.
 *                          zerocm: logical determining if to center snapshot
 *                          quiet: integer how quiet model should be (0=noisy)
 *                 returns: snap: a pointer to the new snapshot, containing
 *			          a Plummer model in which all particles have
 *			          equal masses.
 *                 NOTE: after sprinkling in particles according to a Plummer
 *                       distribution, the whole system is shifted in position
 *                       and velocity so as to put the center of mass at rest
 *                       at the coordinate center.
 *-----------------------------------------------------------------------------
 */

Body *mkplummer(nbody, mlow, mfrac, rfrac, seed, snap_time,zerocm,scale,quiet,mr,mf)
int   nbody;
real  mfrac;
real  mlow;
real  rfrac;
int   seed;
real  snap_time;
bool  zerocm;
real  scale;
int   quiet;
real  mr[2];
rproc mf;
{
    register int  i;
    real  mtot;
    real  radius;		/* absolute value of position vector      */
    real  velocity;		/* absolute value of velocity vector      */
    real  theta, phi;		/* direction angles of above vectors      */
    real  x, y;		        /* for use in rejection technique         */
    real  scalefactor;          /* for converting between different units */
    real  inv_scalefactor;      /* inverse scale factor                   */
    real  sqrt_scalefactor;     /* sqare root of scale factor             */
    real  mrfrac;               /* m( rfrac )                             */
    real  m_min, m_max;         /* mass shell limits for quiet=1          */
    real   m_med;		/* mass shell value for quiet=2           */
    real   w_sum;               /* temporary storage for c.o.m. calc      */
    vector w_pos, w_vel;        /* temporary storage for c.o.m. calc      */
    Body  *btab;                /* pointer to the snapshot                */
    Body  *bp;                  /* pointer to one particle                */

    if (NDIM != 3)
        error("mkplummer: NDIM = %d but should be 3", NDIM);

    btab = (Body *) allocate ((size_t)nbody * sizeof(Body));

/*
 *  Calculating the coordinates is easiest in STRUCTURAL units;
 *  conversion to VIRIAL units will be performed below.
 *
 *    Recipe for scaling to the proper system of units:
 *
 *  Since G = M = 1, if we switch from a coordinate system with
 *  length unit  r_old  to a coordinate system with length unit  r_new ,
 *  the length units simply scale by a factor  C = r_new / r_old .
 *  Consequently, the coordinate values of physical quantities
 *  such as positions should transform inversely to maintain the same
 *  coordinate-invariant meaning. Similarly, the square of the velocities
 *  should transform inversely proportional to the positions,
 *  since  GM = 1  (cf. a relation such as  v*v = G*M/r ).
 *  To summarize: If
 *                       r_unit(new) = C * r_unit(old)  ,
 *                then
 *                       pos(new) = (1/C) * pos(old)
 *                and
 *                       vel(new) = sqrt(C) * vel(old)  .
 */
    scalefactor = (scale < 0 ?  16.0 / (3.0 * PI)  : scale);
    inv_scalefactor = 1.0 / scalefactor;
    sqrt_scalefactor = sqrt( scalefactor );
/*
 *  we can now convert  rfrac  into an equivalent mfrac:
 *  but by allowing rfrac <= 0.0 we can cheat and use mfrac=1 (for example)
 */
    if (rfrac > 0) {
      rfrac *= scalefactor;          /* from VIRIAL to STRUCTURAL units */
      mrfrac = rfrac*rfrac*rfrac / pow(1.0 + rfrac*rfrac, 1.5);
      if (mrfrac < mfrac)
	mfrac = mrfrac;            /* mfrac = min(mfrac, m(rfrac)) */
    } else
      warning("New feature: mfrac=%g\n",mfrac);
      
/*
 *  now we construct the individual particles:
 */
    mtot = 0.0;
    for (i = 0, bp=btab; i < nbody; i++, bp++) {
	if (mf)     /* if mass spectrum given: */
	    Mass(bp) = frandom( mr[0], mr[1], mf );
        else        /* else all stars equal mass */
            Mass(bp) = 1.0/ (real) nbody;
	mtot += Mass(bp);
/*
 *  the position coordinates are determined by inverting the cumulative
 *  mass-radius relation, with the cumulative mass drawn randomly from
 *  [0, mfrac]; cf. Aarseth et al. (1974), eq. (A2).
 */
        if (quiet==0)
	    radius = 1.0 / sqrt( pow (xrandom(mlow,mfrac), -2.0/3.0) - 1.0);
        else if (quiet==1) {
            m_min = (i * mfrac)/nbody;
            m_max = ((i+1) * mfrac)/nbody;
            radius = 1.0 / sqrt( pow (xrandom(m_min,m_max), -2.0/3.0) - 1.0);
        } else if (quiet==2) {
            m_med = ((i+0.5) * mfrac)/nbody;
            radius = 1.0 / sqrt( pow (m_med, -2.0/3.0) - 1.0);
	} else	
	    error("Illegal quiet=%d parameter\n",quiet);
	theta = acos(xrandom(-1.0, 1.0));
	phi = xrandom(0.0, TWO_PI);
	Pos(bp)[0] = radius * sin( theta ) * cos( phi );
	Pos(bp)[1] = radius * sin( theta ) * sin( phi );
        Pos(bp)[2] = radius * cos( theta );
/*
 *  the velocity coordinates are determined using von Neumann's rejection
 *  technique, cf. Aarseth et al. (1974), eq. (A4,5).
 *  First we take initial values for x, the ratio of velocity and escape
 *  velocity (q in Aarseth et al.), and y, as a trick to enter the body of the
 *  while loop.
 */
	x = 0.0;
	y = 0.1;
/*
 *  Then we keep spinning the random number generator until we find a pair
 *  of values (x,y), so that y < g(x) = x*x*pow( 1.0 - x*x, 3.5) . Whenever
 *  an y-value lies above the g(x) curve, the (x,y) pair is discarded, and
 *  a new pair is selected. The value 0.1 is chosen as a good upper limit for
 *  g(x) in [0,1] : 0.1 > max g(x) = 0.092 for 0 < x < 1.
 */
	while (y > x*x*pow( 1.0 - x*x, 3.5)) {
	    x = xrandom(0.0,1.0);
	    y = xrandom(0.0,0.1);
        }
/*
 *  If y < g(x), proceed to calculate the velocity components:
 */
	velocity = x * sqrt(2.0) * pow( 1.0 + radius*radius, -0.25);
	theta = acos(xrandom(-1.0, 1.0));
	phi = xrandom(0.0,TWO_PI);
	Vel(bp)[0] = velocity * sin( theta ) * cos( phi );
	Vel(bp)[1] = velocity * sin( theta ) * sin( phi );
	Vel(bp)[2] = velocity * cos( theta );
    }
    dprintf(1,"Total mass (before scaling) = %g\n",mtot);
/*
 * Now transform to the VIRIAL coordinates by applying
 * the scaling factors to the positions and velocities:
 *
 * Also normalize total mass to 1
 */
    for (i = 0, bp=btab; i < nbody; i++, bp++) {
        Mass(bp) /= mtot;
        MULVS (Pos(bp), Pos(bp), inv_scalefactor);
        MULVS (Vel(bp), Vel(bp), sqrt_scalefactor);
    }

    if (zerocm) {       /* False for Masspectrum */
        w_sum = 0.0;
        CLRV(w_pos);
        CLRV(w_vel);
        for (i = 0, bp = btab; i < nbody; i++, bp++) {
            w_sum = w_sum + 1.0;        /* all bodies same mass */
            ADDV(w_pos, w_pos, Pos(bp));
            ADDV(w_vel, w_vel, Vel(bp));
        }
        DIVVS(w_pos, w_pos, w_sum);
        DIVVS(w_vel, w_vel, w_sum);
        for (i = 0, bp = btab; i < nbody; i++, bp++) {
            SUBV(Pos(bp), Pos(bp), w_pos);
            SUBV(Vel(bp), Vel(bp), w_vel);
        }
    }

    return btab;               /* snap out of it altogether */
}

/* end of: mkplummer.c */




