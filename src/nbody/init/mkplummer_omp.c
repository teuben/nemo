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
 *
 *  OpenMP parallelization notes:
 *
 *  The main particle-generation loop in mkplummer() is embarrassingly
 *  parallel: each particle's position and velocity are sampled independently.
 *  The key challenge is that xrandom() uses global state and is not
 *  thread-safe.  We therefore replace it inside the parallel region with
 *  thread-local calls to drand48_r() (glibc) seeded per-thread.
 *
 *  Sections that remain serial:
 *    - mass-spectrum sampling via frandom() / getrfunc()  (uses global state)
 *    - the quiet==1/2 modes (depend on particle index ordering)
 *    - centre-of-mass correction (trivially parallelisable but tiny cost)
 *    - file I/O
 *
 *  To compile with OpenMP add -fopenmp to CFLAGS, e.g.:
 *    make CFLAGS="-O2 -fopenmp" mkplummer
 *  or edit makedefs.in / the local Makefile and add -fopenmp there.
 *
 *  Without -fopenmp the code compiles and behaves exactly as before
 *  (all #pragma omp directives are ignored by the preprocessor).
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

/* For thread-local RNG (drand48_r is POSIX / glibc) */
#include <stdlib.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

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
    "VERSION=3.0d-omp\n       16-may-2025 PJT + OpenMP",
    NULL,
};

string usage="construct a Plummer model";

/*---------------------------------------------------------------------------
 *  main  --  a tool to make a Plummer model, by invoking mkplummer()
 *            note:  only for 3 dimensions and Cartesian coordinates
 *---------------------------------------------------------------------------
 */
void nemo_main(void)
{
    bool    zerocm;
    int     nbody, seed, bits, quiet, i, j, n, nmodel, mode;
    real    snap_time, rfrac, mfrac, mlow, mrange[2], scale;
    Body    **btab, *bp;
    stream  outstr;
    string  massname;
    char    hisline[80], sseed[80];
    rproc   mfunc;
    rproc_body   r2func, v2func;
    Grid    gr2;
    int     ngr2, ir2;
    Moment  mgv2[MAXNGR2],  mgvz[MAXNGR2];
    Moment  mmgv2[MAXNGR2], mmgvz[MAXNGR2];
    real    v2, vz, r2, r2min, r2max;
      

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

    if (nbody < 1) error("Illegal number of bodies: %d",nbody);
    if (mfrac < 0 || mfrac > 1) error("Illegal mfrac=%g",mfrac);
    if (rfrac <= 0) error("Illegal rfrac=%g",rfrac);
    
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

    if (mode>0) outstr = stropen(getparam("out"), "w");

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
    
    if (mode > 0) {
      bits = (MassBit | PhaseSpaceBit | TimeBit);
      sprintf(hisline,"init_xrandom: seed used %d",seed);
      put_string(outstr, HeadlineTag, hisline);
      put_history (outstr);           /* update history */
      if (*headline)
	put_string (outstr, HeadlineTag, headline);
      for (i=0; i<nmodel; i++)
	put_snap (outstr, &btab[i], &nbody, &snap_time, &bits);
      strclose(outstr);
      if (mode==1) return;
      dprintf(0,"mode=2: data is stored, analysis now following\n");
    } else
      dprintf(0,"mode=0: no data stored, analysis now following\n");      

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
    }

    for (i=0; i<nmodel; i++) {
      for (j=0; j<ngr2; j++) {
	reset_moment(&mgv2[j]);
	reset_moment(&mgvz[j]);
      }

      /* btab[i] is now the root of each model */
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
      } /* j<ngr2 */

    }/* i<nmodel */

    /* final report for the bin statistics */
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


/*-----------------------------------------------------------------------------
 *  mkplummer  --  builds a nbody system according to a Plummer model.
 *
 *  OpenMP strategy
 *  ---------------
 *  The particle-generation loop (position + velocity sampling) is
 *  embarrassingly parallel.  Each iteration is fully independent.
 *
 *  The only complication is the random number generator.  NEMO's xrandom()
 *  uses a single global state and is NOT thread-safe.  We replace it inside
 *  the parallel section with drand48_r(), which takes a per-thread
 *  struct drand48_data buffer.  Each thread is seeded with
 *      base_seed * 1000003 + thread_id
 *  so that different threads produce different (but deterministic) sequences.
 *
 *  Helper macros replacing xrandom(lo,hi) and acos(xrandom(-1,1)):
 *    TL_DRAND(buf)        -> uniform [0,1)  using thread-local buf
 *    TL_XRANDOM(buf,lo,hi)-> uniform [lo,hi)
 *
 *  The mass assignment (frandom / mf path) and quiet==1/2 modes still call
 *  xrandom() and therefore remain serial (guarded with a critical section or
 *  left outside the parallel loop as appropriate).
 *
 *  The scaling pass and centre-of-mass correction are also parallelised with
 *  simple reduction clauses.
 *-----------------------------------------------------------------------------
 */

/* ---- thread-local RNG helpers (fall back to drand48 if not glibc) ---- */
#ifdef _OPENMP
#  ifdef __GLIBC__
     /* drand48_r is reentrant: each thread owns a struct drand48_data */
#    define USE_DRAND48_R 1
#  endif
#endif

#ifdef USE_DRAND48_R
#  define TL_DRAND(buf)          ({ double _r; drand48_r(&(buf), &_r); _r; })
#  define TL_XRANDOM(buf,lo,hi)  ((lo) + ((hi)-(lo)) * TL_DRAND(buf))
#else
   /* Fallback: plain drand48 with a global lock -- still correct but slower */
#  ifdef _OPENMP
#    define TL_DRAND(buf)          ({ double _r; _Pragma("omp critical(rng)") { _r = drand48(); } _r; })
#    define TL_XRANDOM(buf,lo,hi)  ((lo) + ((hi)-(lo)) * TL_DRAND(buf))
#  else
#    define TL_DRAND(buf)          drand48()
#    define TL_XRANDOM(buf,lo,hi)  xrandom(lo, hi)
#  endif
#endif

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
    real  radius=0.0;
    real  scalefactor;
    real  inv_scalefactor;
    real  sqrt_scalefactor;
    real  mrfrac;
    real  m_min, m_max;
    real  m_med;
    real  w_sum;
    vector w_pos, w_vel;
    Body  *btab;
    Body  *bp;
    real  vmax0 = 0.0;

    if (NDIM != 3)
        error("mkplummer: NDIM = %d but should be 3", NDIM);

    btab = (Body *) allocate ((size_t)nbody * sizeof(Body));

    scalefactor     = (scale < 0 ?  16.0 / (3.0 * PI)  : scale);
    inv_scalefactor = 1.0 / scalefactor;
    sqrt_scalefactor = sqrt( scalefactor );

    if (rfrac > 0) {
      rfrac *= scalefactor;
      mrfrac = rfrac*rfrac*rfrac / pow(1.0 + rfrac*rfrac, 1.5);
      if (mrfrac < mfrac)
	mfrac = mrfrac;
    } else
      warning("New feature: mfrac=%g\n",mfrac);

    /* ------------------------------------------------------------------ *
     *  Particle generation loop                                           *
     *                                                                     *
     *  quiet==0 is the only mode that is embarrassingly parallel.        *
     *  quiet==1/2 depend on the particle index for mass-shell bounds so  *
     *  they are left serial.  The mass-spectrum path (mf != NULL) calls  *
     *  frandom() which is also not thread-safe; left serial.             *
     * ------------------------------------------------------------------ */

    mtot = 0.0;

    if (quiet == 0 && mf == NULL) {
        /* ============================================================== *
         *  PARALLEL PATH                                                  *
         *  Each thread owns a drand48_r buffer seeded uniquely.          *
         * ============================================================== */

        /* Accumulate mtot and vmax0 via OpenMP reductions */
        real local_mtot  = 0.0;
        real local_vmax0 = 0.0;

#ifdef _OPENMP
        dprintf(1,"mkplummer: using %d OpenMP threads\n", omp_get_max_threads());
#endif

        #pragma omp parallel reduction(+:local_mtot) reduction(max:local_vmax0)
        {
            /* Per-thread RNG state */
#ifdef USE_DRAND48_R
            struct drand48_data rng_buf;
            int tid = 0;
#  ifdef _OPENMP
            tid = omp_get_thread_num();
#  endif
            /* Each thread gets a deterministic but distinct seed */
            srand48_r((long)(seed) * 1000003L + tid, &rng_buf);
#else
            /* Dummy declaration so TL_DRAND macros compile */
            int rng_buf = 0;
#  ifdef _OPENMP
            #pragma omp critical(rng)
#  endif
            srand48((long)seed);
#endif
            #pragma omp for schedule(dynamic, 64)
            for (i = 0; i < nbody; i++) {
                Body *bpi = btab + i;
                real vel, theta, phi, x, y, r;

                /* Equal mass for all particles (mf==NULL) */
                Mass(bpi) = 1.0 / (real) nbody;
                local_mtot += Mass(bpi);

                /* Position: invert cumulative mass-radius relation */
                r = 1.0 / sqrt( pow( TL_XRANDOM(rng_buf, mlow, mfrac),
                                     -2.0/3.0 ) - 1.0 );
                theta = acos( TL_XRANDOM(rng_buf, -1.0, 1.0) );
                phi   = TL_XRANDOM(rng_buf, 0.0, TWO_PI);
                Pos(bpi)[0] = r * sin(theta) * cos(phi);
                Pos(bpi)[1] = r * sin(theta) * sin(phi);
                Pos(bpi)[2] = r * cos(theta);

                /* Velocity: von Neumann rejection */
                x = 0.0;
                y = 0.1;
                while (y > x*x * pow(1.0 - x*x, 3.5)) {
                    x = TL_XRANDOM(rng_buf, 0.0, 1.0);
                    y = TL_XRANDOM(rng_buf, 0.0, 0.1);
                }
                vel = x * sqrt(2.0) * pow(1.0 + r*r, -0.25);
                if (vel > local_vmax0) local_vmax0 = vel;

                theta = acos( TL_XRANDOM(rng_buf, -1.0, 1.0) );
                phi   = TL_XRANDOM(rng_buf, 0.0, TWO_PI);
                Vel(bpi)[0] = vel * sin(theta) * cos(phi);
                Vel(bpi)[1] = vel * sin(theta) * sin(phi);
                Vel(bpi)[2] = vel * cos(theta);
            } /* end omp for */
        } /* end omp parallel */

        mtot  = local_mtot;
        vmax0 = local_vmax0;

    } else {
        /* ============================================================== *
         *  SERIAL PATH  (quiet!=0  OR  mass-spectrum requested)          *
         *  Identical to original code, uses NEMO's xrandom().            *
         * ============================================================== */
        for (i = 0, bp = btab; i < nbody; i++, bp++) {
            if (mf)
                Mass(bp) = frandom( mr[0], mr[1], mf );
            else
                Mass(bp) = 1.0 / (real) nbody;
            mtot += Mass(bp);

            if (quiet == 0)
                radius = 1.0 / sqrt( pow (xrandom(mlow,mfrac), -2.0/3.0) - 1.0);
            else if (quiet == 1) {
                m_min  = (i * mfrac) / nbody;
                m_max  = ((i+1) * mfrac) / nbody;
                radius = 1.0 / sqrt( pow (xrandom(m_min,m_max), -2.0/3.0) - 1.0);
            } else if (quiet == 2) {
                m_med  = ((i+0.5) * mfrac) / nbody;
                radius = 1.0 / sqrt( pow (m_med, -2.0/3.0) - 1.0);
            } else
                error("Illegal quiet=%d parameter\n", quiet);

            real theta_s = acos(xrandom(-1.0, 1.0));
            real phi_s   = xrandom(0.0, TWO_PI);
            Pos(bp)[0] = radius * sin(theta_s) * cos(phi_s);
            Pos(bp)[1] = radius * sin(theta_s) * sin(phi_s);
            Pos(bp)[2] = radius * cos(theta_s);

            real x_s = 0.0, y_s = 0.1;
            while (y_s > x_s*x_s * pow(1.0 - x_s*x_s, 3.5)) {
                x_s = xrandom(0.0, 1.0);
                y_s = xrandom(0.0, 0.1);
            }
            real velocity = x_s * sqrt(2.0) * pow(1.0 + radius*radius, -0.25);
            if (velocity > vmax0) vmax0 = velocity;

            theta_s = acos(xrandom(-1.0, 1.0));
            phi_s   = xrandom(0.0, TWO_PI);
            Vel(bp)[0] = velocity * sin(theta_s) * cos(phi_s);
            Vel(bp)[1] = velocity * sin(theta_s) * sin(phi_s);
            Vel(bp)[2] = velocity * cos(theta_s);
        }
    }

    dprintf(1,"Total mass (before scaling) = %g  vmax=%g\n", mtot, vmax0);

    /* ------------------------------------------------------------------ *
     *  Scale to VIRIAL units and normalise total mass to 1               *
     *  (trivially parallel -- kept as a separate pass)                   *
     * ------------------------------------------------------------------ */
    #pragma omp parallel for schedule(static)
    for (i = 0; i < nbody; i++) {
        bp = btab + i;
        Mass(bp) /= mtot;
        MULVS (Pos(bp), Pos(bp), inv_scalefactor);
        MULVS (Vel(bp), Vel(bp), sqrt_scalefactor);
    }

    /* ------------------------------------------------------------------ *
     *  Centre-of-mass correction                                         *
     * ------------------------------------------------------------------ */
    if (zerocm) {
        w_sum = 0.0;
        CLRV(w_pos);
        CLRV(w_vel);

        /* Reduction over position and velocity sums */
        real sum_px=0, sum_py=0, sum_pz=0;
        real sum_vx=0, sum_vy=0, sum_vz=0;
        real sum_w=0;

        #pragma omp parallel for schedule(static) \
            reduction(+:sum_px,sum_py,sum_pz,sum_vx,sum_vy,sum_vz,sum_w)
        for (i = 0; i < nbody; i++) {
            bp = btab + i;
            sum_w  += 1.0;
            sum_px += Pos(bp)[0];  sum_py += Pos(bp)[1];  sum_pz += Pos(bp)[2];
            sum_vx += Vel(bp)[0];  sum_vy += Vel(bp)[1];  sum_vz += Vel(bp)[2];
        }
        w_pos[0] = sum_px / sum_w;
        w_pos[1] = sum_py / sum_w;
        w_pos[2] = sum_pz / sum_w;
        w_vel[0] = sum_vx / sum_w;
        w_vel[1] = sum_vy / sum_w;
        w_vel[2] = sum_vz / sum_w;

        #pragma omp parallel for schedule(static)
        for (i = 0; i < nbody; i++) {
            bp = btab + i;
            SUBV(Pos(bp), Pos(bp), w_pos);
            SUBV(Vel(bp), Vel(bp), w_vel);
        }
    }

    return btab;
}

/* end of: mkplummer.c */
