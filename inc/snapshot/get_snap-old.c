/*
 * GET_SNAP.C: generic method for snapshot input.
 *	Note: this file is to be included at the source level 
 *	      it cannot be included, since it needs
 *		<stdinc.h> <snapshot.h> and <body.h>
 * Joshua E. Barnes  23-May-88  Princeton NJ.
 *	18-nov-91 malloc() -> allocate()
 *	 7-mar-92 happy GCC2.0
 *	 2-may-92  fixed allocate() declaration PJT
 *	22-feb-94 ansi headers (w/ allocate)    pjt
 *	26-jun-96 no more local definitions, use extern		PJT
 *       8-oct-01 read eps/dens                 PJT
 *      17-jan-02 detect split Pos/Vel		pjt
 *       2-apr-02 add UdotIntTag for ZENO	pjt
 *      30-may-07 allocate() needs size_t args for > 44.7M      pjt
 *    14-feb-2017 added get_snap_nbody()                        pjt
 */

/*
 * get_snap is a generic method for reading snapshot data from a file,
 * to be #included in the compilation of an application program:
 *
 *	#define Body      ...
 *	#define Mass(b)   ...
 *	#define Phase(b)  ...
 *
 *	#include <snapshot/snapshot.h>
 *	#include <snapshot/get_snap.c>
 *
 * The simplest way to use get_snap is probably as follows:
 *
 *	stream instr;
 *	Body *btab = NULL;
 *	int nbody, bits;
 *	real tsnap;
 *
 *	get_snap(instr, &btab, &nbody, &tsnap, &bits);
 *	if ((bits & <required_bits>) == <required_bits>)
 *	    <process_input_data>;
 *
 * Notes: (1) if *btab == NULL, a new body array of length nbody will be
 * allocated, which the programmer may access when get_snap returns, while
 * if *btab != NULL, then it must point to an array of length nbody, and
 * the new data in the input stream replaces the data in the body table.
 * (2) bits is the logical OR of the bit flags for the individual
 * components, defined in snapshot/snapshot.h. (3) The name get_snap is
 * a macro, defined by get_snap.c, so the include statement must come
 * before the first usage. (4) The vanilla get_snap or any subsidiary
 * routine may be replaced by giving the macro name a definition before
 * including get_snap.c.  For example, to supply your own definition of
 * the mass-input routine, do the following:
 *
 *	#define get_mass  my_get_mass
 *
 *	#include <snapshot/get_snap.c>
 *
 *	local my_get_mass(....) {....}
 *
 * Look at the definition of the standard function(s) you are replacing to
 * find out what arguments to expect.
 */

/*
 * GET_SNAP_PARAMETERS: worker routine to input snapshot parameters.
 */


#ifndef get_snap_parameters

#define get_snap_parameters  _get_snap_parameters

local void
_get_snap_parameters(instr, btptr, nbptr, tsptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
real *tsptr;			/* pointer to time of input */
int *ifptr;			/* pointer to input bit flags */
{
    int nbody;

    if (get_tag_ok(instr, ParametersTag)) {
	get_set(instr, ParametersTag);
	if (get_tag_ok(instr,NobjTag))
	  get_data(instr, NobjTag, IntType, &nbody, 0);
	else if (get_tag_ok(instr,NBodyTag)) {
	  get_data(instr, NBodyTag, IntType, &nbody, 0);
	  warning("Reading a ZENO file with NBody=%d",nbody);
	} else
	  error("Cannot find Nobj or NBody in snapshot");
	if (*btptr != NULL && nbody > *nbptr)	/* bigger than expected? */
	    error("get_snap_parameters: %s = %d is too big now %d\n",
		  NobjTag, nbody, *nbptr);
	*nbptr = nbody;				/* set input value */
	if (get_tag_ok(instr, TimeTag)) {	/* time data specified? */
	    get_data_coerced(instr, TimeTag, RealType, tsptr, 0);
	    *ifptr |= TimeBit;			/*   set time bit */
	}
	get_tes(instr, ParametersTag);
    }
}

#endif

/*
 * GET_SNAP_CSYS: worker routine to check coordinate system.
 */

#ifndef get_snap_csys

#define get_snap_csys  _get_snap_csys

local void
_get_snap_csys(instr, ifptr)
stream instr;			/* input stream, of course */
int *ifptr;			/* pointer to input bit flags */
{
    int cs;

    if (get_tag_ok(instr, CoordSystemTag)) {
	get_data(instr, CoordSystemTag, IntType, &cs, 0);
	if (cs != CSCode(Cartesian, NDIM, 2))
	    error("get_snap_csys: cant handle %s = %#o\n",
		  CoordSystemTag, cs);
    } else
	dprintf(1,"get_snap_csys: assuming %s = %#o\n",
	       CoordSystemTag, CSCode(Cartesian, NDIM, 2));
}

#endif

/*
 * GET_SNAP_MASS: worker routine to input mass data.
 */

#ifndef get_snap_mass

#define get_snap_mass  _get_snap_mass

local void
_get_snap_mass(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
#ifdef Mass
    real *mbuf, *mp;
    Body *bp;

    if (get_tag_ok(instr, MassTag)) {
        mbuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	get_data_coerced(instr, MassTag, RealType, mbuf, *nbptr, 0);
	for (bp = *btptr, mp = mbuf; bp < *btptr + *nbptr; bp++)
	    Mass(bp) = *mp++;
	free(mbuf);
	*ifptr |= MassBit;
    }
#endif
}

#endif

/*
 * GET_SNAP_PHASE: worker routine to input phasespace data.
 */

#ifndef get_snap_phase

#define get_snap_phase  _get_snap_phase

local void
_get_snap_phase(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
#ifdef Phase
    real *rvbuf, *rvp;
    Body *bp;

    if (get_tag_ok(instr, PhaseSpaceTag)) {
        rvbuf = (real *) allocate((size_t)(*nbptr) * 2 * NDIM * sizeof(real));
	get_data_coerced(instr, PhaseSpaceTag, RealType, rvbuf,
			 *nbptr, 2, NDIM, 0);
	for (bp = *btptr, rvp = rvbuf; bp < *btptr + *nbptr; bp++) {
	    SETV(Phase(bp)[0], rvp);
	    rvp += NDIM;
	    SETV(Phase(bp)[1], rvp);
	    rvp += NDIM;
	}
	free(rvbuf);
	*ifptr |= PhaseSpaceBit;
    } else if (get_tag_ok(instr, PosTag) || get_tag_ok(instr, VelTag)) {
      real *rbuf, *vbuf, *rp, *vp;
      rbuf = (real *) allocate((size_t)(*nbptr) * 2 * NDIM * sizeof(real));
      vbuf = (real *) allocate((size_t)(*nbptr) * 2 * NDIM * sizeof(real));
      if (get_tag_ok(instr,PosTag))
	  get_data_coerced(instr, PosTag, RealType, rbuf,
		       *nbptr, NDIM, 0);
      if (get_tag_ok(instr,VelTag))
	  get_data_coerced(instr, VelTag, RealType, vbuf,
		       *nbptr, NDIM, 0);
      for (bp = *btptr, rp=rbuf, vp=vbuf; bp < *btptr + *nbptr; bp++) {
	SETV(Phase(bp)[0], rp);
	rp += NDIM;
	SETV(Phase(bp)[1], vp);
	vp += NDIM;
      }
      free(rbuf);
      free(vbuf);
      *ifptr |= PhaseSpaceBit;
    }
#endif
}

#endif

/*
 * GET_SNAP_PHI: worker routine to input potential data.
 */

#ifndef get_snap_phi

#define get_snap_phi  _get_snap_phi

local void
_get_snap_phi(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
#ifdef Phi
    real *pbuf, *pp;
    Body *bp;

    if (get_tag_ok(instr, PotentialTag)) {
        pbuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	get_data_coerced(instr, PotentialTag, RealType, pbuf, *nbptr, 0);
	for (bp = *btptr, pp = pbuf; bp < *btptr + *nbptr; bp++)
	    Phi(bp) = *pp++;
	free(pbuf);
	*ifptr |= PotentialBit;
    }
#endif
}

#endif

/*
 * GET_SNAP_ACC: worker routine to input accelerations.
 */

#ifndef get_snap_acc

#define get_snap_acc  _get_snap_acc

local void
_get_snap_acc(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
#ifdef Acc
    real *abuf, *ap;
    Body *bp;

    if (get_tag_ok(instr, AccelerationTag)) {
        abuf = (real *) allocate((size_t)(*nbptr) * NDIM * sizeof(real));
	get_data_coerced(instr, AccelerationTag, RealType, abuf,
			 *nbptr, NDIM, 0);
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++) {
	    SETV(Acc(bp), ap);
	    ap += NDIM;
	}
	free(abuf);
	*ifptr |= AccelerationBit;
    }
#endif
}

#endif

/*
 * GET_SNAP_AUX: worker routine to input aux values.
 */

#ifndef get_snap_aux

#define get_snap_aux  _get_snap_aux

local void
_get_snap_aux(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
#ifdef Aux
    real *abuf, *ap;
    Body *bp;

    if (get_tag_ok(instr, AuxTag)) {
        abuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	get_data_coerced(instr, AuxTag, RealType, abuf, *nbptr, 0);
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++)
	    Aux(bp) = *ap++;
	free(abuf);
	*ifptr |= AuxBit;
    }
#endif
}

#endif

/*
 * GET_SNAP_KEY: worker routine to input key values.
 */

#ifndef get_snap_key

#define get_snap_key  _get_snap_key

local void
_get_snap_key(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
#ifdef Key
    int *kbuf, *kp;
    Body *bp;

    if (get_tag_ok(instr, KeyTag)) {
        kbuf = (int *) allocate((size_t)(*nbptr) * sizeof(int));
	get_data(instr, KeyTag, IntType, kbuf, *nbptr, 0);
	for (bp = *btptr, kp = kbuf; bp < *btptr + *nbptr; bp++)
	    Key(bp) = *kp++;
	free(kbuf);
	*ifptr |= KeyBit;
    }
#endif
}

#endif

/*
 * GET_SNAP_DENS: worker routine to input dens values.
 */

#ifndef get_snap_dens

#define get_snap_dens  _get_snap_dens

local void
_get_snap_dens(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
#ifdef Dens
    real *abuf, *ap;
    Body *bp;

    if (get_tag_ok(instr, DensityTag)) {
        abuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	get_data_coerced(instr, DensityTag, RealType, abuf, *nbptr, 0);
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++)
	    Dens(bp) = *ap++;
	free(abuf);
	*ifptr |= DensBit;
    }
#endif
}

#endif


/*
 * GET_SNAP_UINT: worker routine to input UDotInt values.
 */

#ifndef get_snap_uint

#define get_snap_uint  _get_snap_uint

local void
_get_snap_uint(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
#ifdef Uint
    real *abuf, *ap;
    Body *bp;

    if (get_tag_ok(instr, UdotIntTag)) {
        abuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	get_data_coerced(instr, UdotIntTag, RealType, abuf, *nbptr, 0);
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++)
	    Dens(bp) = *ap++;
	free(abuf);
	*ifptr |= UdotIntBit;
    }
#endif
}

#endif

/*
 * GET_SNAP_EPS: worker routine to input eps values.
 */

#ifndef get_snap_eps

#define get_snap_eps  _get_snap_eps

local void
_get_snap_eps(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
#ifdef Eps
    real *abuf, *ap;
    Body *bp;

    if (get_tag_ok(instr, EpsTag)) {
        abuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	get_data_coerced(instr, EpsTag, RealType, abuf, *nbptr, 0);
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++)
	    Eps(bp) = *ap++;
	free(abuf);
	*ifptr |= EpsBit;
    }
#endif
}

#endif


/*
 * GET_SNAP_PARTICLES: managing routine for input of particle data.
 */

#ifndef get_snap_particles

#define get_snap_particles  _get_snap_particles

local void
_get_snap_particles(instr, btptr, nbptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ifptr;			/* pointer to input bit flags */
{
    if (get_tag_ok(instr, ParticlesTag)) {
	if (*btptr == NULL) {
  	  *btptr = (Body *) allocate((size_t)(*nbptr) * sizeof(Body));
	}
	get_set(instr, ParticlesTag);
	get_snap_csys(instr, ifptr);
	get_snap_mass(instr, btptr, nbptr, ifptr);
	get_snap_phase(instr, btptr, nbptr, ifptr);
	get_snap_phi(instr, btptr, nbptr, ifptr);
	get_snap_acc(instr, btptr, nbptr, ifptr);
	get_snap_aux(instr, btptr, nbptr, ifptr);
	get_snap_key(instr, btptr, nbptr, ifptr);
	get_snap_dens(instr, btptr, nbptr, ifptr);
	get_snap_eps(instr, btptr, nbptr, ifptr);
	get_tes(instr, ParticlesTag);
    }
}

#endif

/*
 * GET_SNAP_DIAGNOSTICS: drone routine does nothing.
 */

#ifndef get_snap_diagnostics

#define get_snap_diagnostics  _get_snap_diagnostics

local void
_get_snap_diagnostics(instr, ifptr)
stream instr;			/* input stream, of course */
int *ifptr;			/* pointer to input bit flags */
{
}

#endif

/*
 * GET_SNAP: control routine for snapshot input.
 */

#ifndef get_snap

#define get_snap  _get_snap

local int
_get_snap(instr, btptr, nbptr, tsptr, ifptr)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
real *tsptr;			/* pointer to time of input */
int *ifptr;			/* pointer to input bit flags */
{
    *ifptr = 0;
    if (get_tag_ok(instr, SnapShotTag)) {
	get_set(instr, SnapShotTag);
	get_snap_parameters(instr, btptr, nbptr, tsptr, ifptr);
	get_snap_particles(instr, btptr, nbptr, ifptr);
	get_snap_diagnostics(instr, ifptr);
	get_tes(instr, SnapShotTag);
	return 1;
    } else
	return 0;
}

#endif

/*
 * GET_SNAP_BY_T: control routine for snapshot input, selected by time.
 */

#ifndef get_snap_by_t

#define get_snap_by_t  _get_snap_by_t

#ifndef TimeFuzz
#define TimeFuzz  0.001			/* slop allowed in time comparison  */
#endif

local int
_get_snap_by_t(instr, btptr, nbptr, tsptr, ifptr, times)
stream instr;			/* input stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
real *tsptr;			/* pointer to time of input */
int *ifptr;			/* pointer to input bit flags */
string times;
{
    *ifptr = 0;
    if (get_tag_ok(instr, SnapShotTag)) {
	get_set(instr, SnapShotTag);
	get_snap_parameters(instr, btptr, nbptr, tsptr, ifptr);
	if (streq(times, "all") ||
	      (*ifptr & TimeBit && within(*tsptr, times, TimeFuzz))) {
	    get_snap_particles(instr, btptr, nbptr, ifptr);
	    get_snap_diagnostics(instr, ifptr);
	}
	get_tes(instr, SnapShotTag);
        return 1;
    }
        return 0;
}

#endif

/*
 * GET_SNAP_NBODY: get the number of bodies in a snapshot
 */

#ifndef get_snap_nbody

#define get_snap_nbody  _get_snap_nbody

local int
_get_snap_nbody(instr)
stream instr;
{
  Body *btab = NULL;
  int nbody = 0;
  real tsnap;
  real bits;
  
  get_history(instr);
  if (!get_tag_ok(instr, SnapShotTag))
    return -1;
  get_snap(instr, &btab, &nbody, &tsnap, &bits);
  dprintf(0,"get_nbody: %d %f %d\n",nbody,tsnap,bits);
  free(btab);
  return nbody;
}

#endif
