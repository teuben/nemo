/*
 * PUT_SNAP.C: generic method for snapshot output.
 * Joshua E. Barnes  3-June-87  Princeton NJ.
 *	10 nov-88  put_snap_key/aux implemented	PJT
 *	18-nov-91  malloc -> allocate		PJT
 *	 7-mar-92  happy gcc2.0			pjt
 *	 2-may-92  fixed allocate() declaration PJT
 *	22-feb-94  ansi header (w/ allocate)	PJT
 *	25-mar-97  removed nested decl fflush() PJT
 *       8-oct-01  add some dens/eps if present PJT (for yanc)
 *      29-sep-05  fix for gcc4 supplying default prototypes
 *                 ** only potcode needed this,but we clearly need better solution for this **
 *      30-may-07  allocate() needs size_t argument casting for > 44.7M particles
 */

/*
 * put_snap is a generic method for writing snapshot data to a file,
 * to be #included in the compilation of an application program:
 *
 *	#define Body      ...
 *	#define Mass(b)   ...
 *	#define Phase(b)  ...
 *
 *	#include <snapshot/snapshot.h>
 *	#include <snapshot/put_snap.c>
 *
 * The simplest way to use put_snap is probably as follows:
 *
 *	stream outstr;
 *	Body *btab;
 *	int nbody, bits;
 *	real tsnap;
 *
 *	bits = <required_bits>;
 *	put_snap(outstr, &btab, &nbody, &tsnap, &bits);
 *
 * Notes: (1) put_snap() decides what data to output by looking at its
 * last argument, which is a word with a bit set for each component,
 * defined in snapshot/snapshot.h. (2) The name put_snap is a macro,
 * defined by put_snap.c, so the include statement must come
 * before the first usage. (3) The vanilla put_snap or any subsidiary
 * routine may be replaced by giving the macro name a definition before
 * including put_snap.c.  For example, to supply your own definition of
 * the mass-output routine, do the following:
 *
 *	#define put_mass  my_put_mass
 *
 *	#include <snapshot/put_snap.c>
 *
 *	local my_put_mass(....) {....}
 *
 * Look at the definition of the standard function(s) you are replacing to
 * find out what arguments to expect.
 */

/*
 * PUT_SNAP_PARAM: worker routine to output snapshot parameters.
 */

#ifndef put_snap_param

#define put_snap_param  _put_snap_param

local void
_put_snap_param(outstr, btptr, nbptr, tsptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
real *tsptr;			/* pointer to time of output */
int *ofptr;			/* pointer to output bit flags */
{
    put_set(outstr, ParametersTag);
    put_data(outstr, NobjTag, IntType, nbptr, 0);
    if (*ofptr & TimeBit)
	put_data(outstr, TimeTag, RealType, tsptr, 0);
    put_tes(outstr, ParametersTag);
}

#endif

/*
 * PUT_SNAP_CSYS: worker routine to output coordinate system.
 */

#ifndef put_snap_csys

#define put_snap_csys  _put_snap_csys

local void
_put_snap_csys(outstr, ofptr)
stream outstr;			/* output stream, of course */
int *ofptr;			/* pointer to output bit flags */
{
    int cs = CSCode(Cartesian, NDIM, 2);

    put_data(outstr, CoordSystemTag, IntType, &cs, 0);
}

#endif

/*
 * PUT_SNAP_MASS: worker routine to output mass data.
 */

#ifndef put_snap_mass

#define put_snap_mass  _put_snap_mass

local void
_put_snap_mass(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    real *mbuf, *mp;
    Body *bp;

    if (*ofptr & MassBit) {
#ifdef Mass
        mbuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	for (bp = *btptr, mp = mbuf; bp < *btptr + *nbptr; bp++)
	    *mp++ = Mass(bp);;
	put_data(outstr, MassTag, RealType, mbuf, *nbptr, 0);
	free(mbuf);
#else
	error("put_snap_mass: Mass undefined");
#endif
    }
}

#endif

/*
 * PUT_SNAP_PHASE: worker routine to output phasespace data.
 */

#ifndef put_snap_phase

#define put_snap_phase  _put_snap_phase

local void
_put_snap_phase(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    real *rvbuf, *rvp;
    Body *bp;

    if (*ofptr & PhaseSpaceBit) {
#ifdef Phase
        rvbuf = (real *) allocate((size_t)(*nbptr) * 2 * NDIM * sizeof(real));
	for (bp = *btptr, rvp = rvbuf; bp < *btptr + *nbptr; bp++) {
	    SETV(rvp, Phase(bp)[0]);
	    rvp += NDIM;
	    SETV(rvp, Phase(bp)[1]);
	    rvp += NDIM;
	}
	put_data(outstr, PhaseSpaceTag, RealType, rvbuf, *nbptr, 2, NDIM, 0);
	free(rvbuf);
#else
	error("put_snap_phase: Phase undefined");
#endif
    }
}

#endif

/*
 * PUT_SNAP_PHI: worker routine to output potential.
 */

#ifndef put_snap_phi

#define put_snap_phi  _put_snap_phi

local void
_put_snap_phi(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    real *pbuf, *pp;
    Body *bp;

    if (*ofptr & PotentialBit) {
#ifdef Phi
        pbuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	for (bp = *btptr, pp = pbuf; bp < *btptr + *nbptr; bp++)
	    *pp++ = Phi(bp);
	put_data(outstr, PotentialTag, RealType, pbuf, *nbptr, 0);
	free(pbuf);
#else
	error("put_snap_phi: Potential undefined");
#endif
    }
}

#endif

/*
 * PUT_SNAP_ACC: worker routine to output acceleration.
 */

#ifndef put_snap_acc

#define put_snap_acc  _put_snap_acc

local void
_put_snap_acc(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    real *abuf, *ap;
    Body *bp;

    if (*ofptr & AccelerationBit) {
#ifdef Acc
        abuf = (real *) allocate((size_t)(*nbptr) * NDIM * sizeof(real));
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++) {
	    SETV(ap, Acc(bp));
	    ap += NDIM;
	}
	put_data(outstr, AccelerationTag, RealType, abuf, *nbptr, NDIM, 0);
	free(abuf);
#else
	error("put_snap_acc: Acceleration undefined");
#endif
    }
}

#endif

/*
 * PUT_SNAP_AUX: worker routine to output aux data.
 */

#ifndef put_snap_aux

#define put_snap_aux  _put_snap_aux

local void
_put_snap_aux(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    real *abuf, *ap;
    Body *bp;
    
    if (*ofptr & AuxBit) {
#ifdef Aux
        abuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++)
	    *ap++ = Aux(bp);;
	put_data(outstr, AuxTag, RealType, abuf, *nbptr, 0);
	free(abuf);
#else
	error("put_snap_aux: Aux undefined");
#endif
    }
}

#endif

/*
 * PUT_SNAP_KEY: worker routine to output key data.
 */

#ifndef put_snap_key

#define put_snap_key  _put_snap_key

local void
_put_snap_key(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    int  *kbuf, *kp;
    Body *bp;

    if (*ofptr & KeyBit) {
#ifdef Key
        kbuf = (int *) allocate((size_t)(*nbptr) * sizeof(int));
	for (bp = *btptr, kp = kbuf; bp < *btptr + *nbptr; bp++)
	    *kp++ = Key(bp);;
	put_data(outstr, KeyTag, IntType, kbuf, *nbptr, 0);
	free(kbuf);
#else
	error("put_snap_key: Key undefined");
#endif
    }
}

#endif


/*
 * PUT_SNAP_DENS: worker routine to output dens data.
 */

#ifndef put_snap_dens

#define put_snap_dens  _put_snap_dens

local void
_put_snap_dens(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    real *abuf, *ap;
    Body *bp;
    
    if (*ofptr & DensBit) {
#ifdef Dens
        abuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++)
	    *ap++ = Dens(bp);;
	put_data(outstr, DensityTag, RealType, abuf, *nbptr, 0);
	free(abuf);
#else
	error("put_snap_dens: Dens undefined");
#endif
    }
}

#endif


/*
 * PUT_SNAP_EPS: worker routine to output eps data.
 */

#ifndef put_snap_eps

#define put_snap_eps  _put_snap_eps

local void
_put_snap_eps(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    real *abuf, *ap;
    Body *bp;
    
    if (*ofptr & EpsBit) {
#ifdef Eps
        abuf = (real *) allocate((size_t)(*nbptr) * sizeof(real));
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++)
	    *ap++ = Eps(bp);;
	put_data(outstr, EpsTag, RealType, abuf, *nbptr, 0);
	free(abuf);
#else
	error("put_snap_eps: Eps undefined");
#endif
    }
}

#endif

/*
 * PUT_SNAP_BODY: managing routine for output of particle data.
 */

#ifndef put_snap_body

#define put_snap_body  _put_snap_body

local void
_put_snap_body(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    if (*ofptr & (MassBit | PhaseSpaceBit | PotentialBit |
		    AccelerationBit | AuxBit | KeyBit)) {
	put_set(outstr, ParticlesTag);
	put_snap_csys(outstr, ofptr);
	put_snap_mass(outstr, btptr, nbptr, ofptr);
	put_snap_phase(outstr, btptr, nbptr, ofptr);
	put_snap_phi(outstr, btptr, nbptr, ofptr);
	put_snap_acc(outstr, btptr, nbptr, ofptr);
	put_snap_aux(outstr, btptr, nbptr, ofptr);
	put_snap_key(outstr, btptr, nbptr, ofptr);
	put_snap_dens(outstr, btptr, nbptr, ofptr);
	put_snap_eps(outstr, btptr, nbptr, ofptr);
	put_tes(outstr, ParticlesTag);
    }
}

#endif

/*
 * PUT_SNAP_DIAGNOSTICS: drone routine does nothing.
 */


#ifndef put_snap_diagnostics

#define put_snap_diagnostics  _put_snap_diagnostics

local void 
_put_snap_diagnostics(
stream outstr,			/* output stream, of course */
int *ofptr)			/* pointer to output bit flags */
{
    if (*ofptr & (EnergyBit | KETensorBit | PETensorBit |
		    AMTensorBit | CMPhaseSpaceBit))
	error("put_snap_diagnostics: not implemented");
}
#else
local void put_snap_diagnostics(stream outstr, int *ofptr); 
#endif

/*
 * PUT_SNAP: controling routine for snapshot output.
 */

#ifndef put_snap

#define put_snap  _put_snap

local void
_put_snap(outstr, btptr, nbptr, tsptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
real *tsptr;			/* pointer to time of output */
int *ofptr;			/* pointer to output bit flags */
{
    if (ofptr) {
	put_set(outstr, SnapShotTag);
	put_snap_param(outstr, btptr, nbptr, tsptr, ofptr);
	put_snap_body(outstr, btptr, nbptr, ofptr);
	put_snap_diagnostics(outstr, ofptr);
	put_tes(outstr, SnapShotTag);
	fflush(outstr);
    }
}

#endif
