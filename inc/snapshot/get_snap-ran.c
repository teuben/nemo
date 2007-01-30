// warning: this is still unreliable
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
 *	 8-jan-98 converted to use random + common I/O		pjt
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
static  int first_io_get = 1;

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
	get_data(instr, NobjTag, IntType, &nbody, 0);
	if (*btptr != NULL && nbody > *nbptr)	/* bigger than expected? */
	    error("get_snap_parameters: %s = %d is too big\n",
		  NobjTag, nbody);
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
	printf("get_snap_csys: assuming %s = %#o\n",
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
    int len, offset, ntodo, nbucket;

    if (get_tag_ok(instr, MassTag)) {
        mbuf = (real *) open_common(0);
        nbucket = get_common(0, sizeof(real), 1);

        get_data_set(instr, MassTag, RealType, *nbptr, 0);
        offset = 0;
        ntodo = *nbptr;
        len = MIN(ntodo,nbucket);
        bp = *btptr;
        do {			/* coerced ??? */
            get_data_ran(instr,MassTag,mbuf,offset,len);
            for (mp=mbuf; bp < *btptr + offset + len; bp++) {
                Mass(bp) = *mp++;
            }
            offset += len;
            len = MIN(ntodo-offset, nbucket);
        } while (len > 0);
        get_data_tes(instr, MassTag);
        close_common(0);

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
    int len, offset, ntodo, nbucket, bsize=2*NDIM;

    if (get_tag_ok(instr, PhaseSpaceTag)) {
        rvbuf = (real *) open_common(0);
        nbucket = get_common(0, sizeof(real), bsize);
        get_data_set(instr, PhaseSpaceTag, RealType, *nbptr, 2, NDIM, 0);
        offset = 0;
        ntodo = *nbptr;
        len = MIN(ntodo,nbucket);
        bp = *btptr;
        do {
            get_data_ran(instr,PhaseSpaceTag,rvbuf,offset*bsize,len*bsize);
	    for (rvp = rvbuf; bp < *btptr + offset + len; bp++) {
	        SETV(Phase(bp)[0], rvp);
	        rvp += NDIM;
	        SETV(Phase(bp)[1], rvp);
	        rvp += NDIM;
            }
            offset += len;
            len = MIN(ntodo-offset,nbucket);
	} while (len>0);
        get_data_tes(instr, PhaseSpaceTag);
	close_common(0);

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
    int len, offset, ntodo, nbucket;

    if (get_tag_ok(instr, PotentialTag)) {
        pbuf = (real *) open_common(0);
        nbucket = get_common(0, sizeof(real), 1);

        get_data_set(instr, PotentialTag, RealType, *nbptr, 0);
        offset = 0;
        ntodo = *nbptr;
        len = MIN(ntodo,nbucket);
        bp = *btptr;
        do {
            get_data_ran(instr,PotentialTag,pbuf,offset,len);
            for (pp=pbuf; bp < *btptr + offset + len; bp++)
                Phi(bp) = *pp++;
            offset += len;
            len = MIN(ntodo-offset, nbucket);
        } while (len > 0);
        get_data_tes(instr, PotentialTag);
        close_common(0);

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
    int len, offset, ntodo, nbucket;

    if (get_tag_ok(instr, AccelerationTag)) {
        abuf = (real *) open_common(0);
        nbucket = get_common(0, sizeof(real), NDIM);
        get_data_set(instr, AccelerationTag, RealType, *nbptr, NDIM, 0);
        offset = 0;
        ntodo = *nbptr;
        len = MIN(ntodo,nbucket);
        bp = *btptr;
        do {
            get_data_ran(instr,AccelerationTag,abuf,offset*NDIM,len*NDIM);
	    for (ap = abuf; bp < *btptr + offset + len; bp++) {
	        SETV(Acc(bp), ap);
	        ap += NDIM;
            }
            offset += len;
            len = MIN(ntodo-offset,nbucket);
	} while (len>0);
        get_data_tes(instr, AccelerationTag);
	close_common(0);

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
    int len, offset, ntodo, nbucket;

    if (get_tag_ok(instr, AuxTag)) {
        abuf = (real *) open_common(0);
        nbucket = get_common(0, sizeof(real), 1);

        get_data_set(instr, AuxTag, RealType, *nbptr, 0);
        offset = 0;
        ntodo = *nbptr;
        len = MIN(ntodo,nbucket);
        bp = *btptr;
        do {
            get_data_ran(instr,AuxTag,abuf,offset,len);
            for (ap=abuf; bp < *btptr + offset + len; bp++)
                Aux(bp) = *ap++;
            offset += len;
            len = MIN(ntodo-offset, nbucket);
        } while (len > 0);
        get_data_tes(instr, AuxTag);
        close_common(0);

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
    int len, offset, ntodo, nbucket;

    if (get_tag_ok(instr, KeyTag)) {
        kbuf = (int *) open_common(0);
        nbucket = get_common(0, sizeof(int), 1);

        get_data_set(instr, KeyTag, IntType, *nbptr, 0);
        offset = 0;
        ntodo = *nbptr;
        len = MIN(ntodo,nbucket);
        bp = *btptr;
        do {
            get_data_ran(instr,KeyTag,kbuf,offset,len);
            for (kp=kbuf; bp < *btptr + offset + len; bp++)
                Key(bp) = *kp++;
            offset += len;
            len = MIN(ntodo-offset, nbucket);
        } while (len > 0);
        get_data_tes(instr, KeyTag);
        close_common(0);

	*ifptr |= KeyBit;
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
	    *btptr = (Body *) allocate(*nbptr * sizeof(Body));
	}
	get_set(instr, ParticlesTag);
	get_snap_csys(instr, ifptr);
	get_snap_mass(instr, btptr, nbptr, ifptr);
	get_snap_phase(instr, btptr, nbptr, ifptr);
	get_snap_phi(instr, btptr, nbptr, ifptr);
	get_snap_acc(instr, btptr, nbptr, ifptr);
	get_snap_aux(instr, btptr, nbptr, ifptr);
	get_snap_key(instr, btptr, nbptr, ifptr);
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
#if 1
    if (first_io_get) {
        dprintf(0,"[get_snap: experimental %d buffered random I/O]\n",
			get_common(0,1,1));
        first_io_get = 0;
    }
#endif
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
#if 1
    if (first_io_get) {
        dprintf(0,"[get_snap_by_t: experimental %d buffered random I/O]\n",
			get_common(0,1,1));
        first_io_get = 0;
    }
#endif
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

