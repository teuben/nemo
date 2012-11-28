/*
 * PUT_SNAP.C: generic method for snapshot output.
 * Joshua E. Barnes  3-June-87  Princeton NJ.
 *	10-nov-88  put_snap_key/aux implemented	PJT
 *	12-sep-90  mybody's story; if all masses the same - header write PJT
 *       7-mar-92  Happy gcc2.0                                          pjt
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
 * PUT_SNAP_STORY: worker routine to output snapshot stories.
 */
#ifndef put_snap_story

#define put_snap_story  _put_snap_story

local void
_put_snap_story(outstr, btptr, nbptr, ofptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
{
    char *allocate(), sname[32];
    void free();
    int  *ip, i, nstories;
    Body *bp;

#ifdef Story
	/* first scan the Body, to see how many stories were there */
    for (nstories=0, bp=*btptr; bp < *btptr + *nbptr; bp++)
        if (Story(bp) != NULL) nstories++;
    dprintf(1,"Found %d stories to write\n",nstories);
    if (nstories==0) return;
    ip = (int *) allocate(nstories * sizeof(int));
    for (nstories=0, bp=*btptr, i=0; bp < *btptr + *nbptr; bp++, i++)
        if (Story(bp) != NULL)
            ip[nstories++] = i;
    put_set (outstr,"Story");
    put_data(outstr,"Nstories",IntType,&nstories,0);
    put_data(outstr,"IStories",IntType,ip,nstories,0);
    for (i=0; i<nstories; i++) {
        sprintf(sname,"story_%d",ip[i]);
        put_string(outstr,sname,Story(*btptr+ip[i]));
    }
    free(ip);
    put_tes(outstr,"Story");
#else
    error("put_snap_mass: Stories undefined\n");
#endif
}

#endif

/*
 * PUT_SNAP_PARAM: worker routine to output snapshot parameters.
 */

#ifndef put_snap_param

#define put_snap_param  _put_snap_param

local void
_put_snap_param(outstr, btptr, nbptr, tsptr, ofptr, umptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
real *tsptr;			/* pointer to time of output */
int *ofptr;			/* pointer to output bit flags */
int *umptr;			/* pointer to unique mass flag */
{
    Body *bp;
    real m0;
    int w;

    put_set(outstr, ParametersTag);
    put_data(outstr, NobjTag, IntType, nbptr, 0);
    if (*ofptr & TimeBit)
	put_data(outstr, TimeTag, RealType, tsptr, 0);
    for (bp=*btptr,w=1,m0=Mass(*btptr); bp < *btptr + *nbptr; bp++)
        if(m0 != Mass(bp)) {
            w = 0;
            break;
        }
    if (w) 
        put_data(outstr, MassTag, RealType, &m0, 0);
    *umptr = w;
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
    char *allocate();
    void free();
    real *mbuf, *mp;
    Body *bp;

    if (*ofptr & MassBit) {
#ifdef Mass
	mbuf = (real *) allocate(*nbptr * sizeof(real));
	for (bp = *btptr, mp = mbuf; bp < *btptr + *nbptr; bp++)
	    *mp++ = Mass(bp);;
	put_data(outstr, MassTag, RealType, mbuf, *nbptr, 0);
	free(mbuf);
#else
	error("put_snap_mass: Mass undefined\n");
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
    char *allocate();
    void free();
    real *rvbuf, *rvp;
    Body *bp;

    if (*ofptr & PhaseSpaceBit) {
#ifdef Phase
	rvbuf = (real *) allocate(*nbptr * 2 * NDIM * sizeof(real));
	for (bp = *btptr, rvp = rvbuf; bp < *btptr + *nbptr; bp++) {
	    SETV(rvp, Phase(bp)[0]);
	    rvp += NDIM;
	    SETV(rvp, Phase(bp)[1]);
	    rvp += NDIM;
	}
	put_data(outstr, PhaseSpaceTag, RealType, rvbuf, *nbptr, 2, NDIM, 0);
	free(rvbuf);
#else
	error("put_snap_phase: Phase undefined\n");
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
    char *allocate();
    void free();
    real *pbuf, *pp;
    Body *bp;

    if (*ofptr & PotentialBit) {
#ifdef Phi
	pbuf = (real *) allocate(*nbptr * sizeof(real));
	for (bp = *btptr, pp = pbuf; bp < *btptr + *nbptr; bp++)
	    *pp++ = Phi(bp);
	put_data(outstr, PotentialTag, RealType, pbuf, *nbptr, 0);
	free(pbuf);
#else
	error("put_snap_phi: Potential undefined\n");
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
    char *allocate();
    void free();
    real *abuf, *ap;
    Body *bp;

    if (*ofptr & AccelerationBit) {
#ifdef Acc
	abuf = (real *) allocate(*nbptr * NDIM * sizeof(real));
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++) {
	    SETV(ap, Acc(bp));
	    ap += NDIM;
	}
	put_data(outstr, AccelerationTag, RealType, abuf, *nbptr, NDIM, 0);
	free(abuf);
#else
	error("put_snap_acc: Acceleration undefined\n");
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
    char *allocate();
    void free();
    real *abuf, *ap;
    Body *bp;
    
    if (*ofptr & AuxBit) {
#ifdef Aux
	abuf = (real *) allocate(*nbptr * sizeof(real));
	for (bp = *btptr, ap = abuf; bp < *btptr + *nbptr; bp++)
	    *ap++ = Aux(bp);;
	put_data(outstr, AuxTag, RealType, abuf, *nbptr, 0);
	free(abuf);
#else
	error("put_snap_aux: Aux undefined\n");
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
    char *allocate();
    void free();
    int  *kbuf, *kp;
    Body *bp;

    if (*ofptr & KeyBit) {
#ifdef Key
	kbuf = (int *) allocate(*nbptr * sizeof(int));
	for (bp = *btptr, kp = kbuf; bp < *btptr + *nbptr; bp++)
	    *kp++ = Key(bp);;
	put_data(outstr, KeyTag, IntType, kbuf, *nbptr, 0);
	free(kbuf);
#else
	error("put_snap_key: Key undefined\n");
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
_put_snap_body(outstr, btptr, nbptr, ofptr, umptr)
stream outstr;			/* output stream, of course */
Body **btptr;			/* pointer to body array */
int *nbptr;			/* pointer to number of bodies */
int *ofptr;			/* pointer to output bit flags */
int *umptr;			/* pointer to unique mass flag */
{
    if (*ofptr & (MassBit | PhaseSpaceBit | PotentialBit |
		    AccelerationBit | AuxBit | KeyBit)) {
	put_set(outstr, ParticlesTag);
	put_snap_csys(outstr, ofptr);
        if (! *umptr)
	    put_snap_mass(outstr, btptr, nbptr, ofptr);
	put_snap_phase(outstr, btptr, nbptr, ofptr);
	put_snap_phi(outstr, btptr, nbptr, ofptr);
	put_snap_acc(outstr, btptr, nbptr, ofptr);
	put_snap_aux(outstr, btptr, nbptr, ofptr);
	put_snap_key(outstr, btptr, nbptr, ofptr);
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
_put_snap_diagnostics(outstr, ofptr)
stream outstr;			/* output stream, of course */
int *ofptr;			/* pointer to output bit flags */
{
    if (*ofptr & (EnergyBit | KETensorBit | PETensorBit |
		    AMTensorBit | CMPhaseSpaceBit))
	error("put_snap_diagnostics: not implemented\n");
}

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
int *nbptr;                    /* pointer to number of bodies */
real *tsptr;			/* pointer to time of output */
int *ofptr;			/* pointer to output bit flags */
{
    int um;                     /* unique masses ? */
    int fflush();

    if (ofptr) {
	put_set(outstr, SnapShotTag);
	put_snap_param(outstr, btptr, nbptr, tsptr, ofptr, &um);
	put_snap_body(outstr, btptr, nbptr, ofptr, &um);
	put_snap_story(outstr, btptr, nbptr, ofptr);
	put_snap_diagnostics(outstr, ofptr);
	put_tes(outstr, SnapShotTag);
	fflush(outstr);		/* ensurce all written to disk */
    }
}

#endif

