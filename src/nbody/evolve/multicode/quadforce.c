/*
 * QUADFORCE.C: quadrupole-moment force calculation routine.
 * Defines: quadforce().
 * Requires: Body, Mass(), Pos(), Acc(), Phi(), quadfield, qfld, ...
 *
 *	16-mar-90	PJT	made GCC happy
 *      20-may-94       pjt     remove allocate() decl. into headers
 *      25-dec-02       pjt     better ANSI C
 */

#include "quaddefs.h"

/* SHADOW: local per-body data. */

typedef struct {
    Body *link;					/* pointer to actual body   */
    real rad0, rads1, rads2;			/* exact, softened radii    */
} shadow, *shadowptr;

local void init_quad_field(shadow *, int, real);
local void int_quad_force(shadow *, int);
local void ext_quad_force(shadow *, int);
local int rankrad(shadowptr, shadowptr);

/*
 * QUADFORCE: compute the force-field of a spheroidal distribution.
 */

quadforce(Body *btab, int nb, real eps1, real eps2)
{
    shadowptr shad;
    Body *b;
    int i;
    real rsq;

    shad = (shadowptr) allocate(nb * sizeof(shadow));
    for (b = btab, i = 0; i < nb; b++, i++) {	/* init shadowing array     */
 	shad[i].link = b;			/*   make link to body      */
	rsq = dotvp(Pos(b), Pos(b));		/*   find radius squared    */
	shad[i].rad0 = sqrt(rsq);		/*   set exact radius       */
	shad[i].rads1 = sqrt(rsq + eps1*eps1);	/*   set softened radii     */
	shad[i].rads2 = sqrt(rsq + eps2*eps2);
    }
    qsort(shad, nb, sizeof(shadow), rankrad);	/* sort shadows by radii    */
    init_quad_field(shad, nb, eps1);		/* prepare field tables     */
    int_quad_force(shad, nb);			/* compute internal force   */
    ext_quad_force(shad, nb);			/* and external force       */
    free(shad);
}

/*
 * RANKRAD: ranking function for quicksort.
 */

local int rankrad(shadowptr sp1, shadowptr sp2)
{
    return (sp1->rad0 < sp2->rad0 ? -1 : 1);	/* compare body radii       */
}

/*
 * INIT_QUAD_FIELD: initialize quadrupole-field tables.  If qfld.nqtab
 * is non-zero, assume radii are pretabulated; otherwise use sampled
 * particle radii.
 */

local int jint, jext;
local real rfield;

local void init_quad_field(shadow *shad, int nb, real eps1)
{
    int ktab, j;

    if (qfld.nqtab == 0) {			/* field points not set?    */
	ktab = nb / MIN(MQTAB-2, nb);		/*   set stride thru bodies */
	while (ktab * MIN(MQTAB-2, nb) < nb)	/*   assure nqtab <= MQTAB  */
	    ktab++;
	qfld.radtab[0] = 0.0;			/*   set zeroth radius      */
	for (j = 1; ktab*j <= nb; j++)		/*   loop setting radii     */
	    qfld.radtab[j] = shad[ktab*j - 1].rad0;
 	qfld.radtab[j] = 2.0 * shad[nb-1].rad0;	/*   set final radius       */
	qfld.nqtab = j + 1;			/*   count radii tabulated  */
    }
    jint = 0;
    jext = qfld.nqtab - 1;
    if (qfld.radtab[jext] < shad[nb-1].rad0)	/* check extent of radtab   */
	error("init_quad_field: r < rmax");
}

/*
 * SAVE_INT_FIELD, SAVE_EXT_FIELD: macros for storing field values.
 */

#define SAVE_INT_FIELD() {					\
    qfld.Q00tab[jint] = Q00;					\
    qfld.Q10tab[jint] = Q10;					\
    SETV(qfld.Q11tab[jint], Q11);				\
    SETM(qfld.Q22tab[jint], Q22);				\
    jint++;							\
    if (jint < qfld.nqtab)					\
	rfield = qfld.radtab[jint];				\
}

#define SAVE_EXT_FIELD() {					\
    qfld.P00tab[jext] = P00;					\
    qfld.P10tab[jext] = P10;					\
    SETV(qfld.P11tab[jext], P11);				\
    SETM(qfld.P22tab[jext], P22);				\
    jext--;							\
    if (jext >= 0)						\
	rfield = qfld.radtab[jext];				\
}

/*
 * INT_QUAD_FORCE: compute force due to interior particles.
 */

local void int_quad_force(shadow *shad, int nb)
{
    real Q00, Q10, r1i, r2i, r2is, r2iq, q11r, rq22r, tmp;
    vector Q11, q22r, tmpv;
    matrix Q22, tmpm;
    int i;
    Body *b;

    Q00 = 0.0;					/* initialize moments       */
    Q10 = 0.0;
    CLRV(Q11);
    CLRM(Q22);
    SAVE_INT_FIELD();				/* tabulate field at r = 0  */
    for (i = 0; i < nb; i++) {			/* loop from inside out     */
	if (rfield <= shad[i].rad0)		/*   internal field point?  */
	    SAVE_INT_FIELD();
	b = shad[i].link;			/*   get i'th from center   */
	r1i = 1.0 / shad[i].rads1;		/*   access stored radii    */
	r2i = 1.0 / shad[i].rads2;
	r2is = r2i * r2i;			/*   compute various powers */
	r2iq = r2i * r2is;
	DOTVP(q11r, Q11, Pos(b));		/*   dot Pos into Q11, Q22  */
	MULMV(q22r, Q22, Pos(b));
	DOTVP(rq22r, Pos(b), q22r);
	tmp = Q00 * r1i*r1i*r1i +
	    (3.0*q11r - 1.5*Q10 + 7.5*rq22r*r2is) * r2is*r2iq;
	MULVS(tmpv, Pos(b), tmp);
	SUBV(Acc(b), Acc(b), tmpv);		/*   add radial accel.      */
	MULVS(tmpv, q22r, 3.0*r2is);
	ADDV(tmpv, Q11, tmpv);
	MULVS(tmpv, tmpv, r2iq);
	ADDV(Acc(b), Acc(b), tmpv);		/*   add non-radial accel.  */
	Phi(b) -= Q00 * r1i + (q11r - 0.5*Q10 + 1.5*rq22r * r2is) * r2iq; 
	Q00 += Mass(b);				/*   sum monopole moment    */
	MULVS(tmpv, Pos(b), Mass(b));
	ADDV(Q11, Q11, tmpv);			/*   sum dipole moment      */
	OUTVP(tmpm, tmpv, Pos(b));
	ADDM(Q22, Q22, tmpm);			/*   sum quadpole moment    */
	TRACEM(Q10, Q22);
    }
    SAVE_INT_FIELD();				/* tabulate field at inf    */
    if (jint != qfld.nqtab)			/* check field tabulation   */
	error("int_quad_force: jint = %d != %d", jint, qfld.nqtab);
}

/*
 * EXT_QUAD_FORCE: compute force due to exterior particles.
 */

local void ext_quad_force(shadow *shad, int nb)
{
    real P00, P10, r2i, r2is, tmp;
    vector P11, p22r, tmpv;
    matrix P22, tmpm;
    int i;
    Body *b;

    P00 = 0.0;					/* initialize moments       */
    P10 = 0.0;
    CLRV(P11);
    CLRM(P22);
    SAVE_EXT_FIELD();				/* tabulate field at inf.   */
    for (i = nb - 1; i >= 0; i--) {		/* start from outside       */
	if (shad[i].rad0 <= rfield)		/*   save field values?     */
	    SAVE_EXT_FIELD();
	b = shad[i].link;			/*   get i'th from center   */
	MULVS(tmpv, Pos(b), P10);
	SUBV(Acc(b), Acc(b), tmpv);		/*   add radial accel.      */
	MULMV(p22r, P22, Pos(b));
	MULVS(tmpv, p22r, 3.0);
	ADDV(tmpv, tmpv, P11);
	ADDV(Acc(b), Acc(b), tmpv);		/*   add non-radial accel.  */
	Phi(b) -=
	    P00 + dotvp(Pos(b), P11) - 0.5*dotvp(Pos(b), Pos(b)) * P10 +
		1.5*dotvp(Pos(b), p22r);
	P00 += Mass(b) / shad[i].rads1;		/*   sum monopole moment    */
	r2i = 1.0 / shad[i].rads2;
	r2is = r2i * r2i;
	tmp = Mass(b) * r2i * r2is;
	P10 += tmp;
	MULVS(tmpv, Pos(b), tmp);
	ADDV(P11, P11, tmpv);			/*   sum dipole moment      */
	MULVS(tmpv, tmpv, r2is);
	OUTVP(tmpm, tmpv, Pos(b));
	ADDM(P22, P22, tmpm);			/*   sum quadpole moment    */
    }
    SAVE_EXT_FIELD();				/* tabulate field at r = 0  */
    if (jext != -1)
	error("ext_quad_force: jext = %d != %d", jext, -1);
}
