/*
 * QUADINTER.C: quadrupole-moment force calculation from tabulated field.
 * Defines: quadinter().
 * Requires: Body, Pos(), Acc(), Phi(), quadfield, qfld, ...
 *
 *	16-mar-90	PJT	made GCC happy
 *      25-dec-02       pjt     bit more cleanup
 */

#include "quaddefs.h"

local inter_field(real r);
local force_eval(Body *b, real eps1, real eps2);

/*
 * QUADINTER: interpolate field and evaluate force on array of particles.
 */

quadinter(Body *btab, int nb, real eps1, real eps2)
{
    Body *b;

    for (b = btab; b < btab+nb; b++) {
	inter_field(absv(Pos(b)));
	force_eval(b, eps1, eps2);
    }
}

/*
 * Field moments used for force calculation.
 */

local real   Q00, P00;
local real   Q10, P10;
local vector Q11, P11;
local matrix Q22, P22;

#define INTERS(x,f,y,z)		/* INTERpolate Scalar */		\
{									\
    (x) = (f) * (y) + (1 - (f)) * (z);					\
}

#define INTERV(v,f,u,w)		/* INTERpolate Vector */		\
{									\
    register int _i;							\
    for (_i = 0; _i < NDIM; _i++)					\
	(v)[_i] = (f) * (u)[_i] + (1 - (f)) * (w)[_i];			\
}

#define INTERM(v,f,u,w)		/* INTERpolate Matrix */		\
{									\
    register int _i, _j;						\
    for (_i = 0; _i < NDIM; _i++)					\
	for (_j = 0; _j < NDIM; _j++)					\
	    (v)[_i][_j] = (f) * (u)[_i][_j] + (1 - (f)) * (w)[_i][_j];	\
}

/*
 * INTER_FIELD: interpolate field tables to given radius.
 */

local inter_field(real r)
{
    int i, j, k;
    real f;
    vector tmpv;
    matrix tmpm;

    i = 0;
    j = qfld.nqtab - 1;
    while (j - i > 1) {
	k = i + (j - i) / 2;
	if (qfld.radtab[k] <= r)
	    i = k;
	else
	    j = k;
    }
    if (! (j - i == 1))
	error("inter_field: loop ended with i = %d, j = %d", i, j);
    f = (qfld.radtab[j] - r) / (qfld.radtab[j] - qfld.radtab[i]);
    if (! (0.0 < f && f <= 1.0))
	error("inter_field: f = %f, i = %d, j = %d", f, i, j);
    INTERS(Q00, f, qfld.Q00tab[i], qfld.Q00tab[j]);
    INTERS(Q10, f, qfld.Q10tab[i], qfld.Q10tab[j]);
    INTERV(Q11, f, qfld.Q11tab[i], qfld.Q11tab[j]);
    INTERM(Q22, f, qfld.Q22tab[i], qfld.Q22tab[j]);
    INTERS(P00, f, qfld.P00tab[i], qfld.P00tab[j]);
    INTERS(P10, f, qfld.P10tab[i], qfld.P10tab[j]);
    INTERV(P11, f, qfld.P11tab[i], qfld.P11tab[j]);
    INTERM(P22, f, qfld.P22tab[i], qfld.P22tab[j]);
}

/*
 * FORCE_EVAL: compute force and potential on particle.
 */

local force_eval(Body *b, real eps1, real eps2)
{
    real rsq, r1i, r2i, r2is, r2iq, q11r, rq22r, tmp;
    vector q22r, p22r, tmpv;
    matrix tmpm;

    rsq = dotvp(Pos(b), Pos(b));
    r1i = 1.0 / sqrt(rsq + eps1*eps1);		/* invert softened radii    */
    r2i = 1.0 / sqrt(rsq + eps2*eps2);
    r2is = r2i * r2i;				/* compute various powers   */
    r2iq = r2i * r2is;
    DOTVP(q11r, Q11, Pos(b));			/* dot Pos into Q11, Q22    */
    MULMV(q22r, Q22, Pos(b));
    DOTVP(rq22r, Pos(b), q22r);
    tmp = Q00 * r1i*r1i*r1i +
	(3.0*q11r - 1.5*Q10 + 7.5*rq22r*r2is) * r2is*r2iq;
    MULVS(tmpv, Pos(b), tmp);
    SUBV(Acc(b), Acc(b), tmpv);			/* add radial acceleration  */
    MULVS(tmpv, q22r, 3.0*r2is);
    ADDV(tmpv, Q11, tmpv);
    MULVS(tmpv, tmpv, r2iq);
    ADDV(Acc(b), Acc(b), tmpv);			/* add non-radial accel.    */
    Phi(b) -= Q00 * r1i + (q11r - 0.5*Q10 + 1.5*rq22r * r2is) * r2iq; 
    MULVS(tmpv, Pos(b), P10);
    SUBV(Acc(b), Acc(b), tmpv);			/* add radial acceleration  */
    MULMV(p22r, P22, Pos(b));
    MULVS(tmpv, p22r, 3.0);
    ADDV(tmpv, tmpv, P11);
    ADDV(Acc(b), Acc(b), tmpv);			/* add non-radial accel.    */
    Phi(b) -=
	P00 + dotvp(Pos(b), P11) - 0.5*dotvp(Pos(b), Pos(b)) * P10 +
	    1.5*dotvp(Pos(b), p22r);
}
