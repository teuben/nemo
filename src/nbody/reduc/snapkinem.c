/*
 * SNAPKINEM.C: compute various diagnostic properties for a (typically)
 * lagrangian set of particles with specified weights.
 * 
 *	5-mar-89	V1.2  -- JEB
 *	2-may-92	V1.3  helpvec/usage for new NEMO		    PJT
 *     11-jun-92        V2.0  repaired sign error (?) in angular momentum   PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {		/* DEFAULT INPUT PARAMETERS		    */
    "in=???\n			 input file name (snapshot)",
    "weight=1\n			 weighting for particles",
    "rcut=0.0\n			 cutoff radius effective if > 0.0",
    "times=all\n		 range of times to analyze",
    "VERSION=2.0\n		 11-jun-92 PJT",
    NULL,
};

string usage="compute various diagnostics with specified weights.";


Body *btab = NULL;		/* pointer to array of bodies		    */
int nbody;			/* number of bodies in array		    */
real tsnap = 0.0;		/* time associated with data		    */

rproc weight;			/* weighting function for bodies	    */

real rcut;			/* cutoff to suppress outlying bodies	    */

nemo_main()
{
    stream instr;
    string times;
    rproc btrtrans();
    int bits;

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    weight = btrtrans(getparam("weight"));
    rcut = getdparam("rcut");
    times = getparam("times");
    do {
	get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);
	if (bits & PhaseSpaceBit) {
	    snapkinem();
	    showkinem();
	}
    } while (bits != 0);
}

snapkinem()
{
    roughcenter();
    findcenter();
    findmoment();
}

int n_tot;			/* number of bodies with positive weight    */

vector cm_pos;			/* rough center of mass position	    */

int n_sum;			/* number of bodies contributing	    */

real w_sum;			/* sum of body weights			    */

vector w_pos;			/* weighted center of mass position	    */
vector w_vel;			/* weighted center of mass velocity	    */
vector w_jvec;			/* specific angular momentum		    */

matrix w_qpole;			/* weighted quadrupole moment		    */
matrix w_keten;			/* weighted kinetic energy tensor	    */

roughcenter()
{
    int i;
    Body *b;
    real w_tot, w_b;
    vector tmpv;

    n_tot = 0;
    w_tot = 0.0;
    CLRV(cm_pos);
    for (i = 0, b = btab; i < nbody; i++, b++) {
	w_b = (weight)(b, tsnap, i);
	if (w_b < 0.0)
	    warning("weight[%d] = %g < 0", i, w_b);
	if (w_b > 0.0) {
	    n_tot = n_tot + 1;
	    w_tot = w_tot + w_b;
	    MULVS(tmpv, Pos(b), w_b);
	    ADDV(cm_pos, cm_pos, tmpv);
	}
    }
    if (w_tot == 0.0)
	error("total weight is zero");
    DIVVS(cm_pos, cm_pos, w_tot);
}

findcenter()
{
    int i;
    Body *b;
    real w_b;
    vector tmpv;

    n_sum = 0;
    w_sum = 0.0;
    CLRV(w_pos);
    CLRV(w_vel);
    for (i = 0, b = btab; i < nbody; i++, b++) {
	w_b = (weight)(b, tsnap, i);
	if (w_b > 0.0 && (rcut <= 0.0 || distv(Pos(b), cm_pos) < rcut)) {
	    n_sum = n_sum + 1;
	    w_sum = w_sum + w_b;
	    MULVS(tmpv, Pos(b), w_b);
	    ADDV(w_pos, w_pos, tmpv);
	    MULVS(tmpv, Vel(b), w_b);
	    ADDV(w_vel, w_vel, tmpv);
	}
    }
    if (w_sum == 0.0)
	error("weight within rcut is zero");
    DIVVS(w_pos, w_pos, w_sum);
    DIVVS(w_vel, w_vel, w_sum);
}

findmoment()
{
    int i;
    Body *b;
    real w_b;
    vector tmpv, pos_b, vel_b;
    matrix tmpm;

    CLRV(w_jvec);
    CLRM(w_qpole);
    CLRM(w_keten);
    for (i = 0, b = btab; i < nbody; i++, b++) {
	w_b = (weight)(b, tsnap, i);
	if (w_b > 0.0 && (rcut <= 0.0 || distv(Pos(b), cm_pos) < rcut)) {
	    SUBV(pos_b, Pos(b), w_pos);
	    SUBV(vel_b, Vel(b), w_vel);
	    CROSSVP(tmpv, pos_b, vel_b);        /* jun-92: repaired sign */
	    MULVS(tmpv, tmpv, w_b);
	    ADDV(w_jvec, w_jvec, tmpv);
	    MULVS(tmpv, pos_b, w_b);
	    OUTVP(tmpm, tmpv, pos_b);
	    ADDM(w_qpole, w_qpole, tmpm);
	    MULVS(tmpv, vel_b, w_b);
	    OUTVP(tmpm, tmpv, vel_b);
	    ADDM(w_keten, w_keten, tmpm);
	}
    }
    DIVVS(w_jvec, w_jvec, w_sum);
    DIVMS(w_qpole, w_qpole, w_sum);
    DIVMS(w_keten, w_keten, w_sum);
}

showkinem()
{
    printf("\n");
    printf("time:%7.3f    n_tot:%6d    n_sum:%6d    w_sum:%6g\n",
	   tsnap, n_tot, n_sum, w_sum);
    printf("\n");
    printf("%12s  %10s  %10s  %10s  %10s\n",
	   "            ", "length", "x", "y", "z");
    printvec("pos:", w_pos);
    printvec("vel:", w_vel);
    printvec("jvec:", w_jvec);
    printf("\n");
    printf("%12s  %10s  %10s  %10s  %10s\n",
	   "            ", "eigval", "x", "y", "z");
    printeig("qpole:", w_qpole);
    printf("\n");
    printeig("keten:", w_keten);
}

printvec(name, vec)
string name;
vector vec;
{
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n",
	   name, absv(vec), vec[0], vec[1], vec[2]);
}

#include "nrutil.h"

printeig(name, mat)
string name;
matrix mat;
{
    float **q, *d, **v;
    int i, j, nrot;

    q = fmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++)
	for (j = 1; j <= 3; j++)
	    q[i][j] = mat[i-1][j-1];
    d = fvector(1, 3);
    v = fmatrix(1, 3, 1, 3);
    jacobi(q, 3, d, v, &nrot);
    eigsrt(d, v, 3);
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n", name,
	   d[1], v[1][1], v[2][1], v[3][1]);
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n", "            ",
	   d[2], v[1][2], v[2][2], v[3][2]);
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n", "            ",
	   d[3], v[1][3], v[2][3], v[3][3]);
}
