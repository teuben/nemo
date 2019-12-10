/*
 * SNAPRECT.C: transform snapshot to coordinate system which diagonalizes
 * the moment-of-inertia tensor of a specified set of particles.
 *
 *  5-mar-89    V1.1  some changes      JEB
 *  9-dec-90    V1.2  helpvec           PJT
 *  5-dec-93	V1.3  transform using a test= snapshot (Roald 1yr)	PJT
 *  4-jul-94    V1.3b also output theta,phi of principle axes		pjt
 *			only now implemented test=
 * 22-aug-00        c bit more ansi cc					pjt
 * 15-mar-02        d added option log to allow piping of output        WD 
 * 31-dec-02    V1.4  gcc3/SINGLEPREC
 *  1-mar-2019
 *
 *  TODO:
 *    might need a snapshape program, see e.g.
 *	astro-ph/0505179 (Lee & Kang 2005)
 *      Reconstructing the Triaxial Shapes of Dark Matter Halos from the Anisotropic Spatial Distributions of their Substructures
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>
#include <bodytransc.h>

string defv[] = {
    "in=???\n			  input file name",
    "out=???\n			  output file name",
    "weight=1\n			  weighting for particles",
    "rcut=0.0\n			  cutoff radius effective if > 0.0",
    "times=all\n		  range of times to transform",
    "log=\n                       write to this file instead of stdout",
    "test=\n                      Use this snapshot instead, to transform",
    "pos=t\n                      Use pos, or vel, or rectify",
    "VERSION=1.5a\n		  6-nov-2019 PJT",
    NULL,
};

string usage = "rectify a snapshot using selected particles and weight";

Body *btab = NULL;		/* pointer to array of bodies		    */
int nbody;			/* number of bodies in array		    */
real tsnap = 0.0;		/* time associated with data		    */
Body *testtab = NULL;
bool Qtest;
bool Qpos; 

rproc_body weight;		/* weighting function for bodies	    */

real rcut;			/* cutoff to suppress outlying bodies	    */

void roughcenter(body *btab);
void findcenter(body *btab);
void findmoment(body *btab);
void findmoment_vel(body *btab);
void snaptransform(stream out);
void eigenframe(vector frame[], matrix mat);
void printvec(string name, vector vec, stream out);
void xyz2rtp(vector xyz, vector rtp);


void jacobi(float **a,int n,float d[],float **v,int *nrot);
void eigsrt(float d[],float **v, int n);

void nemo_main(void)
{
    stream instr, outstr, logstr;
    string times;
    int bits;

    weight = btrtrans(getparam("weight"));
    Qpos = getbparam("pos");
    Qtest = hasvalue("test");
    if (Qtest) {
    	instr = stropen(getparam("test"),"r");
	get_history(instr);
    	get_snap(instr, &testtab, &nbody, &tsnap, &bits);
    	strclose(instr);
    	roughcenter(testtab);
    	findcenter(testtab);
    	findmoment(testtab);
        free(testtab);
    }
    instr = stropen(getparam("in"), "r");
    get_history(instr);
    outstr = stropen(getparam("out"), "w");
    logstr = hasvalue("log")? stropen(getparam("log"), "w!") : stdout;
    put_history(outstr);
    rcut = getdparam("rcut");
    times = getparam("times");
    do {
	get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);
	if (bits & PhaseSpaceBit) {
            roughcenter(btab);
            findcenter(btab);
            if (!Qtest) findmoment(btab);
            snaptransform(logstr);
	    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
	}
    } while (bits != 0);
}


vector cm_pos;			/* rough center of mass position	    */

real w_sum;			/* sum of body weights			    */

vector w_pos;			/* weighted center of mass position	    */
vector w_vel;			/* weighted center of mass velocity	    */

matrix w_qpole;			/* weighted quadrupole moment		    */

void roughcenter(Body *btab)
{
    int i;
    real w_tot, w_b;
    vector tmpv;
    Body *b;

    w_tot = 0.0;
    CLRV(cm_pos);
    for (i = 0, b = btab; i < nbody; i++, b++) {
	w_b = (weight)(b, tsnap, i);
	if (w_b < 0.0)
	    warning("weight[%d] = %g < 0\n", i, w_b);
	if (w_b > 0.0) {
	    w_tot = w_tot + w_b;
	    MULVS(tmpv, Pos(b), w_b);
	    ADDV(cm_pos, cm_pos, tmpv);
	}
    }
    if (w_tot == 0.0)
	error("roughcenter: total weight is zero");
    DIVVS(cm_pos, cm_pos, w_tot);
    dprintf(1,"Roughcenter: %g %g %g\n",cm_pos[0], cm_pos[1], cm_pos[2]);
}

void findcenter(Body *btab)
{
    int i;
    Body *b;
    real w_b;
    vector tmpv;

    w_sum = 0.0;
    CLRV(w_pos);
    CLRV(w_vel);
    for (i = 0, b = btab; i < nbody; i++, b++) {
	w_b = (weight)(b, tsnap, i);
	if (w_b > 0.0 && (rcut <= 0.0 || distv(Pos(b), cm_pos) < rcut)) {
	    w_sum = w_sum + w_b;
	    MULVS(tmpv, Pos(b), w_b);
	    ADDV(w_pos, w_pos, tmpv);
	    MULVS(tmpv, Vel(b), w_b);
	    ADDV(w_vel, w_vel, tmpv);
	}
    }
    if (w_sum == 0.0)
	error("findcenter: weight within rcut is zero");
    DIVVS(w_pos, w_pos, w_sum);
    DIVVS(w_vel, w_vel, w_sum);
    dprintf(1,"Findcenter: %g %g %g\n",w_pos[0], w_pos[1], w_pos[2]);
}

void findmoment(Body *btab)
{
    int i;
    Body *b;
    real w_b;
    vector tmpv, pos_b;
    matrix tmpm;

    CLRM(w_qpole);
    for (i = 0, b = btab; i < nbody; i++, b++) {
	w_b = (weight)(b, tsnap, i);
	if (w_b > 0.0 && (rcut <= 0.0 || distv(Pos(b), cm_pos) < rcut)) {
	    SUBV(pos_b, Pos(b), w_pos);
	    MULVS(tmpv, pos_b, w_b);
	    OUTVP(tmpm, tmpv, pos_b);
	    ADDM(w_qpole, w_qpole, tmpm);
	}
    }
    DIVMS(w_qpole, w_qpole, w_sum);
    dprintf(1,"Findmoment:\n");
}

void findmoment_vel(Body *btab)
{
    int i;
    Body *b;
    real w_b;
    vector tmpv, pos_b;
    matrix tmpm;

    CLRM(w_qpole);
    for (i = 0, b = btab; i < nbody; i++, b++) {
	w_b = (weight)(b, tsnap, i);
	if (w_b > 0.0 && (rcut <= 0.0 || distv(Pos(b), cm_pos) < rcut)) {
	    SUBV(pos_b, Pos(b), w_pos);
	    MULVS(tmpv, pos_b, w_b);
	    OUTVP(tmpm, tmpv, pos_b);
	    ADDM(w_qpole, w_qpole, tmpm);
	}
    }
    DIVMS(w_qpole, w_qpole, w_sum);
    dprintf(1,"Findmoment_vel:\n");
}

vector oldframe[3] = {
    { 1.0, 0.0, 0.0, },
    { 0.0, 1.0, 0.0, },
    { 0.0, 0.0, 1.0, },
};

void snaptransform(stream out)
{
    vector frame[3], pos_b, vel_b, acc_b;
    Body *b;
    int i;

    eigenframe(frame, w_qpole);
    if (dotvp(oldframe[0], frame[0]) < 0.0)
	MULVS(frame[0], frame[0], -1.0);
    if (dotvp(oldframe[2], frame[2]) < 0.0)
	MULVS(frame[2], frame[2], -1.0);
    CROSSVP(frame[1], frame[2], frame[0]);
    printvec("e_x:", frame[0], out);
    printvec("e_y:", frame[1], out);
    printvec("e_z:", frame[2], out);
    for (b = btab; b < btab+nbody; b++) {
	SUBV(Pos(b), Pos(b), w_pos);
	SUBV(Vel(b), Vel(b), w_vel);
	for (i = 0; i < NDIM; i++) {
	    pos_b[i] = dotvp(Pos(b), frame[i]);
	    vel_b[i] = dotvp(Vel(b), frame[i]);
	    acc_b[i] = dotvp(Acc(b), frame[i]);
	}
	SETV(Pos(b), pos_b);
	SETV(Vel(b), vel_b);
	SETV(Acc(b), acc_b);
    }
    for (i = 0; i < NDIM; i++)
	SETV(oldframe[i], frame[i]);
}

#include "nrutil.h"

void eigenframe(vector frame[], matrix mat)
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
    for (i = 1; i <= 3; i++)
	for (j = 1; j <= 3; j++)
	    frame[i-1][j-1] = v[j][i];
}

void printvec(string name, vector vec, stream out)
{
    vector rtp;	/* radius - theta - phi */
    xyz2rtp(vec,rtp);
    fprintf(out,"%12s  %10.5f  %10.5f  %10.5f  %10.5f   %5.1f %6.1f\n",
	    name, rtp[0], vec[0], vec[1], vec[2],
	    rtp[1]*180.0/PI, rtp[2]*180.0/PI);
}

void xyz2rtp(vector xyz, vector rtp)
{
    real z = xyz[2];
    real w = sqrt(sqr(xyz[0])+sqr(xyz[1]));
    rtp[1] = atan(w/z);                 /* theta: in range 0 .. PI */
    if (z<0) rtp[1] += PI;
    rtp[2] = atan2(xyz[1], xyz[0]);     /* phi: in range  -PI .. PI */
    rtp[0] = sqrt(z*z+w*w);
    
}

