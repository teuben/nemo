/*
 * SNAPCENTERP: find center of a snapshot with the Cruz et al. method
 *
 *       1-apr-06   0.1 no joke, hotel rembrandt amsterdam,      pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <history.h>
#include <filestruct.h>

#include <snapshot/body.h>
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>
#include <bodytransc.h>

string defv[] = {	
    "in=???\n       input file name ",
    "weight=m\n	    factor used finding center",
    "times=all\n    range of times to process",
    "report=t\n	    report the c.o.m shift",
    "one=f\n        Only output COM as a snapshot?",
    "eps=0.025\n    Softening",
    "eta=0.001\n    Convergence stop criterion",
    "fn=0.5\n       Fraction of particles to consider",
    "iter=10\n      Maximum number of iterations to use",
    "VERSION=0.1\n  1-apr-06 PJT",
    NULL,
};

string usage="Center a snapshot based on iterative Cruz method";


void snapcenter(Body*, int, real, rproc_body, real, vector, vector, bool);

void nemo_main()
{
  stream instr, outstr;
  string times;
  rproc_body weight;
  Body *btab = NULL, *b;
  int i, nbody, bits, iter;
  real tsnap, mass, eps, eta, dr;
  bool Qreport, Qone;
  vector n_pos, n_vel, o_pos, o_vel;
  char line[256];
  
  instr = stropen(getparam("in"), "r");
  weight = btrtrans(getparam("weight"));
  eps = getrparam("eps");
  eta = getrparam("eta");
  iter = getiparam("iter");
  times = getparam("times");
  Qreport = getbparam("report");
  Qone = getbparam("one");
  if (Qreport) dprintf(1,"pos vel of center(s) will be:\n");
  
  get_history(instr);
  
  do {
    get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);
    if (bits & PhaseSpaceBit) {
      CLRV(n_pos);
      CLRV(n_vel);
      for (i=0;i<iter;i++) {
	snapcenter(btab, nbody, tsnap, weight, eps, n_pos, n_vel, Qreport);
	if (i>0) {
	  dr = distv(o_pos,n_pos);
	  dprintf(1,"dr=%g\n",dr);
	  if (dr < eta) break;
	  if (i == iter-1) 
	    warning("eta=%g too small?  dr=%g after iter=%d\n",eta,dr,iter);
	} 
	SETV(o_pos,n_pos);
	SETV(o_vel,n_vel);
      } 
    }
  } while (bits != 0);
}

void snapcenter(
		Body *btab,
		int nbody,
		real tsnap,
		rproc_body weight,
		real eps,
		vector o_pos, 
		vector o_vel,
		bool Qreport)
{
    int i;
    Body *b;
    real s, w_i, w_sum,eps2 = eps*eps;
    vector tmpv, w_pos, w_vel;

    w_sum = 0.0;
    CLRV(w_pos);
    CLRV(w_vel);
    for (i = 0, b = btab; i < nbody; i++, b++) {
	SUBV(tmpv,o_pos,Pos(b));
	DOTVP(s,tmpv,tmpv);
	s += eps2;

	w_i = (weight)(b, tsnap, i);
	if (w_i < 0.0)
	    warning("weight[%d] = %g < 0\n", i, w_i);
	w_i /= s;     /* potentially we could try pow(k) here */

	w_sum += w_i;
	MULVS(tmpv, Pos(b), w_i);
	ADDV(w_pos, w_pos, tmpv);
	MULVS(tmpv, Vel(b), w_i);
	ADDV(w_vel, w_vel, tmpv);
    }
    if (w_sum == 0.0)
	error("total weight is zero");
    SDIVVS(w_pos, w_sum);
    SDIVVS(w_vel, w_sum);

    if (Qreport) {
      for (i=0; i<NDIM; i++)
        printf("%f ",w_pos[i]);
      for (i=0; i<NDIM; i++)
        printf("%f ",w_vel[i]);
      printf("\n");
    } else {
      for (i=0; i<NDIM; i++)
        dprintf(1,"%f ",w_pos[i]);
      for (i=0; i<NDIM; i++)
        dprintf(1,"%f ",w_vel[i]);
      dprintf(1,"\n");
    }
    SETV(o_pos,w_pos);
    SETV(o_vel,w_vel);

#if 0
    for (i = 0, b = btab; i < nbody; i++, b++) {
	SSUBV(Pos(b), w_pos);
	SSUBV(Vel(b), w_vel);
    }
#endif
}
