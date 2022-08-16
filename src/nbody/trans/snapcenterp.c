/*
 * SNAPCENTERP: find center of a snapshot with the Cruz et al. method
 *
 *       1-apr-06   0.1 no joke, hotel rembrandt amsterdam,       pjt
 *      12-aug-22   0.2 cleanup, only report pos now (no vel)     pjt
 *      13-aug-22   0.3 more cleanup, still no proper convergence PJT
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
    "in=???\n       input file name",
    "out=???\n      output file name",
    "weight=m\n	    weight factor used finding center",
    "times=all\n    range of times to process",
    "report=f\n	    report the c.o.m shift",
    "eps=0.025\n    Gravitational softening length",
    "eta=0.001\n    Convergence stop criterion",
    "fn=0.5\n       Fraction of particles to consider (not used yet)",
    "iter=20\n      Maximum number of iterations to use",
    "VERSION=0.3\n  13-aug-2022 PJT",
    NULL,
};

string usage="Center position of a snapshot based on iterative Cruz2002 method";


void snapcenter(Body*, int, real, rproc_body, real, vector, vector, bool);

void nemo_main()
{
  stream instr, outstr;
  string times;
  rproc_body weight;
  Body *btab = NULL;
  int i, j, nbody=0, bits, iter;
  real tsnap, mass, eps, eta, dr;
  bool Qreport;
  vector n_pos, n_vel, o_pos, o_vel;
  
  instr = stropen(getparam("in"), "r");
  outstr = stropen(getparam("out"), "w");
  weight = btrtrans(getparam("weight"));
  eps = getrparam("eps");
  eta = getrparam("eta");
  iter = getiparam("iter");
  times = getparam("times");
  Qreport = getbparam("report");
  if (Qreport) dprintf(1,"pos vel of center(s) will be:\n");
  
  get_history(instr);
  put_history(outstr);
  
  do {
    get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);
    if (bits & PhaseSpaceBit) {
      CLRV(n_pos);
      CLRV(n_vel);
      for (i=0;i<iter;i++) {
	snapcenter(btab, nbody, tsnap, weight, eps, n_pos, n_vel, Qreport);
	if (i>0) {
	  dr = distv(o_pos,n_pos);
	  dprintf(1,"%d ",i);
	  for (j=0; j<NDIM; j++)  dprintf(1,"%f ",n_pos[j]);
	  //for (j=0; j<NDIM; j++)  dprintf(1,"%f ",n_vel[j]);
	  dprintf(1,"%f\n", dr);
	  if (dr < eta) break;
	  if (i == iter-1) 
	    warning("eta=%g too small?  dr=%g after iter=%d\n",eta,dr,iter);
	} 
	SETV(o_pos,n_pos);
	SETV(o_vel,n_vel);
	
      }
      if (Qreport) {
	for (j=0; j<NDIM; j++)  printf("%f ",n_pos[j]);
	//for (j=0; j<NDIM; j++)  printf("%f ",n_vel[j]);
	printf("\n");
      }
      // write output
      // warning("development version; no output data written yet.");
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
    for (i = 0, b = btab; i < nbody; i++, b++) {   // Cruz eq.(4)
	SUBV(tmpv,o_pos,Pos(b));
	DOTVP(s,tmpv,tmpv);
	s += eps2;
	s = s * sqrt(s);  // @todo   could use another power?

	w_i = (weight)(b, tsnap, i);
	if (w_i < 0.0) warning("weight[%d] = %g < 0\n", i, w_i);
	w_i /= s;

	w_sum += w_i;
	MULVS(tmpv, Pos(b), w_i);
	ADDV(w_pos, w_pos, tmpv);
	MULVS(tmpv, Vel(b), w_i);
	ADDV(w_vel, w_vel, tmpv);
    }
    if (w_sum == 0.0) error("total weight is zero");
    SDIVVS(w_pos, w_sum);
    SDIVVS(w_vel, w_sum);

    SETV(o_pos,w_pos);
    SETV(o_vel,w_vel);

#if 0
    for (i = 0, b = btab; i < nbody; i++, b++) {
	SSUBV(Pos(b), w_pos);
	SSUBV(Vel(b), w_vel);
    }
#endif
}
