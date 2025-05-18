/*
 * SNAPCENTERS: shrinking sphere algorithm
 *

 (from Power et. al 2003, S2.5)
 
The centre of each halo is determined using an iterative technique in
which the centre of mass of particles within a shrinking sphere is
computed recursively until a convergence criterion is met. At each
step of the iteration, the centre of the sphere is reset to the last
computed barycentre and the radius of the sphere is reduced by 2.5 per
cent. The iteration is stopped when a specified number of particles
(typically either 1000 particles or 1 per cent of the particles within
the high-resolution region, whichever is smaller) is reached within
the sphere.

Halo centres identified with this procedure are quite independent of
the parameters chosen to initiate the iteration, provided that the
initial sphere is large enough to encompass a large fraction of the
system. In a multi-component system, such as a dark halo with
substructure, this procedure isolates the densest region within the
largest subcomponent. In more regular systems, the centre so obtained
is in good agreement with centres obtained by weighing the centre of
mass by the local density or gravitational potential of each
particle. We have explicitly checked that none of the results
presented here are biased by our particular choice of centering
procedure.

NOTE: the authors make no arugments why this elaborate scheme is
favored to e.g. a weighted density/potential method - see our table in
hackdens.1

 * 
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
    "eta=0.001\n    Convergence stop criterion",
    "fn=0.025\n     Reduction fraction for sphere per iteration",
    "iter=20\n      Maximum number of iterations to use",
    "center=0,0,0\n Initial estimate for the center",
    "nmax=1000\n    Stop iteration for this #particles too",
    "rmax=10\n      Initial radius to shrink from",
    "one=f\n        Only output COM as a snapshot? [not implemented]",
    "VERSION=0.1\n  16-may-2025 PJT",
    NULL,
};

string usage="Center position of a snapshot based on shrinking sphere method (Power2003)";


int snapcenter(Body*, int, real, rproc_body, vector, vector, int, real, bool);

void nemo_main()
{
  stream instr, outstr;
  string times;
  rproc_body weight;
  Body *btab = NULL;
  int i, j, np=0, nbody=0, bits, iter;
  real tsnap, eta, dr;
  bool Qreport, Qone;
  vector n_pos, n_vel, o_pos, o_vel;
  real pos[NDIM];
  int nmax, nleft, nleft0, mode;
  real rmax, fn;

  warning("Under development");
    
  instr = stropen(getparam("in"), "r");
  outstr = stropen(getparam("out"), "w");
  weight = btrtrans(getparam("weight"));
  nmax = getiparam("nmax");
  rmax = getrparam("rmax");
  fn = getrparam("fn");
  eta = getrparam("eta");
  iter = getiparam("iter");
  times = getparam("times");
  Qone = getbparam("one");
  if (Qone) dprintf(1,"Output single particle COM snapshot"); 
  Qreport = getbparam("report");
  if (Qreport) dprintf(1,"pos vel of center(s) will be:\n");
  if (hasvalue("center")) {
    np = nemoinpr(getparam("center"),pos,NDIM);
    if (np != NDIM) error("center=%s needs %d values, got %d", getparam("center"), NDIM, np);
    dprintf(1,"Using center=%g,%g,%g\n", pos[0],pos[1],pos[2]);
  } else
    np = 0;

  
  get_history(instr);
  put_history(outstr);
  
  do {
    get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);
    if (bits & PhaseSpaceBit) {
      CLRV(n_pos);
      CLRV(n_vel);
      if (np > 0) SETV(n_pos, pos);
      nleft0 = nbody;
      mode=0;
      for (i=0;i<iter;i++) {
	nleft = snapcenter(btab, nbody, tsnap, weight, n_pos, n_vel, nmax, rmax, Qreport);
	if (i>0) {
	  dr = distv(o_pos,n_pos);
	  dprintf(1,"%d ",i);
	  for (j=0; j<NDIM; j++)  dprintf(1,"%f ",n_pos[j]);
	  //for (j=0; j<NDIM; j++)  dprintf(1,"%f ",n_vel[j]);
	  dprintf(1,"%f\n", dr);
	  mode=1;
	  if (dr < eta) break;
	  if (i == iter-1) 
	    warning("eta=%g too small?  dr=%g after iter=%d\n",eta,dr,iter);
	  mode=2;	  
	  if (nleft < nmax) break;
	  mode=3;	  
	  if (nleft == nleft0) break;
	  nleft0 = nleft;
	} 
	SETV(o_pos,n_pos);
	SETV(o_vel,n_vel);

	rmax *= (1-fn);
	
      }
      if (Qreport) {
	printf("%d  ",mode);
	for (j=0; j<NDIM; j++)  printf("%f ",n_pos[j]);
	//for (j=0; j<NDIM; j++)  printf("%f ",n_vel[j]);
	printf("\n");
      }
      // write output
      // warning("development version; no output data written yet.");
    }
  } while (bits != 0);
}

int snapcenter(
	       Body *btab,
	       int nbody,
	       real tsnap,
	       rproc_body weight,
	       vector o_pos, 
	       vector o_vel,
	       int nmax,
	       real rmax,
	       bool Qreport)
{
    int i;
    Body *b;
    real s, w_i, w_sum;
    vector tmpv, w_pos, w_vel;

    w_sum = 0.0;
    CLRV(w_pos);
    CLRV(w_vel);
    int cnt=0;
    for (i = 0, b = btab; i < nbody; i++, b++) {   
	SUBV(tmpv,o_pos,Pos(b));
	DOTVP(s,tmpv,tmpv);
	if (s > rmax*rmax) continue;

	cnt++;

	w_i = (weight)(b, tsnap, i);    // should be mass
	if (w_i < 0.0) warning("weight[%d] = %g < 0\n", i, w_i);
	w_i /= s;

	w_sum += w_i;
	MULVS(tmpv, Pos(b), w_i);
	SADDV(w_pos, tmpv);
	MULVS(tmpv, Vel(b), w_i);  // has no meaning
	SADDV(w_vel, tmpv);
    }
    if (w_sum == 0.0) error("total weight is zero");
    dprintf(0,"Rmax=%g n=%d\n",rmax,cnt);
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
    return cnt;
}
