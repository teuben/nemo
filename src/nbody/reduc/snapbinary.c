/*
 * SNAPBINARY.C: analyse a pair of body's as a binary, purely based on pos and vel
 *
 * See also $NEMO/src/nbody/evolve/aarseth/nbody1/source/bodies.f
 *
 * V0.1  2-mar-2019   created - PJT
 * V0.3  7-mar-2019   added list option to i1,i2   - PJT
 */

#include <stdinc.h>
#include <vectmath.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>
#include <vectmath.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n			input file name",
    "i1=0\n                     First star (or list of stars)",
    "i2=1\n                     Second star, must be > i1 (or list of stars > i2)",
    "bound=t\n                  Only show bound stars?",
    "out=\n                     Optional out with strongest bound pair sink'd",
    "VERSION=0.5\n	        9-apr-2019 PJT",
    NULL,
};

string usage="analyse binaries in a snapshot, optionally merge a selected pair";


#define MAXNBODY 10000



void nemo_main(void)
{
  stream instr, outstr;
  Body *btab = NULL, *bp, *bp1, *bp2;
  int n1,i1[MAXNBODY],j1, n2,i2[MAXNBODY],j2;
  int nbody, bits, i1_min,i2_min, nch=0, nbn=0;
  real tsnap, m1, m2, mu, T,W,E,a,b,L,period, w_sum, p, e, w_min;
  vector x, v, w_pos, w_vel, H;
  bool Qbound = getbparam("bound");

  instr = stropen(getparam("in"), "r");                  // get snapshot
  get_history(instr);
  get_snap(instr, &btab, &nbody, &tsnap, &bits);

  if (hasvalue("out")) outstr = stropen(getparam("out"), "w");

  n1 = nemoinpi(getparam("i1"),i1,MAXNBODY);             // get list of stars-1
  n2 = nemoinpi(getparam("i2"),i2,MAXNBODY);             // get list of stars-2

  w_min  = 0.0;
  i1_min = -1;
  i2_min = -1;

  for (j2=0; j2<n2; j2++)                                // loop over the list
    for (j1=0; j1<n1; j1++) {
      bp1 = &btab[i1[j1]];                               // pointer to body
      bp2 = &btab[i2[j2]];
      dprintf(2,"i1,i2:    %d %d\n",i1[j1],i2[j2]);
      if (i1[j1] >= i2[j2])continue;                     // skip
      if (i2[j2] >= nbody) continue;      
      if (i1[j1] >= nbody) {                             // or break when illegal
	error("nbody=%d\n",nbody);
	break;
      }
      nch++;                                             // count #checks
      
      m1 = Mass(bp1);
      m2 = Mass(bp2);
      mu = m1+m2;                                        // not reduced mass here!

      SUBV(x,Pos(bp1),Pos(bp2));
      SUBV(v,Vel(bp1),Vel(bp2));

      T = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);     // kinetic
      W = -mu/sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);   // potential
      E = T + W;
      a = -mu/2/E;                                       // semi major axis
      if (a>0 && W < w_min) {
	w_min = W;                   // find max bound cluster
	i1_min = i1[j1];
	i2_min = i2[j2];
	dprintf(1,"W_min:    %g %d %d\n",w_min, i1[j1],i2[j2]);
      } else 
	dprintf(2,"Unbound a=%g %d %d\n",a,i1[j1],i2[j2]);

      CROSSVP(H,x,v);
      ABSV(L,H);

      p = L*L/mu;
      e = sqrt(1 - p/a);                                 // eccentricity
      if (e<1)
	b = a * sqrt(1-e*e);                             // ellipse
      else 
	b = a * sqrt(e*e-1);                             // hyperbola
      
      if (a>0 || !Qbound) {
	printf("i1,i2:    %d %d\n",i1[j1],i2[j2]);
	printf("m1,m2,mu: %g %g %g\n",m1,m2,mu);
	printf("pos,vel:  %g %g %g %g %g %g\n", x[0],x[1],x[2],v[0],v[1],v[2]);
	printf("T,W,E:    %g %g %g \n",T,W,E);
	printf("H,|H|:    %g %g %g %g\n",H[0],H[1],H[2],L);
	printf("a,b,p,e:  %g %g %g %g\n",a,b,p,e);
      }
      
      if (a>0) {
	period = TWO_PI * a * sqrt(a/mu);
	printf("period:   %g\n",period);
	nbn++;
      } else if (!Qbound) !
	printf("period:   Inf (Not a binary)!\n");
    }
  dprintf(1,"Checked %d pairs, %d were bound\n",nch,nbn);
  printf("W_min:    %g for %d %d (found %d bound pairs)\n",w_min,i1_min,i2_min,nbn);
  
  if (hasvalue("out")) {
    bp1 = &btab[i1_min];
    bp2 = &btab[i2_min];

    w_sum = Mass(bp1) + Mass(bp2);
    CLRV(w_pos);
    CLRV(w_vel);
    MULVS(x,Pos(bp1),Mass(bp1));
    ADDV(w_pos,w_pos,x);
    MULVS(v,Vel(bp1),Mass(bp1));
    ADDV(w_vel,w_vel,v);
    MULVS(x,Pos(bp2),Mass(bp2));
    ADDV(w_pos,w_pos,x);
    MULVS(v,Vel(bp2),Mass(bp2));
    ADDV(w_vel,w_vel,v);
    SDIVVS(w_pos,w_sum);
    SDIVVS(w_vel,w_sum);
    dprintf(1,"COM: %g %g %g %g %g %g\n" , w_pos[0],w_pos[1],w_pos[2],w_vel[0],w_vel[1],w_vel[2]);

    // shift over the masses so bp2 is removed
    Mass(bp1) = w_sum;
    SETV(Pos(bp1),w_pos);
    SETV(Vel(bp1),w_vel);
    for (bp=bp2; bp<btab+nbody-1; bp++) {
      bp1 = bp+1;
      dprintf(1,"Fixing tail end %g <- %g\n",Mass(bp),Mass(bp1));
      Mass(bp) = Mass(bp1);
      SETV(Pos(bp),Pos(bp1));
      SETV(Vel(bp),Vel(bp1));
    }
    nbody--;
    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    strclose(outstr);
  }
}
