/*
 * SNAPBINARY.C: analyse a pair of body's as a binary
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

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n			input file name",
    "i1=0\n                     First star (or list of stars)",
    "i2=1\n                     Second star (or list of stars > i2)",
    "bound=t\n                  Only show bound stars?",
    "VERSION=0.3\n	        7-mar-2019 PJT",
    NULL,
};

string usage="analyse binaries in a snapshot";


#define MAXNBODY 1000



void nemo_main(void)
{
  stream instr;
  Body *btab = NULL, *bp1, *bp2;
  int n1,i1[MAXNBODY],j1, n2,i2[MAXNBODY],j2;
  int nbody, bits;
  real tsnap, m1, m2, mu, T,W,E,a,b,L,period, w_sum, p, e;
  vector x, v, w_pos, w_vel, H;
  bool Qbound = getbparam("bound");

  instr = stropen(getparam("in"), "r");                  // get snapshot
  get_history(instr);
  get_snap(instr, &btab, &nbody, &tsnap, &bits);

  n1 = nemoinpi(getparam("i1"),i1,MAXNBODY);             // get list of stars-1
  n2 = nemoinpi(getparam("i2"),i2,MAXNBODY);             // get list of stars-2

  for (j2=0; j2<n2; j2++)                                // loop over the list
    for (j1=0; j1<n1; j1++) {
      bp1 = &btab[i1[j1]];                               // pointer to body
      bp2 = &btab[i2[j2]];
      dprintf(1,"i1,i2:    %d %d\n",i1[j1],i2[j2]);      
      if (i1[j1] > i2[j2]) continue;                     // skip
      if (i2[j2] >= nbody) continue;      
      if (i1[j1] >= nbody) {                             // or break when illegal
	error("nbody=%d\n",nbody);
	break;
      }
      
  
      m1 = Mass(bp1);
      m2 = Mass(bp2);
      mu = m1+m2;                                        // not reduced mass here!

      SUBV(x,Pos(bp1),Pos(bp2));
      SUBV(v,Vel(bp1),Vel(bp2));

      T = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);     // kinetic
      W = -mu/sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);   // potential
      E = T + W;
      a = -mu/2/E;                                       // semi major axis
      
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
      } else if (!Qbound) !
	printf("period:   Inf (Not a binary)!\n");
    }
}
