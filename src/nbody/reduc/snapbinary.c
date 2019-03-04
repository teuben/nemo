/*
 * SNAPBINARY.C: analyse a pair of body's as a binary
 *
 * V0.1  2-mar-2019   PJT
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
    "i1=0\n                     First star",
    "i2=1\n                     Second star",
    "zerocm=f\n                 Center of mass correction?",
    "VERSION=0.2\n	        4-mar-2019 PJT",
    NULL,
};

string usage="analyse binaries in a snapshot";


void nemo_main(void)
{
  stream instr;
  Body *btab = NULL, *bp1, *bp2;
  int i1 = getiparam("i1");
  int i2 = getiparam("i2");
  bool Qcm = getbparam("zerocm");
  int nbody, bits;
  real tsnap, m1, m2, mu, T,W,E,a,b,L,period, w_sum, p, e;
  vector tmpx, tmpv, w_pos, w_vel, H;

  instr = stropen(getparam("in"), "r");
  get_history(instr);
  get_snap(instr, &btab, &nbody, &tsnap, &bits);
  bp1 = &btab[i1];
  bp2 = &btab[i2];
  if (Qcm) {
    w_sum = Mass(bp1) + Mass(bp2);
    CLRV(w_pos);
    CLRV(w_vel);
    MULVS(tmpx,Pos(bp1),Mass(bp1));
    ADDV(w_pos,w_pos,tmpx);
    MULVS(tmpv,Vel(bp1),Mass(bp1));
    ADDV(w_vel,w_vel,tmpv);
    MULVS(tmpx,Pos(bp2),Mass(bp2));
    ADDV(w_pos,w_pos,tmpx);
    MULVS(tmpv,Vel(bp2),Mass(bp2));
    ADDV(w_vel,w_vel,tmpv);
    SDIVVS(w_pos,w_sum);
    SDIVVS(w_vel,w_sum);
    printf("COM: %g %g %g %g %g %g\n" , w_pos[0],w_pos[1],w_pos[2],w_vel[0],w_vel[1],w_vel[2]);
    SUBV(Pos(bp1),Pos(bp1),w_pos);
    SUBV(Vel(bp1),Vel(bp1),w_vel);    
    SUBV(Pos(bp2),Pos(bp2),w_pos);
    SUBV(Vel(bp2),Vel(bp2),w_vel);    
  }

  m1 = Mass(bp1);
  m2 = Mass(bp2);
  mu = m1+m2;
  printf("m1,m2,mu: %g %g %g\n",m1,m2,mu);
  SUBV(tmpx,Pos(bp1),Pos(bp2));
  SUBV(tmpv,Vel(bp1),Vel(bp2));
  printf("pos,vel:  %g %g %g %g %g %g\n", tmpx[0],tmpx[1],tmpx[2],tmpv[0],tmpv[1],tmpv[2]);
  T = 0.5 * (tmpv[0]*tmpv[0] + tmpv[1]*tmpv[1] + tmpv[2]*tmpv[2]);
  W = -mu/sqrt(tmpx[0]*tmpx[0] + tmpx[1]*tmpx[1] + tmpx[2]*tmpx[2]);
  E = T + W;
  a = -mu/2/E;
  printf("T,W,E:    %g %g %g \n",T,W,E);

  CROSSVP(H,tmpx,tmpv);
  ABSV(L,H);
  printf("H,|H|:    %g %g %g %g\n",H[0],H[1],H[2],L);
  p = L*L/mu;
  e = sqrt(1 - p/a);
  if (e<1)
    b = a * sqrt(1-e*e);
  else
    b = a * sqrt(e*e-1);
  printf("a,b,p,e:  %g %g %g %g\n",a,b,p,e);
  
  if (a>0) {
    period = TWO_PI * a * sqrt(a/mu);
    printf("period:   %g\n",period);
  } else
    printf("period:   Inf (Not a binary)!\n");
    

}
