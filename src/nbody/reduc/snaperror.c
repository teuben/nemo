/*
 *  SNAPERROR: potential/force error calculations
 *
 *      24-mar-2005    Created              Peter Teuben
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n		   Input file (snapshot)",
    "model=0\n             Dehnen's MOD number",
    "gamma=1\n             Parameter of model (if any)",
    "VERSION=0.1\n	   24-mar-2005 PJT",
    NULL,
};

string usage="potential and force error calculations where analytical form is known";

string cvsid="$Id$";

void exact(real *pos, int model, real gamma, real *pot, real *acc);

void nemo_main(void)
{
    stream instr;
    real   tsnap;
    real gamma;
    bool Qpot;
    string headline=NULL;
    Body *btab = NULL, *bp;
    int i, j, nbody, bits, model;

    instr = stropen(getparam("in"), "r");	/* open input file */
    model = getiparam("model");
    gamma = getrparam("gamma");

    get_history(instr);                 /* read history */
    if (!get_tag_ok(instr, SnapShotTag))
      return;
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if ( (bits & PhaseSpaceBit) == 0)
      return;
    dprintf(1,"Processing snapshot @ 0x%x time=%g\n",btab,tsnap);
    Qpot = bits & PotentialBit;

    for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
      exact(Pos(bp),model,gamma,&Phi(bp),Acc(bp));
    }

    strclose(instr);
}


void exact(real *pos, int model, real gamma, real *pot, real *acc)
{
  real r, tm, P, fr, rh;

  r  = sqrt(sqr(pos[0]) + sqr(pos[1]) + sqr(pos[2]));
  switch(model) {
  case 0:  /* homogeneous */
    P  = 0.5*(3*r*r-1);
    fr =-1.;
    rh = 0.75/PI;
    break;
  case 1:  /* Plummer */
    tm = r*r;
    P  =-1./sqrt(tm+1);
    fr = P*P*P;
    rh = 0.75/PI*powi(-P,5);
    break;
#if 0
  case 2:  /* gamma-model */
    tm = pow(drand48(),ig3);
    r  = tm / (1.-tm);
    P  = (GAM==2.)?  log(tm) : -ig2*(1-pow(tm,g2));
    fr =-pow(tm,g1) * (1-tm) * (1-tm) / r;
    rh = 0.25*g3/PI * pow(1.-tm,4)/pow(tm,GAM);
    break;
  case 3:  /* Kuzmin disk */
    tm = drand48();
    r  = sqrt((2-tm)*tm)/(1.-tm);
    P  =-1./sqrt(r*r+1);
    fr = P*P*P;
    rh =-0.5/PI*fr;
    break;
#endif
  case 4: default: /* homogeneous disk */
    tm = r*r;
    P  = 2*(tm/3.-2);
    fr =-4*r/3.;
    rh = 1.0/PI;
    break;
  }
#if 0
  cth      = (MOD>=3) ? 0. : 2*drand48()-1;
  R        = r*sqrt(1.-cth*cth);
  phi      = TPi*drand48();
  F[i]     = 1;
  if (MOD >= 0) {
    X[0][i]  = R * sin(phi);
    X[1][i]  = R * cos(phi);
    X[2][i]  = r * cth;
  }
  M[i]     = mi;
  for(j=0; j<3; j++) AT[j][i] = fr*X[j][i]; /* "true" force   */
  RH[i] = rh;                               /* "true" density */
#endif
}
