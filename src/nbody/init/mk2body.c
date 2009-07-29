/* 
 * MK2BODY: make a 2 body system
 *
 *	29-jul-09	PJT     why did this take so long....
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <mdarray.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "out=???\n		output file name",
    "m=0.1\n            mass of 2nd particle",
    "x=1\n              X position of 2nd particle",
    "y=0\n              Y position of 2nd particle",
    "z=0\n              Z position of 2nd particle",
    "vx=0\n             VX velocity of 2nd particle",
    "vy=1\n             VY velocity of 2nd particle",
    "vz=0\n             VY velocity of 2nd particle",
    "zerocm=true\n	if true, zero the center of mass",
    "headline=\n	verbiage for output",
    "VERSION=1.0\n	29-jul-09 PJT",
    NULL,
};

string usage = "make a 2 body system";

string cvsid="$Id$";


local int nobj;
local real *m;
local mdarray3 phase;
local double radius;	                /* must be double !! */

extern double xrandom(double, double);



void nemo_main()
{
  int i,j;
  real rsq, pot, kin;

  nobj = 2;
      
  m = (real *) allocate(nobj * sizeof(real));
  phase = allocate_mdarray3(nobj,2,NDIM);

  m[0] = 1.0;             /* central particle, mass=1 at center */
  for (j=0; j<2; j++)
    for (i=0; i<NDIM; i++)
      phase[0][j][i] = 0.0;

  m[1] = getdparam("m");
  phase[1][0][0] = getdparam("x");
  phase[1][0][1] = getdparam("y");
  phase[1][0][2] = getdparam("z");
  phase[1][1][0] = getdparam("vx");
  phase[1][1][1] = getdparam("vy");
  phase[1][1][2] = getdparam("vz");

  if (getbparam("zerocm"))
    zerocms(phase[0][0], 2 * NDIM, m, nobj, nobj);

  rsq = sqrt( sqr(phase[0][0][0]-phase[1][0][0]) + 
	      sqr(phase[0][0][1]-phase[1][0][1]) + 
	      sqr(phase[0][0][2]-phase[1][0][2]));
  pot = -m[0]*m[1]/rsq;
  kin = 0.5*m[0]*(sqr(phase[0][1][0]) + sqr(phase[0][1][1]) + sqr(phase[0][1][2])) + 
        0.5*m[1]*(sqr(phase[1][1][0]) + sqr(phase[1][1][1]) + sqr(phase[1][1][2]));

  dprintf(0,"U: %g   T: %g   T+U: %g  T/U: %g\n", pot,kin,pot+kin,kin/pot);

  writesnap();
}


writesnap()
{
    stream outstr;
    string headline = getparam("headline");
    int cs = CSCode(Cartesian, NDIM, 2);

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(getparam("out"), "w");

    put_history(outstr);
    put_set(outstr, SnapShotTag);
     put_set(outstr, ParametersTag);
      put_data(outstr, NobjTag, IntType, &nobj, 0);
     put_tes(outstr, ParametersTag);
     put_set(outstr, ParticlesTag);
      put_data(outstr, CoordSystemTag, IntType, &cs, 0);
      put_data(outstr, MassTag, RealType, m, nobj, 0);
      put_data(outstr, PhaseSpaceTag, RealType, phase[0][0], nobj, 2, NDIM, 0);
     put_tes(outstr, ParticlesTag);
    put_tes(outstr, SnapShotTag);
    strclose(outstr);
}
