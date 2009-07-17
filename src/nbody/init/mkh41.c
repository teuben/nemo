/* 
 * MKH41 - each ring in its own snapshot, use snapmerge to merge them
 *
 * Some units from Holmberg's paper:
 *   mass     = 10^11 M_solar
 *   diameter = 2500 pc
 *
 *   15-jul-2009     1.0   Created at PiTP - Alar Toomre & Peter Teuben
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <mdarray.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "out=???\n		          output file name",
    "nbody=1,6,8,10,12\n          number of particles per ring",
    "radius=0,1,2,3,4\n	          radii of rings",
    "mass=1.0,1.0,1.0,0.7,0.3\n   masses of each particle",
    "phi=0,0,0,0,0\n              angles of first particle (deg)",
    "candlepower=1\n              Alar's candlepower scaling factor",
    "headline=\n	          verbiage for output",
    "VERSION=1.0\n                15-jul-09 PJT+AT",
    NULL,
};

string usage = "Create a Holmberg 1941 disk";

string cvsid="$Id$";


#define MAXRAD 1024


local int nobj, nobj_max, ntot = 0;
local mdarray3 pphase;
local real     *pmass;
local int    nbody[MAXRAD];
local real   radius[MAXRAD];
local real   mass[MAXRAD];
local real   phi[MAXRAD];
local real   cp;

local stream outstr;
local string headline;

local bool Qgrow;



void nemo_main()
{
    int i, nrad, n;

    cp = getdparam("candlepower");

    nrad = nemoinpi(getparam("nbody"),nbody,MAXRAD);
    n = nemoinpd(getparam("radius"),radius,MAXRAD);
    if (n!=nrad) error("radius=");
    n = nemoinpd(getparam("mass"),mass,MAXRAD);
    if (n!=nrad) error("mass=");
    n = nemoinpd(getparam("phi"),phi,MAXRAD);
    if (n!=nrad) error("phi=");

    nobj_max = nbody[0];
    for (i=1; i<nrad; i++)
      if (nbody[i] > nobj_max) nobj_max = nbody[i];

    pmass = (real *) allocate(sizeof(real)*nobj_max);
    pphase = allocate_mdarray3(nobj_max,2,NDIM);
    headline = getparam("headline");

    for (i=0; i<nrad; i++) {
      makering(nbody[i],mass[i],radius[i],phi[i]);
      writesnap(nbody[i]);
    }
    strclose(outstr);
    nemo_dprintf(1,"Total number of particles written: %d\n",ntot);
}

makering(int n, real m, real r, real p)
{
  int i;
  real theta, v = 0;

  if (r>0) v = r*pow(r*r,-0.75);         /* some fake keplerian vel */

  for (i = 0; i < n; i++) {
    pmass[i] = m;
    theta = TWO_PI * ( ((real) i) / n   + p/360.0 + 0.25);
    CLRV(pphase[i][0]);
    pphase[i][0][0] = r * cos(theta);
    pphase[i][0][1] = r * sin(theta);
    CLRV(pphase[i][1]);
    pphase[i][1][0] = -v * sin(theta);
    pphase[i][1][1] =  v * cos(theta);
  }
}


writesnap(int n)
{
    int cs = CSCode(Cartesian, NDIM, 2);
    static bool first = TRUE;

    if (n==0) return;

    if (first) {
        if (! streq(headline, ""))
            set_headline(headline);
        outstr = stropen(getparam("out"), "w");
        put_history(outstr);
	first = FALSE;
    }

    put_set(outstr, SnapShotTag);
     put_set(outstr, ParametersTag);
      put_data(outstr, NobjTag, IntType, &n, 0);
     put_tes(outstr, ParametersTag);
     put_set(outstr, ParticlesTag);
      put_data(outstr, CoordSystemTag, IntType, &cs, 0);
      put_data(outstr, MassTag, RealType, pmass, n, 0);
      put_data(outstr, PhaseSpaceTag, RealType, pphase[0][0], n, 2, NDIM, 0);
     put_tes(outstr, ParticlesTag);
    put_tes(outstr, SnapShotTag);
    ntot += n;
}
