/* 
 * MKTT77 - each ring in its own snapshot, use snapmerge to merge them
 *
 *   19-nov-2002	0.1 Created - Q&D
 *    4-dec-2002        0.2 Added mass=, eps= (mode= not implemetned)
 *    5-dec-2002        0.3 add central=, grow=
 *   19-feb-2006        0.4 dynamic memory using mdarray
 *    3-jun-2010        0.5 add ring number as Key
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <mdarray.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "out=???\n		output file name",
    "nbody=100\n	number of particles per (first) ring",
    "radius=1:6:1\n	radii of rings",
    "mass=1.0\n         central mass",
    "eps=0.0\n          softening length for central particle",
    "central=f\n        add the central mass also as 1st point?",
    "grow=t\n           grow number of points per ring to keep a constant density",
    "headline=\n	verbiage for output",
    "VERSION=0.5a\n	9-apr-2010 PJT",
    NULL,
};

string usage = "Create a Toomre & Toomre 1972 test disk centered around a point mass";

string cvsid="$Id$";


#define MAXRAD 1024


local int nobj, nobj_max, ntot = 0;
local real *mass;
local mdarray3 phase;
local int *key;
local double radius[MAXRAD];

local real eps2;
local real sqrtm;

local stream outstr;
local string headline;

local bool Qgrow;

void makering(int n, int k, real radius);
void makecenter(real m);
void writesnap(int n);


void nemo_main(void)
{
    int i, nrad, n;

    nrad = nemoinpd(getparam("radius"),radius,MAXRAD);
    nobj = getiparam("nbody");

    Qgrow = getbparam("grow");
    if (Qgrow) {
      nobj_max = (int) (nobj * (radius[nrad-1]/radius[0]));
      dprintf(0,"Grow nbody= %d -> %d\n",nobj,nobj_max);
    } else
      nobj_max = nobj;

    mass = (real *) allocate(sizeof(real)*nobj_max);
    phase = allocate_mdarray3(nobj_max,2,NDIM);
    key = (int *) allocate(sizeof(int)*nobj_max);
    headline = getparam("headline");
    sqrtm = getdparam("mass");
    if (getbparam("central")) {
      makecenter(sqrtm);
      writesnap(1);
    } else if (nobj==0)
      error("Cannot produce models with no points");
    sqrtm = sqrt(sqrtm);
    eps2 = getdparam("eps");
    eps2 = eps2*eps2;

    for (i=0; i<nrad; i++) {
      if (Qgrow)
	n = (int)(nobj*(radius[i]/radius[0]));
      else
	n = nobj;
      makering(n,i,radius[i]);
      writesnap(n);
      nemo_dprintf(1,"Ring %d %d:%d\n",i,ntot-n,ntot-1);
    }
    strclose(outstr);
    nemo_dprintf(1,"Total number of particles written: %d\n",ntot);
}

void makering(int n, int k, real radius)
{
  int i;
  real theta, velo = sqrtm*radius*pow(radius*radius+eps2,-0.75);

  for (i = 0; i < n; i++) {
    mass[i] = 0.0;
    theta = TWO_PI * ((real) i) / n;
    CLRV(phase[i][0]);
    phase[i][0][0] = radius * cos(theta);
    phase[i][0][1] = radius * sin(theta);
    CLRV(phase[i][1]);
    phase[i][1][0] = -velo * sin(theta);
    phase[i][1][1] =  velo * cos(theta);
    key[i] = k;
  }
}

void makecenter(real m)
{
  mass[0] = m;
  CLRV(phase[0][0]);
  CLRV(phase[0][1]);
}

void writesnap(int n)
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
      put_data(outstr, MassTag, RealType, mass, n, 0);
      put_data(outstr, PhaseSpaceTag, RealType, phase[0][0], n, 2, NDIM, 0);
      put_data(outstr, KeyTag, IntType, key, n, 0);
     put_tes(outstr, ParticlesTag);
    put_tes(outstr, SnapShotTag);
    ntot += n;
}
