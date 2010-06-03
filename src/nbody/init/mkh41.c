/* 
 * MKH41 - each ring in its own snapshot, use snapmerge to merge them
 *
 * Some units from Holmberg's paper:
 *   mass     = 10^11 M_solar
 *   diameter = 2500 pc
 *
 *   15-jul-2009     1.0   Created at PiTP - Alar Toomre & Peter Teuben
 *   22-jul-2009     1.1   added sign=
 *   13-apr-2010     1.2   added Dubinski's version as an option
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
    "sign=1\n                     Sign of angular momentum",
    "mode=0\n                     0=nemo 1=dubinski",
    "headline=\n	          verbiage for output",
    "VERSION=1.2\n                13-apr-10 PJT",
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
local real   sign_L;

local stream outstr;
local string headline;

local bool Qgrow;


typedef struct {
	float m, x, y, z, vx, vy, vz, ax, ay, az;
} phase;

local phase r[37], r0[38], r1[37];




void nemo_main()
{
    int i, nrad, n;
    int mode = getiparam("mode");

    if (mode==1) {
      dubinski();
      return;
    }

    sign_L = getdparam("sign");

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
    pphase[i][1][0] = -v * sin(theta) * sign_L;
    pphase[i][1][1] =  v * cos(theta) * sign_L;
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

dubinski(void)
{
  int i, j, nobj;
  double a, da, rad;
  double fx, fy, fz, fr;
  double dx, dy, dz, dr2, dr, vc;
  double mtot, epssq=0.25;
  
  nobj = 37;

  r[0].x = 0; r[0].y = 0; r[0].z = 0;
  r[0].vx = 0; r[0].vy = 0; r[0].vz = 0;

  a = M_PI/2.0;
  da = 2.0*M_PI/6.;
  rad = 1.0;
  for(i=1; i<7; i++) {
    r[i].x = rad*cos(a);
    r[i].y = rad*sin(a);
    a += da;
  }
  a = M_PI/2.0;
  da = 2.0*M_PI/8.;
  rad = 2.0;
  for(i=7; i<15; i++) {
    r[i].x = rad*cos(a);
    r[i].y = rad*sin(a);
    a += da;
  }
  a = M_PI/2.0;
  da = 2.0*M_PI/10.;
  rad = 3.0;
  for(i=15; i<25; i++) {
    r[i].x = rad*cos(a);
    r[i].y = rad*sin(a);
    a += da;
  }
  a = M_PI/2.0;
  da = 2.0*M_PI/12.;
  rad = 4.0;
  for(i=25; i<37; i++) {
    r[i].x = rad*cos(a);
    r[i].y = rad*sin(a);
    a += da;
  }
  for(i=0; i<37; i++) {
    r[i].m = 1.0/nobj;
    r[i].z = 0;
    r[i].vx = 0;
    r[i].vy = 0;
    r[i].vz = 0;
  }
  mtot = 0;
  for(i=0; i<15; i++) {
    r[i].m = 1.0;
    mtot += r[i].m;
  }
  for(i=15; i<25; i++) {
    r[i].m = 0.7;
    mtot += r[i].m;
  }
  for(i=25; i<37; i++) {
    r[i].m = 0.3;
    mtot += r[i].m;
  }
  for(i=0; i<37; i++) {
    r[i].m /= mtot;
  }
  
  for(i=0; i<nobj; i++) {
    fx = fy = fz = 0;
    for(j=0; j<nobj; j++) {
      if( j == i ) continue;
      dx = r[j].x - r[i].x;
      dy = r[j].y - r[i].y;
      dz = r[j].z - r[i].z;
      dr2 = dx*dx + dy*dy + dz*dz + epssq;
      dr = sqrt(dr2);
      fx += r[j].m*dx/(dr*dr2);
      fy += r[j].m*dy/(dr*dr2);
      fz += r[j].m*dz/(dr*dr2);
    }
    r[i].ax = fx; r[i].ay = fy; r[i].az = fz;
  }
  r[0].vx = r[0].vy = r[0].vz = 0;
  
  fr = 0;
  for(i=1; i<7; i++) {
    fx = r[i].ax; fy = r[i].ay; fz = r[i].az;
    dr = sqrt(r[i].x*r[i].x + r[i].y*r[i].y + r[i].z*r[i].z);
    fr += (fx*r[i].x + fy*r[i].y + fz*r[i].z)/dr;
  }
  fr /= 6;
  vc = sqrt(-fr);
  for(i=1; i<7; i++) {
    r[i].vx = -vc*r[i].y/dr;
    r[i].vy = vc*r[i].x/dr;
    r[i].vz = 0.0;
  }
  
  fr = 0;
  for(i=7; i<15; i++) {
    fx = r[i].ax; fy = r[i].ay; fz = r[i].az;
    dr = sqrt(r[i].x*r[i].x + r[i].y*r[i].y + r[i].z*r[i].z);
    fr += (fx*r[i].x + fy*r[i].y + fz*r[i].z)/dr;
  }
  fr /= 8;
  vc = sqrt(-2*fr);
  for(i=7; i<15; i++) {
    r[i].vx = -vc*r[i].y/dr;
    r[i].vy = vc*r[i].x/dr;
    r[i].vz = 0.0;
  }
  
  fr = 0;
  for(i=15; i<25; i++) {
    fx = r[i].ax; fy = r[i].ay; fz = r[i].az;
    dr = sqrt(r[i].x*r[i].x + r[i].y*r[i].y + r[i].z*r[i].z);
    fr += (fx*r[i].x + fy*r[i].y + fz*r[i].z)/dr;
  }
  fr /= 10;
  vc = sqrt(-3*fr);
  for(i=15; i<25; i++) {
    r[i].vx = -vc*r[i].y/dr;
    r[i].vy = vc*r[i].x/dr;
    r[i].vz = 0.0;
  }
  
  fr = 0;
  for(i=25; i<37; i++) {
    fx = r[i].ax; fy = r[i].ay; fz = r[i].az;
    dr = sqrt(r[i].x*r[i].x + r[i].y*r[i].y + r[i].z*r[i].z);
    fr += (fx*r[i].x + fy*r[i].y + fz*r[i].z)/dr;
  }
  fr /= 12;
  vc = sqrt(-4*fr);
  for(i=25; i<37; i++) {
    r[i].vx = -vc*r[i].y/dr;
    r[i].vy = vc*r[i].x/dr;
    r[i].vz = 0.0;
  }
  
/*
	for(i=0; i<37; i++) {
		fprintf(stdout,"%.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
			r[i].m,r[i].x,r[i].y,r[i].z,r[i].vx,r[i].vy,r[i].vz);
	}
*/


  
  {
    float m1, m2, mu1, mu2, vp, b, x, y, vx, vy, rsep;
    
    m1 = m2 = 1;
    rsep = 12;
    b = 9;
    mu1 = m2/(m1+m2);
    mu2 = -m1/(m1+m2);
    vp = sqrt(2*(m1 + m2)/b);
    
    x = (2*b - rsep);  y = 2*sqrt(b*(rsep-b));
    vx = -sqrt(b*(rsep-b))*vp/rsep; vy = +b*vp/rsep;
    fprintf(stderr,"%g %g %g %g\n",x,y,vx,vy);
    
    for(i=0; i<nobj; i++) {
      r0[i] = r[i];
      r1[i] = r[i];
      r0[i].x = r[i].x - mu1*x;
      r0[i].y = r[i].y + mu1*y;
      r0[i].vx = r[i].vx + mu1*vx;
      r0[i].vy = r[i].vy - mu1*vy;
      r1[i].x = -r[i].x - mu2*x;
      r1[i].y = -r[i].y + mu2*y;
      r1[i].vx = -r[i].vx + mu2*vx;
      r1[i].vy = -r[i].vy - mu2*vy;
    }
  }
  
/*
	fprintf(stdout,"74\n");
*/
  for(i=0; i<37; i++)
/*
		fwrite(r0+i,28,1,stdout);
*/
    fprintf(stdout,"%.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
	    r0[i].m,r0[i].x,r0[i].y,r0[i].z,r0[i].vx,r0[i].vy,r0[i].vz);
  for(i=0; i<37; i++) 
/*
		fwrite(r1+i,28,1,stdout);
*/
    fprintf(stdout,"%.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
	    r1[i].m,r1[i].x,r1[i].y,r1[i].z,r1[i].vx,r1[i].vy,r1[i].vz);
/*
	for(i=0; i<37; i++)
		fwrite(r0+i,28,1,stdout);
	for(i=0; i<37; i++)
		fwrite(r1+i,28,1,stdout);
*/
  
}
