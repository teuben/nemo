/*
 *  Compare distance between points in two tables
 *
 *      13-feb-2013     0.1    Created, Q&D
 *
 *  See also the SIRTF/map2 project in 'c2d'
 *
 *  @todo   consider allowing scaling/offset between two datasets 
 *          Y = aX + b
 *          use k-d Match?   http://arxiv.org/abs/1304.0838
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "in1=???\n		First (ascii table) dataset",
    "in2=???\n		Second (ascii table) dataset",
    "col1=1,2\n         X (optionally Y and Z) Column from 1st dataset",
    "col2=1,2\n         X (optionally Y and Z) Column from 2nd dataset",
    "radec=false\n      Interpret X and Y as RA/DEC angles on sky",
    "VERSION=0.3\n	14-nov-2013 PJT",
    NULL,
};

string usage="Compare distances between points in two tables";
string cvsid="$Id$";

void compare_1d(int npt1, int npt2, real *x1, real *x2);
void compare_2d(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2);
void compare_2da(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2);
void compare_3d(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, real *z1, real *z2);

inline real dist1(real x1,real x2) {
  real d = x1-x2;
  if (d > 0.0) return d;
  return -d;
}

inline real dist2(real x1,real y1,real x2,real y2) {
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

inline real dist3(real x1,real y1,real z1,real x2,real y2,real z2) {
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}

inline real dist2a(real x1,real y1,real x2,real y2) {
  return sin(y1)*sin(y2)+cos(y1)*cos(y2)*cos(x1-x2);
}

#define NDIM   3

nemo_main()
{
  int i, npt, npt1, npt2, ncol1, ncol2;
  real *coldat1[NDIM], *coldat2[NDIM];
  int colnr1[NDIM], colnr2[NDIM];
  real *x1 = NULL, *x2 = NULL;
  real *y1 = NULL, *y2 = NULL;
  real *z1 = NULL, *z2 = NULL;
  string input1 = getparam("in1");
  string input2 = getparam("in2");
  stream instr1, instr2;
  bool Qradec = getbparam("radec");

  if (Qradec) warning("new radec mode ");

  ncol1 = nemoinpi(getparam("col1"),colnr1,NDIM);
  ncol2 = nemoinpi(getparam("col2"),colnr2,NDIM);

  if (ncol1 != ncol2) error("column count error: ncol1=%d and ncol2=%d not same",ncol1,ncol2);
  if (ncol1 < 1 || ncol1 > NDIM) error("ncol1=%d not supported",ncol1);
  
  instr1 = stropen(input1,"r");
  instr2 = stropen(input2,"r");

  npt1 = nemo_file_lines(input1,0);
  npt2 = nemo_file_lines(input2,0);

  if (ncol1 > 0) {
    x1 = (real *) allocate(npt1 * sizeof(real));
    x2 = (real *) allocate(npt2 * sizeof(real));
  }
  if (ncol1 > 1) {
    y1 = (real *) allocate(npt1 * sizeof(real));
    y2 = (real *) allocate(npt2 * sizeof(real));
  }
  if (ncol1 > 2) {
    z1 = (real *) allocate(npt1 * sizeof(real));
    z2 = (real *) allocate(npt2 * sizeof(real));
  }
  
  coldat1[0] = x1;     coldat2[0] = x2;   
  coldat1[1] = y1;     coldat2[1] = y2;   
  coldat1[2] = z1;     coldat2[2] = z2;   

  npt1 = get_atable(instr1,ncol1,colnr1,coldat1,npt1);
  npt2 = get_atable(instr2,ncol2,colnr2,coldat2,npt2);

  dprintf(0,"Npt1=%d Npt2=%d\n",npt1,npt2);
  
  strclose(instr1);
  strclose(instr2);

  if (Qradec && ncol1==2)
    compare_2da(npt1, npt2, x1, x2, y1, y2);
  else if (ncol1 == 1) 
    compare_1d(npt1, npt2, x1, x2);
  else if (ncol1 == 2) 
    compare_2d(npt1, npt2, x1, x2, y1, y2);
  else if (ncol1 == 3)
    compare_3d(npt1, npt2, x1, x2, y1, y2, z1, z2);
  else
    warning("%d columns not supported",ncol1);
}




void compare_1d(int npt1, int npt2, real *x1, real *x2)
{
  int i, j, jmin;
  real d, dmin;
  
  for (i=0; i<npt1; i++) {
    dmin = dist1(x1[i],x2[0]);
    jmin = 0;
    for (j=1; j<npt2; j++) {
      d = dist1(x1[i],x2[j]);
      if (d<dmin) {
	dmin = d;
	jmin = j;
      }
    }
    printf("%d %g   %d %g   %g\n",i+1,x1[i],jmin+1,x2[jmin], dmin);
  }
}

void compare_2d(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2)
{  
  int i, j, jmin;
  real d, dmin;
  
  for (i=0; i<npt1; i++) {
    dmin = dist2(x1[i],y1[i],x2[0],y2[0]);
    jmin = 0;
    for (j=1; j<npt2; j++) {
      d = dist2(x1[i],y1[i],x2[j],y2[j]);
      if (d<dmin) {
	dmin = d;
	jmin = j;
      }
    }
    dmin = sqrt(dmin);
    printf("%d %g %g   %d %g %g   %g\n",i+1,x1[i],y1[i],jmin+1,x2[jmin], y2[jmin],dmin);
  }
}

void compare_2da(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2)
{  
  int i, j, jmin;
  real d, dmax, dmin;
  
  for (i=0; i<npt1; i++) {
    dmax = dist2a(x1[i],y1[i],x2[0],y2[0]);
    jmin = 0;
    for (j=1; j<npt2; j++) {
      d = dist2a(x1[i],y1[i],x2[j],y2[j]);
      if (d>dmax) {
	dmax = d;
	jmin = j;
      }
    }
    dmin = acos(dmax)*206265;  /* in arcsec */
    printf("%d %g %g   %d %g %g   %g\n",i+1,x1[i],y1[i],jmin+1,x2[jmin], y2[jmin],dmin);
  }
}

void compare_3d(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, real *z1, real *z2)
{  
  int i, j, jmin;
  real d, dmin;
  
  for (i=0; i<npt1; i++) {
    dmin = dist3(x1[i],y1[i],z1[i],x2[0],y2[0],z2[0]);
    jmin = 0;
    for (j=1; j<npt2; j++) {
      d = dist3(x1[i],y1[i],z1[i],x2[j],y2[j],z2[j]);
      if (d<dmin) {
	dmin = d;
	jmin = j;
      }
    }
    dmin = sqrt(dmin);
    printf("%d %g %g %g   %d %g %g %g   %g\n",i+1,x1[i],y1[i],z1[i],jmin+1,x2[jmin],y2[jmin],z2[jmin],dmin);
  }
}
