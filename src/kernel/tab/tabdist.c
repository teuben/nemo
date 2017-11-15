/*
 *  Compare distance between points in two tables
 *
 *      13-feb-2013     0.1    Created, Q&D
 *      24-jun-2016     0.7    allow radec in degrees; fix C language issues, introduce dmin=,
 *                             faster precompute
 *
 *  See also the SIRTF/map2 project in 'c2d'
 *
 *  @todo   consider allowing scaling/offset between two datasets 
 *          Y = aX + b
 *          use k-d Match?   http://arxiv.org/abs/1304.0838
 *  @todo   ID column
 *
 *  @todo   Faster lookup using 2D grid with linked list?
 *
 *  @todo   Generalize to N tables?
 */

#include <stdinc.h>
#include <getparam.h>
#include <grid.h>

string defv[] = {
    "in1=???\n		First (ascii table) dataset",
    "in2=???\n		Second (ascii table) dataset",
    "col1=1,2\n         X (optionally Y and Z) Column from 1st dataset",
    "col2=1,2\n         X (optionally Y and Z) Column from 2nd dataset",
    "id1=0\n            Column in first data reprenting the ID",
    "id2=0\n            Column in second data reprenting the ID",
    "radec=false\n      Interpret X and Y as RA/DEC angles on sky",
    "units=r\n          radians or degree",
    "xgrid=\n           Grid to limit nearest neighbors for special 2D case",
    "ygrid=\n           Grid to limit nearest neighbors for special 2D case",
    "dmin=-1\n          If positive, only print distances less than this value",
    "VERSION=0.7\n	24-jun-2016 PJT",
    NULL,
};

string usage="Compare distances between points in two tables";
string cvsid="$Id$";

#define MAXGRID  1024

void compare_1d (int npt1, int npt2, real *x1, real *x2, real dmin);
void compare_2d (int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, real dmin);
void compare_2da(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, real scale, string *name1, string *name2, real dmin);
void compare_2dg(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, Grid *g, real dmin);
void compare_3d (int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, real *z1, real *z2, real dmin);

void get_atable_a(stream instr, int id, int npt, string *name);

inline static real dist1(real x1,real x2) {
  real d = x1-x2;
  if (d > 0.0) return d;
  return -d;
}

inline static real dist2(real x1,real y1,real x2,real y2) {
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

inline static real dist3(real x1,real y1,real z1,real x2,real y2,real z2) {
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}

/* slow version : compute sin() and cos() in the fly */

inline static real dist2a(real x1,real y1,real x2,real y2) {
  return sin(y1)*sin(y2)+cos(y1)*cos(y2)*cos(x1-x2);
}

/* fast version : stores and precomputes */
   
static real *X1, *X2;
static real *Sinx1, *Sinx2, *Siny1, *Siny2;
static real *Cosx1, *Cosx2, *Cosy1, *Cosy2;

inline static real dist2ai(i,j) {
  return Siny1[i]*Siny2[j]+Cosy1[i]*Cosy2[j]*cos(X1[i]-X2[j]);
}



void key2grid(string key, Grid *g)
{

}


#define NDIM   3

nemo_main()
{
  int i, npt, npt1, npt2, ncol1, ncol2, id1, id2;
  real *coldat1[NDIM], *coldat2[NDIM];
  int colnr1[NDIM], colnr2[NDIM];
  real *x1 = NULL, *x2 = NULL;
  real *y1 = NULL, *y2 = NULL;
  real *z1 = NULL, *z2 = NULL;
  real ascale = 1.0;
  real dmin = getrparam("dmin");
  string input1 = getparam("in1");
  string input2 = getparam("in2");
  stream instr1, instr2;
  string *name1, *name2;
  bool Qradec = getbparam("radec");
  string units = getparam("units");
  Grid gx, gy;

  if (Qradec) warning("new radec mode ");
  if (*units == 'd') ascale = PI/180.0;

  ncol1 = nemoinpi(getparam("col1"),colnr1,NDIM);
  ncol2 = nemoinpi(getparam("col2"),colnr2,NDIM);

  key2grid("xgrid",&gx);
  key2grid("ygrid",&gy);

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

  if (Qradec) {   /* fast precompute for angular case */
    warning("precompute, ascale=%g",ascale);
    X1 = (real *) allocate(npt1 * sizeof(real));
    X2 = (real *) allocate(npt2 * sizeof(real));
    Sinx1 = (real *) allocate(npt1 * sizeof(real));
    Siny1 = (real *) allocate(npt1 * sizeof(real));
    Cosx1 = (real *) allocate(npt1 * sizeof(real));
    Cosy1 = (real *) allocate(npt1 * sizeof(real));
    Sinx2 = (real *) allocate(npt2 * sizeof(real));
    Siny2 = (real *) allocate(npt2 * sizeof(real));
    Cosx2 = (real *) allocate(npt2 * sizeof(real));
    Cosy2 = (real *) allocate(npt2 * sizeof(real));
    for (i=0; i<npt1; i++) {
      X1[i]    = x1[i]*ascale;
      Sinx1[i] = sin(x1[i]*ascale);
      Siny1[i] = sin(y1[i]*ascale);
      Cosx1[i] = cos(x1[i]*ascale);
      Cosy1[i] = cos(y1[i]*ascale);
    }
    for (i=0; i<npt2; i++) {
      X2[i]    = x2[i]*ascale;
      Sinx2[i] = sin(x2[i]*ascale);
      Siny2[i] = sin(y2[i]*ascale);
      Cosx2[i] = cos(x2[i]*ascale);
      Cosy2[i] = cos(y2[i]*ascale);
    }
  }

  id1 = getiparam("id1");
  id2 = getiparam("id2");
  if (id1 > 0) {
    rewind(instr1);
    name1 = (string *) allocate(npt1 * sizeof(string));
    get_atable_a(instr1,id1,npt1,name1);
  } else
    name1 = NULL;
  if (id2 > 0) {
    rewind(instr2);
    name2 = (string *) allocate(npt2 * sizeof(string));
    get_atable_a(instr2,id2,npt2,name2);
  } else
    name2 = NULL;
  
  strclose(instr1);
  strclose(instr2);

  if (Qradec && ncol1==2)
    compare_2da(npt1, npt2, x1, x2, y1, y2, ascale, name1, name2, dmin);
  else if (ncol1 == 1) 
    compare_1d(npt1, npt2, x1, x2, dmin);
  else if (ncol1 == 2) 
    compare_2d(npt1, npt2, x1, x2, y1, y2, dmin);
  else if (ncol1 == 3)
    compare_3d(npt1, npt2, x1, x2, y1, y2, z1, z2, dmin);
  else
    warning("%d columns not supported",ncol1);
}




void compare_1d(int npt1, int npt2, real *x1, real *x2, real dmin_out)
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
    if (dmin_out>0 && dmin<dmin_out)
      printf("%d %g   %d %g   %g\n",i+1,x1[i],jmin+1,x2[jmin], dmin);
  }
}

void compare_2d(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, real dmin_out)
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
    if (dmin_out>0 && dmin<dmin_out)
      printf("%d %g %g   %d %g %g   %g\n",i+1,x1[i],y1[i],jmin+1,x2[jmin], y2[jmin],dmin);
  }
}

/* fast version - doesn't work yet 
 * Npt1=1493 Npt2=12681
 * cartesian:  0.2"
 * slow:       8.0"
 * fast:       2.1"
 */

/* fast */
void compare_2da(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, real s, string *name1, string *name2, real dmin_out)
{  
  int i, j, jmin;
  real d, dmax, dmin;
  string n1, n2;
  char ord1[16], ord2[16];

  
  for (i=0; i<npt1; i++) {
    if (name1) dprintf(2,"NAME1: %s\n",name1[i]);
    dmax = dist2ai(i,0);
    jmin = 0;
    for (j=1; j<npt2; j++) {
      d = dist2ai(i,j);
      if (d>dmax) {
	dmax = d;
	jmin = j;
      }
    }
    dmin = acos(dmax)*DR2AS;  /* radians to arcsec */

    if (name1==NULL) {
      sprintf(ord1,"%d",i+1);
      n1 = ord1;
    } else
      n1 = name1[i];

    if (name2==NULL) {
      sprintf(ord2,"%d",jmin+1);
      n2 = ord2;
    } else
      n2 = name2[jmin];

    if (dmin_out>0 && dmin<dmin_out)
      printf("%s %g %g   %s %g %g   %g\n",n1,x1[i],y1[i],n2,x2[jmin],y2[jmin],dmin);
  }
}

/* slow */
void compare_2da_slow(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, real s, string *name1, string *name2, real dmin_out)
{  
  int i, j, jmin;
  real d, dmax, dmin;
  string n1, n2;
  char ord1[16], ord2[16];

  
  for (i=0; i<npt1; i++) {
    if (name1) dprintf(2,"NAME1: %s\n",name1[i]);
    dmax = dist2a(x1[i]*s,y1[i]*s,x2[0]*s,y2[0]*s);
    jmin = 0;
    for (j=1; j<npt2; j++) {
      d = dist2a(x1[i]*s,y1[i]*s,x2[j]*s,y2[j]*s);
      if (d>dmax) {
	dmax = d;
	jmin = j;
      }
    }
    dmin = acos(dmax)*DR2AS;  /* radians to arcsec */

    if (name1==NULL) {
      sprintf(ord1,"%d",i+1);
      n1 = ord1;
    } else
      n1 = name1[i];

    if (name2==NULL) {
      sprintf(ord2,"%d",jmin+1);
      n2 = ord2;
    } else
      n2 = name2[jmin];

    if (dmin_out>0 && dmin<dmin_out)
      printf("%s %g %g   %s %g %g   %g\n",n1,x1[i],y1[i],n2,x2[jmin],y2[jmin],dmin);
  }
}

void compare_2g(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, Grid *gx, Grid *gy, real dmin_out)
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
    if (dmin_out>0 && dmin<dmin_out)
      printf("%d %g %g   %d %g %g   %g\n",i+1,x1[i],y1[i],jmin+1,x2[jmin], y2[jmin],dmin);
  }
}


void compare_3d(int npt1, int npt2, real *x1, real *x2, real *y1, real *y2, real *z1, real *z2, real dmin_out)
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

/* 
 *  read lines and pick a column as the ID.
 *  skip lines that start with '#'
 */

void get_atable_a(stream instr, int id, int npt, string *name)
{
  static char line[MAX_LINELEN];
  string *sp;
  int i;

  for (i=0; i<npt; ) {
    if (fgets(line,MAX_LINELEN,instr) == NULL)
      error("unexpected EOF on table read (%d,%d)",id,npt);
    if (line[0] == '#') continue;
    sp = burststring(line,", \t\r");
    if (xstrlen(sp,sizeof(string))<id-1) {
      warning("skipping line, not enough columns");
      continue;
    }
    name[i] = strdup(sp[id-1]);
    freestrings(sp);
    i++;
  }
  dprintf(1,"First and last ID: %s %s\n",name[0],name[npt-1]);
}
