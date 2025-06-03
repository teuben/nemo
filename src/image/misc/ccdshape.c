/* 
 * CCDSHAPE: shape of a 2D or 3D distribution based on moments of inertia
 *
 *      (based off snapkinem/snaprect)
 *	quick and dirty:  16-may-2013		pjt
 */



#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image file",
  "clip=\n        Use only values above clip",
  "wcs=t\n        Use WCS of the cube (else use integer 0-based coordinates)",
  "radecvel=f\n   Split the RA/DEC from VEL",
  "weight=t\n     Weights by intensity",
  "cross=t\n      Use cross correlations between X and Y to get angles",
  "scale=1\n      Scale factor to be applied to radii",
  "VERSION=0.6n   29-may-2025 PJT",
  NULL,
};

string usage = "shape of a 2D or 3D distribution based on moments of inertia";

vector oldframe[3] = {
    { 1.0, 0.0, 0.0, },
    { 0.0, 1.0, 0.0, },
    { 0.0, 0.0, 1.0, },
};

real printeig(string name, matrix mat, real *a, real *b, real *c);
real printvec(string name, vector vec);
void xyz2rtp(vector xyz, vector rtp);
void eigsrt(float *d,float **v, int n);
void eigenframe(vector frame[], matrix mat);
void jacobi(float **a,int n,float *d,float **v,int *nrot);


void nemo_main()
{
  stream  instr;
  string  oper;
  int     i,j,k,nx, ny, nz, nx1, ny1, nz1, mom;
  int     nclip, apeak, apeak1, cnt;
  imageptr iptr=NULL;             /* pointer to image */
  real    tmp0, tmp1, tmp2, tmp00, newvalue, peakvalue, offset;
  real    *spec, cv, clip[2];
  real    scale = getdparam("scale");  
  bool    Qclip = hasvalue("clip");
  bool    Qwcs = getbparam("wcs");
  bool    Qrdv = getbparam("radecvel");
  bool    Qiwm = getbparam("weight");
  bool    Qcross = getbparam("cross");
  vector  tmpv, w_pos, pos, pos_b, ds, frame[3];
  matrix  tmpm, w_qpole;
  real    w_sum, dmin, dmax;
  real    inc, pa_k, pa_m, dPA, a_m, b_m, c_m, a_k, b_k, c_k, dvdr;

  instr = stropen(getparam("in"), "r");
  if (Qclip) {
    nclip = nemoinpr(getparam("clip"),clip,2);
    if (nclip<1) error("error parsing clip=%s",getparam("clip"));
    if (nclip==1) {
      clip[1] =  clip[0];
      clip[0] = -clip[1];
    }
  }
  read_image( instr, &iptr);            /* read the cube */
  nx = Nx(iptr);  ny = Ny(iptr);  nz = Nz(iptr);
  
  /* loop over all relevant points and compute a rough center */
  cnt = 0;
  w_sum = 0.0;
  CLRV(w_pos);
  for (k=0; k<nz; k++) {
    pos[2] = Qwcs ? k*Dz(iptr) + Zmin(iptr)  :  k;
    for (j=0; j<ny; j++) {
      pos[1] = Qwcs ? j*Dy(iptr) + Ymin(iptr)  :  j;
      for (i=0; i<nx; i++) {
	pos[0] = Qwcs ? i*Dx(iptr) + Xmin(iptr)  :  i;
	cv = CubeValue(iptr,i,j,k);
	if (Qclip && (clip[0]<=cv && cv<=clip[1])) continue;
	if (cnt==0) {
	  dmin = dmax = cv;
	} else {
	  dmin = MIN(dmin, cv);
	  dmax = MAX(dmax, cv);
	}
	cnt++;
	if (!Qiwm) cv = 1.0;
	w_sum += cv;
	MULVS(tmpv, pos, cv);
	ADDV(w_pos, w_pos, tmpv);
      }
    }
  }
  DIVVS(w_pos,w_pos,w_sum);
  printf("Npoints:    %d\n",cnt);
  printf("DataMinMax: %g %g\n",dmin,dmax);
  printf("DataSum:    %g\n",w_sum);
  printf("Size:       %d %d %d\n",nx,ny,nz);
  printf("Center:     %g %g %g %s\n",w_pos[0], w_pos[1], w_pos[2],
	 Qwcs ? "[wcs]" : "[grid]");

  /* based on this center, compute quadrupole moments */

  CLRM(w_qpole);
  for (k=0; k<nz; k++) {
    pos[2] = Qwcs ? k*Dz(iptr) + Zmin(iptr)  :  k;
    for (j=0; j<ny; j++) {
      pos[1] = Qwcs ? j*Dy(iptr) + Ymin(iptr)  :  j;
      for (i=0; i<nx; i++) {
	pos[0] = Qwcs ? i*Dx(iptr) + Xmin(iptr)  :  i;
	cv = CubeValue(iptr,i,j,k);
	if (Qclip && (clip[0]<=cv && cv<=clip[1])) continue;
	cnt++;
	if (!Qiwm) cv = 1.0;
	SUBV(pos_b, pos, w_pos);
	MULVS(tmpv, pos_b, cv);
	OUTVP(tmpm, tmpv, pos_b);
	ADDM(w_qpole, w_qpole, tmpm);
      }
    }
  }
  DIVMS(w_qpole, w_qpole, w_sum);
  if (!Qcross) {
    w_qpole[0][1] = w_qpole[0][2] = 0.0;
    w_qpole[1][0] = w_qpole[1][2] = 0.0;
    w_qpole[2][0] = w_qpole[2][1] = 0.0;
  }
  
  /* get the meat */

  eigenframe(frame, w_qpole);
  if (dotvp(oldframe[0], frame[0]) < 0.0)
    MULVS(frame[0], frame[0], -1.0);
  if (dotvp(oldframe[2], frame[2]) < 0.0)
    MULVS(frame[2], frame[2], -1.0);
  CROSSVP(frame[1], frame[2], frame[0]);
  pa_k = printvec("e_x:", frame[0]) + 90.0;
  printvec("e_y:", frame[1]);
  printvec("e_z:", frame[2]);
  inc = printeig("qpole:",w_qpole, &a_k, &b_k, &c_k);
  if (Qcross) {
    printf("a,b:  %g %g\n",a_k*scale,b_k*scale);
    printf("inc:  %g (meaningless without radecvel)\n",inc);
    printf("pa:   %g (meaningless without radecvel)\n",pa_k);
  } else {
    printf("x,y,z: %g %g %g\n", sqrt(w_qpole[0][0]), 
	   sqrt(w_qpole[1][1]), sqrt(w_qpole[2][2]));
  }
  

  if (Qrdv) {
    warning("RA-DEC-VEL cube assumed. Now presenting decoupled geometry");

    w_pos[2] = 0.0;
    
    CLRM(w_qpole);
    for (k=0; k<nz; k++) {
      pos[2] = 0.0;
      for (j=0; j<ny; j++) {
	pos[1] = Qwcs ? j*Dy(iptr) + Ymin(iptr)  :  j;
	for (i=0; i<nx; i++) {
	  pos[0] = Qwcs ? i*Dx(iptr) + Xmin(iptr)  :  i;
	  cv = CubeValue(iptr,i,j,k);
	  if (Qclip && (clip[0]<=cv && cv<=clip[1])) continue;
	  cnt++;
	  if (!Qiwm) cv = 1.0;
	  SUBV(pos_b, pos, w_pos);
	  MULVS(tmpv, pos_b, cv);
	  OUTVP(tmpm, tmpv, pos_b);
	  ADDM(w_qpole, w_qpole, tmpm);
	}
      }
    }
    DIVMS(w_qpole, w_qpole, w_sum);
  
    /* get the meat */

    eigenframe(frame, w_qpole);
    if (dotvp(oldframe[0], frame[0]) < 0.0)
      MULVS(frame[0], frame[0], -1.0);
    if (dotvp(oldframe[2], frame[2]) < 0.0)
      MULVS(frame[2], frame[2], -1.0);
    CROSSVP(frame[1], frame[2], frame[0]);
    pa_m = printvec("e_x:", frame[0]) + 90.0;
    printvec("e_y:", frame[1]);
    printvec("e_z:", frame[2]);
    inc = printeig("qpole:",w_qpole, &a_m, &b_m, &c_m);

    printf("a,b:   %g %g\n",a_m*scale,b_m*scale);
    printf("inc:   %g\n",inc);
    printf("pa_m:  %g\n",pa_m);
    printf("pa_k:  %g\n",pa_k);
    dPA = ABS(pa_k-pa_m);
    if (dPA > 90) dPA = 180-dPA;
    printf("dPA:   %g\n",dPA);
    printf("dV/dR: %g\n",sqrt(sqr(a_k)-sqr(a_m))/a_m);
    dvdr = sqrt(sqr(a_k)-sqr(a_m))/a_m/sin(inc*PI/180);
    printf("dV/dR: %g sini corrected\n",dvdr);
    /* dV/dR includes sin(i) for rotation */
    /* dV/dR = {O_r*cos(dPA) ,O_e*sin(dPA)} */
   }
}


/* see also snaprect/snapkinem */

#include "nrutil.h"

void eigenframe(vector frame[], matrix mat)
{
    float **q, *d, **v;
    int i, j, nrot;

    q = fmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++)
	for (j = 1; j <= 3; j++)
	    q[i][j] = mat[i-1][j-1];
    d = fvector(1, 3);
    v = fmatrix(1, 3, 1, 3);
    jacobi(q, 3, d, v, &nrot);
    eigsrt(d, v, 3);
    for (i = 1; i <= 3; i++)
	for (j = 1; j <= 3; j++)
	    frame[i-1][j-1] = v[j][i];
}

real printeig(string name, matrix mat,  real *a, real *b, real *c)
{
    float **q, *d, **v;
    int i, j, nrot;
    real inc;

    q = fmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++)
	for (j = 1; j <= 3; j++)
	    q[i][j] = mat[i-1][j-1];
    d = fvector(1, 3);
    v = fmatrix(1, 3, 1, 3);
    jacobi(q, 3, d, v, &nrot);
    eigsrt(d, v, 3);
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n", name,
	   d[1], v[1][1], v[2][1], v[3][1]);
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n", "            ",
	   d[2], v[1][2], v[2][2], v[3][2]);
    printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f\n", "            ",
	   d[3], v[1][3], v[2][3], v[3][3]);

    inc = acos(sqrt(d[2]/d[1]))*180.0/PI;
    *a = sqrt(d[1]);
    *b = sqrt(d[2]);
    *c = sqrt(d[3]);
    return inc;
}




real printvec(string name, vector vec)
{
  vector rtp;	/* radius - theta - phi */
  real pa;
  xyz2rtp(vec,rtp);
  printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f   %5.1f %6.1f\n",
	  name, rtp[0], vec[0], vec[1], vec[2],
	  rtp[1]*180.0/PI, rtp[2]*180.0/PI);
  pa = rtp[2]*180.0/PI;
  return pa;
}


void xyz2rtp(vector xyz, vector rtp)
{
  real z = xyz[2];
  real w = sqrt(sqr(xyz[0])+sqr(xyz[1]));
  rtp[1] = atan(w/z);                 /* theta: in range 0 .. PI */
  if (z<0) rtp[1] += PI;
  rtp[2] = atan2(xyz[1], xyz[0]);     /* phi: in range  -PI .. PI */
  rtp[0] = sqrt(z*z+w*w);
}

