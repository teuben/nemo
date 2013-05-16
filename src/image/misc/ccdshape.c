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
  "wcs=f\n        Use WCS of the cube (else use integer 0-based coordinates)",
  "VERSION=0.1\n  16-may-2013 PJT",
  NULL,
};

string usage = "shape of a 2D or 3D distribution based on moments of inertia";
string cvsid="$Id$";

vector oldframe[3] = {
    { 1.0, 0.0, 0.0, },
    { 0.0, 1.0, 0.0, },
    { 0.0, 0.0, 1.0, },
};



void nemo_main()
{
  stream  instr;
  string  oper;
  int     i,j,k,nx, ny, nz, nx1, ny1, nz1;
  int     mom;
  int     nclip, apeak, apeak1, cnt;
  imageptr iptr=NULL;             /* pointer to images */
  real    tmp0, tmp1, tmp2, tmp00, newvalue, peakvalue, scale, offset;
  real    *spec, cv, clip[2];
  bool    Qclip = hasvalue("clip");
  bool    Qwcs = getbparam("wcs");
  vector  tmpv, w_pos, pos, pos_b, frame[3];
  matrix  tmpm, w_qpole;
  real    w_sum;

  instr = stropen(getparam("in"), "r");
  if (Qclip) {
    nclip = nemoinpr(getparam("clip"),clip,2);
    if (nclip<1) error("error parsing clip=%s",getparam("clip"));
    if (nclip==1) {
      clip[1] =  clip[0];
      clip[0] = -clip[1];
    }
  }

  /* read the cube */

  read_image( instr, &iptr);
  nx = Nx(iptr);	
  ny = Ny(iptr);
  nz = Nz(iptr);
  
  /* loop over all relevant points and computer a rough center */

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
	cnt++;
	w_sum += cv;
	MULVS(tmpv, pos, cv);
	ADDV(w_pos, w_pos, tmpv);
      }
    }
  }
  DIVVS(w_pos,w_pos,w_sum);
  dprintf(0,"Size:   %d %d %d\n",nx,ny,nz);
  dprintf(0,"Center: %g %g %g (%d points)\n",w_pos[0], w_pos[1], w_pos[2],cnt);


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
  printvec("e_x:", frame[0]);
  printvec("e_y:", frame[1]);
  printvec("e_z:", frame[2]);
}


#include "nrutil.h"

eigenframe(vector frame[], matrix mat)
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

printvec(string name, vector vec, stream out)
{
  vector rtp;	/* radius - theta - phi */
  xyz2rtp(vec,rtp);
  printf("%12s  %10.5f  %10.5f  %10.5f  %10.5f   %5.1f %6.1f\n",
	  name, rtp[0], vec[0], vec[1], vec[2],
	  rtp[1]*180.0/PI, rtp[2]*180.0/PI);
}


xyz2rtp(vector xyz, vector rtp)
{
  real z = xyz[2];
  real w = sqrt(sqr(xyz[0])+sqr(xyz[1]));
  rtp[1] = atan(w/z);                 /* theta: in range 0 .. PI */
  if (z<0) rtp[1] += PI;
  rtp[2] = atan2(xyz[1], xyz[0]);     /* phi: in range  -PI .. PI */
  rtp[0] = sqrt(z*z+w*w);
}
