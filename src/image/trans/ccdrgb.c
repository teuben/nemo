/* 
 * CCDRGB: read an R, G and B image, and create composite
 *		** unfinished **
 *
 *	 8-jun-98 	written	- very quick and dirty       pjt
 *      10-jun-10       3rd image, and changed keyword names pjt
 *      13-jul-10       playing with OpenMP                  pjt
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in1=???\n      Input image in band 1",
        "in2=???\n      Input image in band 2",
        "in3=\n         Input image in band 3",
	"out=???\n      Output combined image",
	"c=0.5\n	color composition index RGB=(A, B+c*A, B)",
	"range1=\n      Range in image 1",
	"range2=\n      Range in image 2",
	"range3=\n      Range in image 3",
	"bw=\n          If given, scale factors for BW image creation (0.299,0.587,0.114)",
	"nbench=1\n     Benching OpenMP",
	"VERSION=0.3\n  13-jul-10 PJT",
	NULL,
};

string usage = "create composite RGB color images from multiple input images";

string cvsid="$Id$";


void nemo_main(void)
{
  stream   instr1, instr2, instr3, outstr;
  int      nx, ny, n1, n2, n3, i, j, nbw, nxy;
  imageptr iptr1=NULL, iptr2=NULL, iptr3=NULL, optr=NULL;      /* pointer to images */
  real     a, b, r,g, o, c = getdparam("c");
  real     range1[2], range2[2], range3[2], 
           scale1, offset1, scale2, offset2, scale3, offset3,
           *d1, *d2, *d3;
  real     bw[3], r_fac, g_fac, b_fac, gray;
  int      nbench = getiparam("nbench");
  
  instr1 = stropen(getparam("in1"), "r");
  n1 = nemoinpd(getparam("range1"),range1,2);
  read_image( instr1, &iptr1);
  if (n1==0) {
    scale1 = 1.0;
    offset1 = 0.0;
  } else if (n1==2) {
    scale1 = 1.0/(range1[1]-range1[0]);
    offset1 = range1[0];
  } else
    error("range1= needs two values");

  instr2 = stropen(getparam("in2"), "r");
  n2 = nemoinpd(getparam("range2"),range2,2);
  read_image( instr2, &iptr2);
  if (n2==0) {
    scale2 = 1.0;
    offset2 = 0.0;
  } else if (n2==2) {
    scale2 = 1.0/(range2[1]-range2[0]);
    offset2 = range2[0];
  } else
    error("range2= needs two values");
  
  
  if (hasvalue("in3")) {
    instr3 = stropen(getparam("in3"), "r");
    n3 = nemoinpd(getparam("range3"),range3,2);
    read_image( instr3, &iptr3);
    if (n3==0) {
      scale3 = 1.0;
      offset3 = 0.0;
    } else if (n3==2) {
      scale3 = 1.0/(range3[1]-range3[0]);
      offset3 = range3[0];
    } else
      error("range3= needs two values");
  } 

  if (hasvalue("bw")) {
    nbw = nemoinpd(getparam("bw"),bw,3);
  } else
    nbw = 0;
  

  outstr = stropen(getparam("out"), "w");

  nx = Nx(iptr1);
  ny = Ny(iptr1);
  if (iptr2) {
    if (nx!=Nx(iptr2)) error("Bad X-size of 2nd imagine",Nx(iptr2));
    if (ny!=Ny(iptr2)) error("Bad Y-size of 2nd imagine",Ny(iptr2));
  }
  if (iptr3) {
    if (nx!=Nx(iptr3)) error("Bad X-size of 3rd imagine",Nx(iptr3));
    if (ny!=Ny(iptr3)) error("Bad Y-size of 3rd imagine",Ny(iptr3));
  }

  if (nbw>0) {
    warning("OpenMP B/W image creation; nbench=%d",nbench);
    nxy = nx*ny;
    d1 = Frame(iptr1);
    d2 = Frame(iptr2);
    d3 = Frame(iptr3);
    r_fac = bw[0];
    g_fac = bw[1];
    b_fac = bw[2];
    for (j=0; j<nbench; j++) {
#pragma omp parallel for
      for (i=0; i<nxy; i++) {
	gray = r_fac*d1[i] + g_fac*d2[i] + b_fac*d3[i];
	d1[i] = gray;
      }
    }
      
  } else {
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	a = MapValue(iptr1,i,j);
	b = MapValue(iptr2,i,j);
	
	r = scale1*a + offset1;
	b = scale2*b + offset2;
	g = b+c*r;
	// (A, B+c*A, B)
	MapValue(iptr1,i,j) = g;
	// commonly used formulae for luminance
	// luminance =  0.3 * R + 0.59 * G + 0.11 * B
      }
    }
  }
  write_image(outstr, iptr1);
}


