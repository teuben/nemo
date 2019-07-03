/* 
 * TESTBSC: test Beam Smearing Correction code
 *
 *      2-jul-2019    Begeman's appendix B
 *
 *
 * We will be using our (D,V,S) and (d,v,s) notation, not Begeman's (N,V,S) and (n,w,.)
 *
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "inD=\n        Input Model Density image file",
        "inV=\n        Input Model Velocity image file",
        "ind=\n        Input Observed Density image file",
        "inv=\n        Input Observed Velocity image file",
	"out=\n        Output (TBD)",
	"beam=\n       Override beam/pixel ratio from header value",
	"diff=f\n      Difference Map? (as opposed to Observed)",
	"beeg=t\n      Use the Begeman correction terms to V",
	"VERSION=0.4\n 3-jul-2019 PJT",
	NULL,
};

string usage = "testing Beam Smearing Corrections";

#define CV01(x,y)  MapValue(iptr01,x,y)
#define CV02(x,y)  MapValue(iptr02,x,y)
#define CV11(x,y)  MapValue(iptr11,x,y)
#define CV12(x,y)  MapValue(iptr12,x,y)
#define CVO(x,y)   MapValue(optr,x,y)


void nemo_main()
{
    stream   instr01,instr02,instr11,instr12,outstr;
    int      nx, ny;
    int      i,j;
    imageptr iptr01=NULL, iptr02=NULL;      /* pointer to images */
    imageptr iptr11=NULL, iptr12=NULL;
    imageptr optr = NULL;
    imageptr iptr = NULL;
    real     d1, d2, d3, d4, dx, dy, b;
    real     v1, n1, v2, n2;
    bool     Qdiff = getbparam("diff");
    bool     Qbeeg = getbparam("beeg");

    if (hasvalue("inD")) {
      instr01 = stropen(getparam("inD"),"r");
      read_image( instr01, &iptr01);
      iptr = iptr01;
    }
    if (hasvalue("inV")) {
      instr02 = stropen(getparam("inV"),"r");
      read_image( instr02, &iptr02);
      iptr = iptr02;      
    }
    if (hasvalue("ind")) {
      instr11 = stropen(getparam("ind"),"r");
      read_image( instr11, &iptr11);
      iptr = iptr11;      
    }
    if (hasvalue("inv")) {
      instr12 = stropen(getparam("inv"),"r");
      read_image( instr12, &iptr12);
      iptr = iptr12;      
    }
    if (iptr == NULL) error("No input image(s) given");
    
    nx = Nx(iptr);	
    ny = Ny(iptr);
    dx = Dx(iptr);
    dy = Dy(iptr);
    /* @todo    should check if maps are same size */

    b = Beamy(iptr) / Dy(iptr) / 2.355;    // 2.355 = 2*sqrt(2*ln(2))
    if (hasvalue("beam"))
      b = getrparam("beam");

    dprintf(0,"Maps size %d x %d assumed,   beam/pixel = %g\n",nx,ny, b);

    outstr = stropen(getparam("out"), "w");
    create_cube(&optr,nx,ny,1);
    Dx(optr) = Dx(iptr);
    Dy(optr) = Dy(iptr);
    Xmin(optr) = Xmin(iptr);
    Ymin(optr) = Ymin(iptr);
    /* should do the others too */

    for (j=1; j<ny-1; j++) {
      for (i=1; i<nx-1; i++) {
	if (CV01(i,j) == 0) {
	  CVO(i,j) = 0;        // blank out where model is 0
	  continue;
	}
	// 2nd order terms (Laplacian)
	d1 = CV02(i-1,j) - CV02(i,j);
	d2 = CV02(i+1,j) - CV02(i,j); 
	d3 = CV02(i,j-1) - CV02(i,j);
	d4 = CV02(i,j+1) - CV02(i,j); 
	// linear terms
	n1 = 0.5*(CV01(i+1,j) - CV01(i-1,j));    // dN/dx
	v1 = 0.5*(CV02(i+1,j) - CV02(i-1,j));    // dV/dx
	n2 = 0.5*(CV01(i,j+1) - CV01(i,j-1));    // dN/dy
	v2 = 0.5*(CV02(i,j+1) - CV02(i,j-1));    // dV/dy
	//
	CVO(i,j)  = CV02(i,j);
	if (Qbeeg) {
	  CVO(i,j) += (v1*n1 + v2*n2) * b * b / CV11(i,j);                   // linear
	  CVO(i,j) += (d1+d2+d3+d4) * CV01(i,j) * b * b / ( 2 * CV11(i,j));  // 2nd order
	}
	if (Qdiff)
	  CVO(i,j) -= CV12(i,j);
      }
    }
    write_image(outstr, optr);
}
