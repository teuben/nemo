/* 
 * TESTBSC: test Beam Smearing Correction code
 *
 *      2-jul-2019  Begeman's appendix B using our notation
 *
 *
 * We will be using our (D,V,S) and (d,v,s) notation, not Begeman's (N,V,S) and (n,w,.)
 *
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
	"beam=\n       Override beam/pixel ratio",
	"diff=f\n      Difference Map? (as opposed to Observed)",
	"beeg=t\n      Use the Begeman correction terms to V",
	"VERSION=0.2\n 2-jul-2019 PJT",
	NULL,
};

string usage = "testing Beam Smearing Corrections";



#define CV01(x,y)  MapValue(iptr01,x,y)
#define CV02(x,y)  MapValue(iptr02,x,y)
#define CV11(x,y)  MapValue(iptr11,x,y)
#define CV12(x,y)  MapValue(iptr12,x,y)
#define CVO(x,y)   MapValue(optr,x,y)

local string valid_modes = "laplace,lapabs,aregan,pregan,divergence,vorticity";

#define MODE_LAPLACE (1<<0)
#define MODE_LAPABS  (1<<1)
#define MODE_AREGAN  (1<<2)
#define MODE_PREGAN  (1<<3)
#define MODE_DIV     (1<<4)
#define MODE_VORT    (1<<5)

extern int match(string, string, int *);

void nemo_main()
{
    stream  instr01,instr02,instr11,instr12,outstr;
    int     nx, ny, nz, mode;
    int     i,j,k;
    imageptr iptr01=NULL, iptr02=NULL;      /* pointer to images */
    imageptr iptr11=NULL, iptr12=NULL;
    imageptr optr = NULL;
    imageptr iptr = NULL;
    real    d1, d2, d3, d4, dx, dy, b;
    real    v1, n1, v2, n2;
    bool    Qsym = TRUE;            /* symmetric derivates w.r.t. pixel point */
    bool    Qdiff = getbparam("diff");
    bool    Qbeeg = getbparam("beeg");

#if 0    
    match(getparam("mode"),valid_modes,&mode);
    if (mode==0) error("Not a valid mode; valid:%s",valid_modes);
    dprintf(0,"Image sharpening method #%d\n",mode);
#endif    

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
	
	// 2nd order terms
	d1 = CV02(i-1,j) - CV02(i,j);
	d2 = CV02(i+1,j) - CV02(i,j); // d2V/dx2 = d1+d2
	d3 = CV02(i,j-1) - CV02(i,j);
	d4 = CV02(i,j+1) - CV02(i,j); // d2V/dy2 = d3+d4
	// linear terms
	n1 = 0.5*(CV01(i+1,j) - CV01(i-1,j));    // dN/dx
	v1 = 0.5*(CV02(i+1,j) - CV02(i-1,j));    // dV/dx
	n2 = 0.5*(CV01(i,j+1) - CV01(i,j-1));    // dN/dy
	v2 = 0.5*(CV02(i,j+1) - CV02(i,j-1));    // dV/dy
	//
	CVO(i,j)  = CV02(i,j);
	if (Qbeeg) {
	  CVO(i,j) += (d1+d2+d3+d4) * CV01(i,j) * b * b / ( 2 * CV11(i,j));
	  CVO(i,j) += (v1*n1 + v2*n2) * b * b / CV11(i,j);
	}
	if (Qdiff)
	  CVO(i,j) -= CV12(i,j);
      }
    }

    write_image(outstr, optr);
    


#if 0
    } else if (mode & MODE_DIV || mode & MODE_VORT) {
        dprintf(1,"mode=div/vort\n");
        for (k=0; k<nz; k++) {
            for (j=0; j<ny-1; j++) {
                for (i=0; i<nx-1; i++) {
                    if (Qsym) {
                        if (i>0 && j>0) {
                            d1 = 0.5*(CV1(i+1,j,k) - CV1(i-1,j,k));         /* dv_x/dx */
                            d2 = 0.5*(CV1(i,j+1,k) - CV1(i,j-1,k));         /* dv_x/dy */
                            d3 = 0.5*(CV2(i+1,j,k) - CV2(i-1,j,k));         /* dv_y/dx */
                            d4 = 0.5*(CV2(i,j+1,k) - CV2(i,j-1,k));         /* dv_y/dy */
                        } else
                            d1 = d2 = d3 = d4 = 0.0;
                    } else {
                        d1 = CV1(i+1,j,k) - CV1(i,j,k);         /* dv_x/dx */
                        d2 = CV1(i,j+1,k) - CV1(i,j,k);         /* dv_x/dy */
                        d3 = CV2(i+1,j,k) - CV2(i,j,k);         /* dv_y/dx */
                        d4 = CV2(i,j+1,k) - CV2(i,j,k);         /* dv_y/dy */
                    }
                    if (mode&MODE_DIV)
                        CVO(i,j,k) = d1/dx + d4/dy;
                    else if (mode&MODE_VORT)
                        CVO(i,j,k) = d3/dx - d2/dy;
                }
                CVO(nx-1,j,k) = 0.0;
            }
            for (i=0; i<nx; i++) {
                CVO(i,ny-1,k) = 0.0;
            }
        }
    } else if (mode & MODE_AREGAN || mode & MODE_PREGAN) {
        dprintf(1,"mode=aregan/pregan\n");      
        for (k=0; k<nz; k++) {
            for (j=0; j<ny-1; j++) {
                for (i=0; i<nx-1; i++) {
                    d1 = CV1(i,j,k)   - CV1(i+1,j,k);
                    d2 = CV1(i,j+1,k) - CV1(i+1,j+1,k);
                    d3 = CV1(i,j,k)   - CV1(i,j+1,k);
                    d4 = CV1(i+1,j,k) - CV1(i+1,j+1,k);
                    if (mode&MODE_AREGAN)
                        CVO(i,j,k) = sqrt(sqr(d1+d2)+sqr(d3+d4))/2;
                    else {
                        if (d3+d4==0.0 && d1+d2==0.0)
                            CVO(i,j,k) = 0.0;
                        else
                            CVO(i,j,k) = atan2(d3+d4,d1+d2) * 180 / PI;
                    }
                }
                CVO(nx-1,j,k) = 0.0;
            }
            for (i=0; i<nx; i++) {
                CVO(i,ny-1,k) = 0.0;
            }
        }
    }
    
#endif
}
