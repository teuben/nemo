/* 
 * CCDPOT: potential of an infinitesimally thin disk
 *
 *	26-jul-02   q&d, from Gipsy's potential.dc1  (the slow coffee way)  pjt
 *                  dumb coding:  128*128 takes 47.8" on a P600 (pjt's laptop)
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in=???\n       Input image file",
	"out=???\n      Output image file",
	"VERSION=0.1\n  26-jul-02 PJT",
	NULL,
};

string usage = "potential of a density disk";



#define CVI(x,y)  MapValue(iptr,x,y)
#define CVO(x,y)  MapValue(optr,x,y)
#define DIS(x,y)  MapValue(dptr,x,y)

#define QABS(a,b) (a>b ? a-b : b-a)

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny;
    int     i,j,k,l,k1,l1;
    real    sum, val, add, d;
    imageptr iptr=NULL, optr, dptr;

    instr = stropen(getparam("in"), "r");
    read_image( instr, &iptr);
    nx = Nx(iptr);	
    ny = Ny(iptr);

    outstr = stropen(getparam("out"), "w");
    create_image(&optr,nx,ny);
    Dx(optr) = Dx(iptr);
    Dy(optr) = Dy(iptr);
    Xmin(optr) = Xmin(iptr);
    Ymin(optr) = Ymin(iptr);

#if 0
    create_image(&dptr,nx,ny);
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++)
	if (i>0 || j>0) 
	  DIS(i,j) = 1.0/sqrt((double)(i**2 + j**2));
	else
	  DIS(0,0) = 3.54;
#endif


    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	sum = 0.0;
	for (l=0; l<ny; l++) {
	  l1 = QABS(j,l);
	  for (k=0; k<nx; k++) {
	    k1 = QABS(k,i);
	    val = CVI(k,l);
	    if (k==i && l==j)
	      add = 3.54 * val;
	    else {
	      d = (i-k)*(i-k) + (j-l)*(j-l);
	      add = val / sqrt(d);
	    }
	    sum += add;
	  }
	}
	CVO(i,j) = sum;
      }
    } 
    write_image(outstr, optr);
}
