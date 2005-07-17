/* 
 * CCDPOT: potential of an infinitesimally thin disk - the slow way 
 *
 *	26-jul-02   q&d, from Gipsy's potential.dc1  (the slow coffee way)  pjt
 *      28-feb-03   changed sign to make potentials most negative in center, use G
 *      13-feb-05   0.5 nbench=
 *
 *                  dumb coding:  128*128 takes 47.8" on a P600 (pjt's laptop)
 *                  using dptr    128*128 takes  9.7" (speedup 5)
 *        On P1.6:  
 *          8^2         0.000505 *
 *         16^2         0.00337  *
 *         32^2         0.0479   *
 *         64^2         0.742    *
 *        128^2 map:    4.0"     11.76
 *        256^2 map:  410.2"     579.
 *        512^2 map:
 */
  

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in=???\n       Input image file",
	"out=???\n      Output image file",
	"gravc=1\n      Gravitational Constant",
	"report=0\n     report if this number cells done (0=none)",
	"nbench=1\n     benchmark number for N^4 convolution",
	"VERSION=0.5\n  13-apr-05 PJT",
	NULL,
};

string usage = "potential of a thin disk - the slow way";



#define CVI(x,y)  MapValue(iptr,x,y)
#define CVO(x,y)  MapValue(optr,x,y)
#define DIS(x,y)  MapValue(dptr,x,y)

#define QABS(a,b) (a>b ? a-b : b-a)

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny;
    int     i,j,k,l,k1,l1;
    real    sum, val, add, d, dx,dy,gravc = getdparam("gravc");
    imageptr iptr=NULL, optr, dptr;
    int     report = getiparam("report");
    int     nbench = getiparam("nbench");
    int     count = 0;

    if (nbench < 1) error("Bad value nbench=%d",nbench);

    /* read input image */

    instr = stropen(getparam("in"), "r");          
    read_image( instr, &iptr);                     
    nx = Nx(iptr);	
    ny = Ny(iptr);
    dx = ABS(Dx(iptr));
    dy = ABS(Dy(iptr));
    gravc *= sqrt(dx*dy);
    if (dx != dy) 
      warning("Pixel size Dx and Dy are not equal: %g != %g\n",dx,dy);

    /* create output image */

    outstr = stropen(getparam("out"), "w");
    create_image(&optr,nx,ny);                     
    Dx(optr) = Dx(iptr);
    Dy(optr) = Dy(iptr);
    Xmin(optr) = Xmin(iptr);
    Ymin(optr) = Ymin(iptr);

    /* create and set kernel image (only 1 quadrant needed) */

    create_image(&dptr,nx,ny); 
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++)
	if (i>0 || j>0) 
	  DIS(i,j) = 1.0/sqrt((double)(i*i + j*j));
	else
	  DIS(0,0) = 3.54;

    /* convolve input with kernel 
       this is a very expensive operation, so if requested, 
       keep the user informed, so he can abort if taking too long
    */

    while (nbench--) {
      for (j=0; j<ny; j++) {
	for (i=0; i<nx; i++) {
	  sum = 0.0;
	  count++;
	  if (report && count % report == 0) {
	    printf("%3d%% done\r", (int) (100*count/(nx*ny)));
	    fflush(stdout);
	  }
	  for (l=0; l<ny; l++) {
	    l1 = QABS(j,l);
	    for (k=0; k<nx; k++) {
	      k1 = QABS(k,i);
	      sum += CVI(k,l)*DIS(k1,l1); 
	    }
	  }
	  CVO(i,j) = -sum*gravc;     /* note that this is now in the correct units */
	}
      } 
      write_image(outstr, optr);
    }
}
