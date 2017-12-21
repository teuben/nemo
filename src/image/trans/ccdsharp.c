/* 
 * CCDSHARP: image 'enhancement' tricks
 *
 *	19-sep-96   q&d                             pjt
 *      5-apr-96    added divergence and vorticity  PJT
 *	11-apr-96   forgot to copy proper header elements   PJT
 *      17-apr-2017   laplace (classic) vs. lapabs
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
	"mode=laplace\n	Modes (laplace, lapabs, aregan, pregan, divergence, vorticity)",
	"VERSION=0.3\n  17-apr-2017 PJT",
	NULL,
};

string usage = "enhance/sharpen an image";



#define CV1(x,y,z)  CubeValue(iptr1,x,y,z)
#define CV2(x,y,z)  CubeValue(iptr2,x,y,z)
#define CVO(x,y,z)  CubeValue(optr,x,y,z)

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
    stream  instr, outstr;
    int     nx, ny, nz, mode;
    int     i,j,k;
    imageptr iptr1=NULL, iptr2=NULL, optr;      /* pointer to images */
    real    d1, d2, d3, d4, dx, dy;
    bool    Qsym = TRUE;            /* symmetric derivates w.r.t. pixel point */

    match(getparam("mode"),valid_modes,&mode);
    if (mode==0) error("Not a valid mode; valid:%s",valid_modes);
    dprintf(0,"Image sharpening method #%d\n",mode);

    instr = stropen(getparam("in"), "r");
    read_image( instr, &iptr1);
    nx = Nx(iptr1);	
    ny = Ny(iptr1);
    nz = Nz(iptr1);
    dx = Dx(iptr1);
    dy = Dy(iptr1);
    if (mode & MODE_DIV || mode & MODE_VORT) {
        if (read_image(instr,&iptr2) == 0)
            error("No second image found in %s\n",getparam("in"));
        if (nx != Nx(iptr2))  
            error("Second image doesn't match in NX: %d <> %d\n",Nx(iptr2),nx);
        if (ny != Ny(iptr2))  
            error("Second image doesn't match in NY: %d <> %d\n",Ny(iptr2),ny);
        if (nz != Nz(iptr2))  
            error("Second image doesn't match in NZ: %d <> %d\n",Nz(iptr2),nz);
    }

    outstr = stropen(getparam("out"), "w");
    create_cube(&optr,nx,ny,nz);
    Dx(optr) = Dx(iptr1);
    Dy(optr) = Dy(iptr1);
    Dz(optr) = Dz(iptr1);
    Xmin(optr) = Xmin(iptr1);
    Ymin(optr) = Ymin(iptr1);
    Zmin(optr) = Zmin(iptr1);
    /* should do the others too */

    if (mode & MODE_LAPLACE) {
        dprintf(1,"mode=laplace\n");
        for (k=0; k<nz; k++) {
            for (j=1; j<ny-1; j++) {
                for (i=1; i<nx-1; i++) {
		    d1 = CV1(i-1,j,k) - CV1(i,j,k);
		    d2 = CV1(i+1,j,k) - CV1(i,j,k);
		    d3 = CV1(i,j-1,k) - CV1(i,j,k);
		    d4 = CV1(i,j+1,k) - CV1(i,j,k);
                    CVO(i,j,k) = d1+d2+d3+d4;
                }
                CVO(0,j,k) = 0.0;
                CVO(nx-1,j,k) = 0.0;
            }
            for (i=0; i<nx; i++) {
                CVO(i,0,k) = 0.0;
                CVO(i,ny-1,k) = 0.0;
            }
        }
    } else if (mode & MODE_LAPABS) {
        dprintf(1,"mode=lapabs\n");
        for (k=0; k<nz; k++) {
            for (j=1; j<ny-1; j++) {
                for (i=1; i<nx-1; i++) {
                    d1 = CV1(i,j,k) - CV1(i-1,j,k);
                    d2 = CV1(i,j,k) - CV1(i+1,j,k);
                    d3 = CV1(i,j,k) - CV1(i,j-1,k);
                    d4 = CV1(i,j,k) - CV1(i,j+1,k);
                    CVO(i,j,k) = sqrt(d1*d1+d2*d2+d3*d3+d4*d4);
                }
                CVO(0,j,k) = 0.0;
                CVO(nx-1,j,k) = 0.0;
            }
            for (i=0; i<nx; i++) {
                CVO(i,0,k) = 0.0;
                CVO(i,ny-1,k) = 0.0;
            }
        }
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
    
    write_image(outstr, optr);

}
