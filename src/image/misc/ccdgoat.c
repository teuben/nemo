/* 
 * CCDGOAT:   LV diagram diagnostics
 *
 *	31-oct-90	created		PJT
 *	30-mar-97	usage
 *	 7-may-01    	cleaned up superfluous #define's	PJT
 *
 *                      
 */


#include <stdinc.h>             /* also gets <stdio.h>  */
#include <getparam.h>
#include <vectmath.h>           /* otherwise NDIM undefined */
#include <filestruct.h>
#include <image.h>

char *defv[] = {
        "in=???\n       Input image file - must be an LV diagram",
        "vel=0.0\n      Velocity cutoff for Ndiff table",
        "temp=1.0\n     Temperature to measure Isotherms at",
        "mode=goat\n    Mode {goat, int, vdiff}",
        "VERSION=1.0c\n 7-may-01 PJT",
        NULL,
};

string usage="LV diagram diagnostics (PJT)";



nemo_main()
{
    stream  instr;
    imageptr iptr = NULL;
    real vel, temp;
    string mode;
    bool scanopt();

    instr = stropen(getparam("in"), "r");     /* get file name and open file */
    read_image( instr, &iptr);                /* read image */
    strclose(instr);                          /* close file */
    check_image(instr,iptr);
    mode = getparam("mode");
    if (scanopt(mode,"goat")) {
        vel = getdparam("vel");
        goat(iptr, vel);                         /* make table */
    } else if (scanopt(mode,"int")) {
        integral(iptr);
    } else if (scanopt(mode,"vdiff")) {
        temp = getdparam("temp");
        vdiff(iptr,temp);
    } else
        error("Wrong mode: choose any of {goat, int, vdiff}");
}
check_image(iptr)
imageptr iptr;
{
}

/*
 * GOAT:   in 1st and 4th quadrant only take emission below -vel and
 *	   above +vel, in 2nd and 3rd take it all, since its outside
 *	   solar circle
 */
goat(iptr, vel)
imageptr iptr;
real vel;
{
    int  ix, iy, nx, ny, ivel, itest;
    real phi, diff, sumN, sumS;
    
    nx = Nx(iptr);
    ny = Ny(iptr);
    if (nx%4)
        warning("Longitude is not multiple of 4 pixels");
    itest = -1;             /* no testing */

    ivel = (int)  ( (vel-Ymin(iptr))/Dy(iptr) );        /* Vel cutoff pixel */
    if (ivel<0 || ivel>=ny) {
        warning("vel=%f not inside range of Y-axis(ivel=%d)",vel,ivel);
        return(0);
    } else
	dprintf(0,"ivel = %d\n",ivel);
    phi = 0.5 * nx * Dx(iptr) + Xmin(iptr);             /* first long. point */
    dprintf(0,"phi (sumS-sumN)/sumS\n");
    for (ix=nx/2; ix<3*nx/4; ix++) {            /* 1st and 4th quadrant */
        sumN = sumS = 0.0;
        for(iy=ivel; iy<ny; iy++) {             
            sumN += MapValue(iptr,ix,ny-1-iy);  /* v<-vel in 1st quadrant */
            sumS += MapValue(iptr,nx-1-ix,iy);  /* v>vel in 4th quadrant */
            if (ix==itest)  dprintf(1,"iy : %d %f %f %f %f\n",
                iy, MapValue(iptr,ix,iy),MapValue(iptr,nx-1-ix,ny-1-iy),
                    sumN, sumS);
        }
        if (sumS == 0.0)
            diff = 0.0;
        else
            diff = (sumS-sumN)/sumS;
        printf("%g %g %g %g\n",phi,diff,sumS,sumN);
        phi += Dx(iptr);
    }
    for (ix=3*nx/4; ix<nx; ix++) {            /* 2nd and 3rd quadrant */
        sumN = sumS = 0.0;
        for(iy=0; iy<ny; iy++) {
            sumN += MapValue(iptr,ix,iy);
            sumS += MapValue(iptr,nx-1-ix,ny-1-iy);
            if (ix==itest)  dprintf(1,"iy : %d %f %f %f %f\n",
                iy, MapValue(iptr,ix,iy),MapValue(iptr,nx-1-ix,ny-1-iy),
                    sumN, sumS);
        }
        if (sumS == 0.0)
            diff = 0.0;
        else
            diff = (sumS-sumN)/sumS;
        printf("%g %g %g %g\n",phi,diff,sumS,sumN);
        phi += Dx(iptr);
    }
}

/*
 * Integral:  compute total column densities as function of longitude
 *            It only does 0..pi, and gives both the North and South
 *            column densities for direct comparison of asymmetries
 */

integral(iptr)
imageptr iptr;
{
    int  ix, iy, nx, ny, ivel, ix1,iy1;
    real phi, diff, sumN, sumS;
    
    nx = Nx(iptr);      
    ny = Ny(iptr);

    phi = 0.5 * nx * Dx(iptr) + Xmin(iptr);             /* first long. point */
    dprintf(0,"phi sumS sumN\n");
    for (ix=nx/2; ix<nx; ix++) {                    /* loop each longitude */
        sumN = sumS = 0.0;
        ix1=nx-1-ix;
        for(iy=0; iy<ny; iy++) {
            iy1 = ny-1-iy;
            sumN += MapValue(iptr,ix,iy);
            sumS += MapValue(iptr,ix1,iy1);
            if (ix==nx/2)  dprintf(1,"iy : %d %f %f %f %f\n",
                iy, MapValue(iptr,ix,iy),MapValue(iptr,ix1,iy1),
                    sumN, sumS);
        }
        printf("%g %g %g\n",phi,sumS,sumN);
        phi += Dx(iptr);
    }
}

/*
 *   VDIFF: compute the vdiff's of N and S
 */

vdiff(iptr,t)
imageptr iptr;
real t;
{
    int  ix, iy, nx, ny, ivel, ix1,iy1;
    real phi, vdiff, vN, vS, oldt, newt;
    
    nx = Nx(iptr);      
    ny = Ny(iptr);

    phi = 0.5 * nx * Dx(iptr) + Xmin(iptr);             /* first long. point */
    dprintf(0,"phi sumS sumN\n");
    for (ix=nx/2; ix<nx; ix++) {                    /* loop each longitude */
        ix1 = ix;
        oldt = MapValue(iptr,ix,0);
        for (iy=1; iy<ny; iy++) {
            newt = MapValue(iptr,ix,iy);
            if (oldt < t && t < newt) {
                vN = Ymin(iptr) + Dy(iptr) *
                     ( iy-1.0 + Dy(iptr)/(newt-oldt) * (t-oldt));
                break;
            }
            oldt = newt;
        }
        if (iy>=ny) {
            warning("No N isotherm found at ix=%d",ix);
            continue;
        }

        ix1 = nx-1-ix;
        oldt = MapValue(iptr,ix1,ny-1);
        for (iy=ny-2; iy>=0; iy--) {
            newt = MapValue(iptr,ix1,iy);
            if (oldt < t && t < newt) {
                vS = Ymin(iptr) + Dy(iptr) *
                     ( iy + Dy(iptr)/(newt-oldt) * (t-newt));
                break;
            }
            oldt = newt;
        }
        if (iy<0) {
            warning("No S isotherm found at ix=%d",ix);
            continue;
        }
        printf("%f %f %f %f\n",phi,vN,vS,vS+vN);
        phi += Dx(iptr);
    }
}

