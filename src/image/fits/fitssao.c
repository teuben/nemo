/*
 *   FITSSAO:    dump a fits file onto an SAO-IIS server 
 *		(e.g. ximtool, -- tested
 *		      saoimage -- untested
 *
 *	20-mar-95   cloned off fitsgids             pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <image.h>
#include <fitsio.h>
#include "sao-iis.h"

string defv[] = {
    "in=???\n		Input fits file",
    "datamin=\n         Datamin,if to override",
    "datamax=\n         Datamax,if to override",
    "x=\n               Range in X to display [all]",
    "y=\n               Range in Y to display [all]",
    "z=\n               Range in Y to display [first]",
    "frame=1\n          Frame to store image in",
    "VERSION=1.0\n	20-mar-95 PJT",
    NULL,
};

string usage = "display fits images on SAO-IIS display server";

#ifndef MAXSIZE
#define MAXSIZE  2048	/* max size in each dimension of an image */
#endif

local int make_idx(string, int *, int, int);

void nemo_main()
{
    FITS *fitsfile;
    string fitsname;
    int ndim=3, naxis[3], nx, ny, nz, i, j, k, p,
	idx[MAXSIZE], idy[MAXSIZE], idz[MAXSIZE];
    int nbval=0, plane, bitpix, frame;
    real bval_out, rmin, rmax, tmp;
    FLOAT *bufin, buf;
    Rawpix_t *bufout, *bp;
    float datamin, datamax, datarange, bscale, bzero;
    bool Qshort;

    fitsname = getparam("in");
    fitsfile = fitopen(fitsname,"old",ndim,naxis);
    if (fitsfile==NULL) error("Could not open FITS file %s",getparam("in"));
    dprintf(0,"%s: %d * %d * %d\n",
    		getparam("in"),naxis[0],naxis[1],naxis[2]);


    nx = make_idx(getparam("x"),idx,MAXSIZE,naxis[0]);
    ny = make_idx(getparam("y"),idy,MAXSIZE,naxis[1]);
    nz = make_idx(getparam("z"),idz,MAXSIZE,naxis[2]);

    if (nx < 1 || ny < 1 || nz < 1 )
	error("Cannot deal with (%d,%d,%d)",nx,ny,nz);

    if (hasvalue("datamin"))
        datamin=getdparam("datamin");
    else
        if (fitexhd(fitsfile,"DATAMIN"))
            fitrdhdr(fitsfile,"DATAMIN",&datamin,0.0); 
        else
            datamin=getdparam("datamin");
    if (hasvalue("datamax"))
        datamax=getdparam("datamax");
    else
        if (fitexhd(fitsfile,"DATAMAX"))
            fitrdhdr(fitsfile,"DATAMAX",&datamax,0.0); 
        else
            datamax=getdparam("datamax");

    dprintf(0,"Displaying (%d,%d,%d) with datamin=%g datamax=%g\n",
            nx, ny, nz, datamin, datamax);
    if (datamin >= datamax) error("Cannot display with datamin/max=%g %g",
    	    datamin, datamax);
    bscale = (datamax-datamin)/256.0;
    bzero = datamin;

    frame = getiparam("frame");

    /* buffers to read and write */
    bufin  = (FLOAT *) allocate(naxis[0]*sizeof(FLOAT));
    bufout = (Rawpix_t *) allocate(nx*ny*sizeof(Rawpix_t));
    rmin = HUGE;
    rmax = -HUGE;
    if (si_open() == FAIL) 
        error("si_open() failed");
    if (si_selectframe(frame) == FAIL)
        error("si_selectframe(%d) failed",frame);

    bp = bufout;
    for (k=0; k<nz; k++) {     
        p = idz[k]-1;
        fitsetpl(fitsfile,1,&p);
        for (j=0; j<ny; j++) {      /* loop over all rows */
            fitread(fitsfile,idy[j]-1,bufin);     /* read it from fits file */
            for (i=0; i<nx; i++) {                /* stuff it in memory */
                buf = bufin[idx[i]-1];
                rmin=MIN(rmin,buf);
                rmax=MAX(rmax,buf);
                /* clip datamin:datamax */
                if (buf < datamin)
                    *bp++ = 0;
                else if (buf > datamax)
                    *bp++ = 255;
                else
                    *bp++ = (Rawpix_t) ((buf-datamin)*255/(datamax-datamin));
            }

        }
        if (si_send_image(nx,ny,nx/2,ny/2,bufout) == FAIL)
            error("si_send_image() failed");
    }
    fitclose(fitsfile);

    dprintf(0,"Clipped %d pixels; Found datamin=%g datamax=%g Last REC=%d\n",
            nbval,rmin,rmax,plane-1);
    si_close();
}

local int make_idx(string axis, int *idx, int maxlen, int axlen)
{
    int i, j, n;
    if (axlen>maxlen)
        warning("Axis %s is too long (%d) to hold %d pixels",axis,axlen,maxlen);
    n = nemoinpi(axis,idx,maxlen);
    if (n>0) {
    	for (i=0, j=0; i<n; i++) {
            if (idx[i]<1 || idx[i]>axlen) continue; /* don't copy: */
            idx[j++] = idx[i];
    	}
    	if (j!=n)
            warning("Only %d/%d elements out of range %s used",j,n,axis);
        return j;
    }
    if (n<0) error("Parsing error %s",axis);
    for (i=0; i<axlen; i++) idx[i] = i+1;
    return axlen;
}
