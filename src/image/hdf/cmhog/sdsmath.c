/*
 *  SDSMATH:	fiddle with a cmhog produced set of 3 SDS files,
 *		Created for Witold Maciejewski - but unfinished
 *
 *       7-apr-98       V1.0  created, clones of hdfgrid actually
 *                            and writing out using the 'rules' in hdfall.src
 */

 
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <history.h>
#include <image.h>

#ifdef INC_HDF
#include <hdf.h> 	/* some conflicts with nemo include files */
#endif

string defv[] = {
    "in=???\n			Input cmhog file (HDF SD)",
    "out=???\n                  Output cmhog file (HDF SD)",
    "ome1=\n                    Some kind of ......",
    "VERSION=1.0\n		7-apr-98 PJT",
    NULL,
};

string usage="**unfinished**: Manipulate a CMHOG polar HDF SDS (den,vr,vt) dataset";

#define MAXRANK 2

local int rank, shape[MAXRANK], run[MAXRANK];
local char label[256], unit[256], fmt[256], coordsys[256];


extern string *burststring(string, string);

void nemo_main()
{
    float **image, **coord, *buffer, *rads, *phis, rad, phi, phi_orig, 
         dmin, dmax, phi1, phi2, rad1, rad2;
    float *buffer1, *buffer2, *buffer3;
    float **image1, **image2, **image3;
    int i, j, k, ret, size, type, nsds, isel, nx, ny;
    int ir, ip, jr, jp, nr, np, n1, n2;
    char **output, *cbuff;
    char ntype[32];
    string filter, zvar, infile = getparam("in");
    string outfile = getparam("out");
    real xrange[3], yrange[3], cosp, sinp;
    real x, y, a1, a2, c1,c2,c3,c4, cmin, cmax, dcon, tmp, vr, vt;
    imageptr iptr;
    bool mirror, first = TRUE, both = FALSE, flip=FALSE;

    nsds = DFSDndatasets(infile);
    if (nsds<0) 
        error("%s is probably not an HDF scientific dataset",infile);
    dprintf(1,"Found %d scientific data set%s in %s\n",
            nsds,(nsds > 1 ? "s" : ""),infile);

    
    for (k=0; k<nsds; k++) {        /* read until we have the right SDS */

  	ret = DFSDgetdims(infile,&rank, shape, MAXRANK);
    	if (ret < 0) error("Problem getting rank/shape at SDS #%d",k+1);

    	label[0] = unit[0] = fmt[0] = coordsys[0] = 0;
        ret = DFSDgetdatastrs(label, unit, fmt, coordsys);
        ret = DFSDgetNT(&type);

    	if (k==0) {				/* first time: allocate */
            coord = (float **) allocate(rank*sizeof(float *));
            for (i=0, size=1; i<rank; i++) {
    	        size *= shape[i];
                coord[i] = (float *) allocate(shape[i] * sizeof(float));
                ret = DFSDgetdimscale(i+1, shape[i],coord[i]);
                if (ret<0) error("getting shape[%d]",i+1);
                dprintf(0,"Dimension %d  Size %d\n",i+1,shape[i]);
            }
            buffer1 = (float *) allocate(size * sizeof(float));
            buffer2 = (float *) allocate(size * sizeof(float));
            buffer3 = (float *) allocate(size * sizeof(float));
        }
	if (k==0) buffer = buffer1;
	if (k==1) buffer = buffer2;
	if (k==2) buffer = buffer3;

        ret = DFSDgetdata(infile,rank,shape,buffer);
        dmin = dmax = buffer[0];
        for (i=1; i<size; i++) {
            dmin = MIN(dmin,buffer[i]);
            dmax = MAX(dmax,buffer[i]);
        }
        dprintf(0,"%d Datamin/max read: %g %g\n",k+1,dmin,dmax);
    }

    /* use some convenient set of pointers to handle the data in 2D format */

    image1 = (float **) allocate(shape[0] *sizeof(float *));        /* VR */
    for (i=0; i<shape[0]; i++) {
    	image1[i] = &buffer1[i*shape[1]];
    }
    image2 = (float **) allocate(shape[0] *sizeof(float *));        /* VT */
    for (i=0; i<shape[0]; i++) {
        image2[i] = &buffer2[i*shape[1]];
    }
    image3 = (float **) allocate(shape[0] *sizeof(float *));        /* DEN */
    for (i=0; i<shape[0]; i++) {
        image3[i] = &buffer3[i*shape[1]];
    }

        
    nr = shape[1];
    np = shape[0];
    rads = coord[1];
    phis = coord[0];


    /* the image can now be referred to as in:  image[phi][rad]  */
    /* with phi=0..np-1 and rad=0..nr-1                          */

    dprintf(0,"NR=%d NP=%d\n",nr,np);

#if 0    
c       external : dssdims,dssdast,dssdisc,dsadata
c
    ret = dssdims(rank,shape)                               DFSDsetdims    
    ret = dssdisc(1,shape(1),yscale)                        DFSDsetdimscale
    ret = dssdisc(2,shape(2),zscale)                        
c
    write(string,"('R-VELOCITY AT TIME=',1pe8.2)") time
    ret = dssdast(string,'km/sec   ','1pe8.2','Polar')      DFSDsetdatastrs
    ret = dsadata(filename,rank,shape,data)                 DFSDadddata    
c
    write(string,"('PHI-VELOCITY AT TIME=',1pe8.2)") time
    ret = dssdast(string,'km/sec   ','1pe8.2','Polar')
    ret = dsadata(filename,rank,shape,data)
c
    write(string,"('DENSITY AT TIME=',1pe8.2)") time
    ret = dssdast(string,'Msolar/pc**2','1pe8.2','Polar')
    ret = dsadata(filename,rank,shape,data)
    


    ret = DFSDadddata(outfile, rank, dimsizes, data);                
    ret = DFSDsetdatastr(label, unit, format, coordsys);
    ret = DFSDsetdims(rank, dimsizes)
    ret = DFSDsetdimscale(dim, dimsize, scale)
    
#endif

    ret = DFSDsetdims(rank, shape);
/*    ret = DFSDsetdimscale(1,			*/

    ret = DFSDgetdatastrs(label, unit, fmt, coordsys);
    
}
