/* 
 * CCDPRINT: print values at gridpoints, with optional numeric labels
 *
 *	 23-jan-89  V1.0 created        Peter Teuben
 *	 15-jul-89  V1.0a   - nemoinpX instead of nemoinp	PJT
 *	  8-jul-93  V1.1 blank default for x=,y=,z=             pjt
 *		    and  fixep various old style coding
 *	 14-sep-93  V1.2 Added index mode keyword offset=	pjt
 *	 14-oct-99  V1.2b fixed bug printing WCS labels		pjt
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
	"in=???\n  Input filename",
	"x=0\n              Pixels in X to print",
	"y=0\n              Pixels in Y to print",
	"z=0\n              Pixels in Z to print",
	"scale=1.0\n        Scale factor for printout",
	"format=%g\n        Format specification for output",
	"newline=f\n        Force newline between each number?",
	"label=\n	    Add x, y and or z labels add appropriate labels",
	"offset=0\n         Offset (0 or 1) to index coordinates X,Y,Z",
	"VERSION=1.2b\n     14-oct-99 PJT",
	NULL,
};

string usage = "print values at gridpoints of an image";


nemo_main()
{
    int     *ix, *iy, *iz, i, j, k, nx, ny, nz, nxpos, nypos, nzpos, offset;
    bool    newline, xlabel, ylabel, zlabel;
    real    f;
    string  infile;			        /* file name */
    stream  instr;				/* file stream */
    imageptr iptr=NULL;			      /* allocated dynamically */
    string   fmt, label;
    real     scale_factor, x, y, z;
    
    int      ini_array();
    bool     scanopt();

    scale_factor = getdparam("scale");
    fmt = getparam("format");
    instr = stropen (getparam("in"), "r");
    newline = getbparam("newline");
    label = getparam("label");
    xlabel = scanopt(label,"x");
    ylabel = scanopt(label,"y");
    zlabel = scanopt(label,"z");
    offset = getiparam("offset");
    if (read_image (instr,&iptr) == 0)
    	error("Problem reading image");
    
    nx = Nx(iptr);	                        /* cube dimensions */
    ny = Ny(iptr);
    nz = Nz(iptr);
    ix = (int *) allocate(nx*sizeof(int));        /* allocate position arrays */
    iy = (int *) allocate(ny*sizeof(int));
    iz = (int *) allocate(nz*sizeof(int));
    nxpos = ini_array("x",ix,nx,offset);
    nypos = ini_array("y",iy,ny,offset);
    nzpos = ini_array("z",iz,nz,offset);

    for (k=0; k<nzpos; k++) {				/* over all planes */
        z = Zmin(iptr) + iz[k] * Dz(iptr);
        if (!newline && zlabel) {
            printf("plane Z = ");
            myprintf(fmt,z);
            printf("\n\n");
        }
        for (j=nypos-1; j>=0; j--) {			/* over all rows */
            y = Ymin(iptr) + iy[j] * Dy(iptr);
            if (j==(nypos-1) && !newline && xlabel) {    /* pr. row of X coords */
                if (ylabel) printf(" Y\\X "); /* ? how to get correct length ? */
                for (i=0; i<nxpos; i++) {
                    x = Xmin(iptr) + ix[i] * Dx(iptr);
                    myprintf(fmt,x);
                }
                printf("\n\n");
            }
            if (!newline && ylabel) {   /* print first column of Y coord */
                myprintf(fmt,y);
                printf(" ");
            }
            for (i=0; i<nxpos; i++) {			/* over all columns */
                f = CubeValue(iptr,ix[i],iy[j],iz[k]);
                if (newline) {
                    x = Xmin(iptr) + ix[i] * Dx(iptr);
                    if (xlabel) myprintf(fmt,x);
                    if (ylabel) myprintf(fmt,y);
                    if (zlabel) myprintf(fmt,z);
                } 
                printf (fmt,f*scale_factor);
                if (newline)
                    printf("\n");
                else
                    printf (" ");
            }
            if (!newline) printf("\n");
        }
        if (!newline) printf("\n\n");
    }
    strclose(instr);
}

ini_array(key, dat, ndat, offset)
string key;             /* keyword */
int *dat;               /* array */
int ndat;               /* length of array */
int offset;             /* offset applied to index array */
{
    int i, n;

    if (ndat <= 0) return ndat;
    n = hasvalue(key) ? nemoinpi(getparam(key),dat,ndat) : 0;
    if (n > 0) {
        for (i=0; i<n; i++) {
            dat[i] -= offset;
            if (dat[i] < 0 || dat[i] >= ndat)
                error("Illegal boundaries in %s",key);
        }
    } else if (n==0) {
        for (i=0; i<ndat; i++)
            dat[i] = i;
        n = ndat;
    } else
        error("Error parsing %s=%s",key,getparam(key));
    return n;
}

myprintf(fmt,v)
string fmt;
real v;
{
    printf(fmt,v);
    printf(" ");
}
