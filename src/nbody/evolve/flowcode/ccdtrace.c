/* 
 * CCDTRACE:  trace a set of coordinates in an image and return
 *            the interpolated values from the image. Doesn't work
 *            for 3D cubes (yet).
 *
 *
 *  31-jul-96  Created (for use by Mike Regan to aid printing
 *             out additional values from "orbits" produced by flowcode)
 *  26-sep-02  file_lines fix
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
    "in=???\n       Input Image(5NEMO) file",
    "out=???\n      Output table of X,Y,I coordinates",
    "xytab=???\n    Input table of X,Y coordinates",
    "nmax=10000\n   Allocation space for input XY table, if needed",
    "VERSION=1.0b\n  10-oct-02 PJT",
    NULL,
};

string usage = "interpolate a set of coordinates in an image";

nemo_main()
{
    get_image();
    get_table();
    trace_image();
}


local imageptr iptr = NULL;
local double   idx, idy, xmin, ymin;
local int      nx, ny;

get_image()
{
    string input = getparam("in");
    stream instr = stropen(input,"r");

    read_image (instr,&iptr);

    idx = 1.0/Dx(iptr);          /* set some constants to remember */
    idy = 1.0/Dy(iptr);
    xmin = Xmin(iptr);
    ymin = Ymin(iptr);
    nx = Nx(iptr);
    ny = Ny(iptr);

    dprintf(0,"Read %d x %d image %s\n",nx,ny,input);
    dprintf(0,"Image Xrange: %g %g\n",xmin-0.5/idx,xmin+(nx-0.5)/idx);
    dprintf(0,"Image Yrange: %g %g\n",ymin-0.5/idy,ymin+(ny-0.5)/idy);

    if (Nz(iptr) > 1) error("Cannot handle 3D cubes yet");
#if 0
    if (idx != idy) {
        if (idx == -idy && xmin == -ymin) {     /* try and patch it */
            idx = -idx;
            xmin = -xmin;
            warning("Astronomical coordinate convention assumed");
        } else
            warning("Possible bug when using dx != dy");
    }
#endif
    if (idx<0) warning("1/Dx=%f",idx);
    if (idy<0) warning("1/Dy=%f",idy);
    dprintf(1,"Offset and scale factors: xmin,ymin,1/dx,1/dy=%f %f %f %f\n",
            xmin,ymin,idx,idy);
}


local real  *x, *y; 			/* XY data from table */
local int    npt;			/* number of XY data points */

get_table()
{
    string input;
    stream instr;
    real *coldat[2];
    int i, j, colnr[2];
    real xmint, xmaxt, ymint, ymaxt;
    int xcol=1, ycol=2;
    int nmax=getiparam("nmax");

    input = getparam("xytab");
    nmax = nemo_file_lines(input,nmax);
    dprintf(0,"Allocated %d lines for table %s\n",nmax,input);
    if (nmax < 1) error("Empty table");

    colnr[0] = xcol;
    colnr[1] = ycol;    
    coldat[0] = x = (real *) allocate(sizeof(real) * nmax);
    coldat[1] = y = (real *) allocate(sizeof(real) * nmax);

    instr = stropen(input,"r");
    npt = get_atable(instr,2,colnr,coldat,nmax);    /* get data */
    if (npt < 0) {
    	npt = -npt;
    	warning("Could only read %d data",npt);
    }

    /* go through the table, find and report min & max in X and Y */

    xmint = ymint =  HUGE;
    xmaxt = ymaxt = -HUGE;
    for (i=0; i<npt; i++) {
        xmaxt=MAX(x[i],xmaxt);
        xmint=MIN(x[i],xmint);
        ymaxt=MAX(y[i],ymaxt);
	ymint=MIN(y[i],ymint);
    }
    dprintf(0,"Table Xrange: %g %g\n",xmint,xmaxt);
    dprintf(0,"Table Yrange: %g %g\n",ymint,ymaxt);
}


#define F(ix,iy)	MapValue(iptr,ix,iy)

trace_image()
{
    double xp, yp, dx, dy, outval;
    int i, ix, iy, nout = 0;
    stream outstr;

    outstr = stropen(getparam("out"),"w");

    for (i=0; i<npt;i++) {
        xp = (x[i]-xmin)*idx + 0.5;        /* fractional cell index  0..nx */
        yp = (y[i]-ymin)*idy + 0.5;
        ix = (int) floor(xp);              /* cell index:   0 .. nx-1 */
        iy = (int) floor(yp);
        dx = xp-ix;			   /* index in a cell: 0.0 .. 1.0 */
        dy = yp-iy;
        dprintf(1,"%g %g %d %d %g %g\n",
            x[i],y[i],ix,iy,xp,yp);
        if (ix<1 || ix>nx-2 || iy<1 || iy>ny-2) {           /* outside grid */
            nout++;
            fprintf(outstr,"# %g %g Outside\n",x[i],y[i]);
        } else {                                            /* inside grid */
            outval = (1-dx)*(1-dy)*F(ix,iy)     +      /* See e.g.         */
	                 dx*(1-dy)*F(ix+1,iy)   +      /* Abram.& Steguhn  */
	                 (1-dx)*dy*F(ix,iy+1)   +      /* par. 25.2.66     */
			     dx*dy*F(ix+1,iy+1);       /* Four Foint Formu */
            fprintf(outstr,"%g %g %g\n",x[i],y[i],outval);
        }
    }
    dprintf(0,"Trace: %d/%d outside grid\n",nout,npt);
}



