/*
 *  SDSFITS: convert HDF SDS into a FITS image/cube (2d/3d only)
 *
 *      19-jan-95	V1.0  toy model, for Jim Stone		Peter Teuben
 *      13-apr-95       V1.1  also made it work for 2D          pjt
 *	13-apr-96	V1.2  optionally dump out axis descriptor pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <fitsio.h>

/* #include <hdf.h> */

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (SDS HDF file)",
    "out=\n			Optional output FITS file",
    "select=1\n			Select which SDS (1=first)",
    "axis=0\n                   Select this axis for 1D fits output",
    "VERSION=1.2\n		13-apr-96 PJT",
    NULL,
};

string usage="convert SDS HDF to FITS";

#define MAXRANK    3

void nemo_main()
{
    string infile, outfile;
    char label[128], unit[128], format[128], coordsys[128];
    float *buff, *fp, fmin, fmax, *coord = 0;
    int i, j, k, nsds, ret, size, rank, dimsizes[MAXRANK];
    int nx, ny, nz, num_type, select, axis, axlen, count = 0;
    FITS *ff;

    select = getiparam("select");
    axis = getiparam("axis");

    /*
     *  Open SDS and read and display it's header
     */

    infile = getparam("in");
    nsds = DFSDndatasets(infile);
    if (nsds<0) error("Error interpreting SDS file %s",infile);
    dprintf(0,"%s: has %d scientific data sets\n",infile,nsds);


  for (count=1; count<=select; count++) {
    
    ret = DFSDgetdims(infile, &rank, dimsizes, MAXRANK);
    if (ret)
        error("Error obtaining rank and dimensions (MAXRANK=%d)",MAXRANK);
    dprintf(0,"%s: found SDS rank %d: (",infile, rank);
    for (i=rank-1, size=1; i>=0; i--) {		/* reverse order !! */
    	dprintf(0,"%d%s",dimsizes[i], i==0 ? ")" : ",");
    	size *= dimsizes[i];
    }
    dprintf(0," - total of %d \"pixels\"\n",size);
    
    for (i=rank-1, j=1; i>=0; i--, j++) {		/* reverse order !! */
    	ret = DFSDgetdimstrs(i+1,label,unit,format);
    	if (ret) {
            warning("Problem getdimstr for axis %s -- ignoring",j);
            continue;
    	}
    	dprintf(0,"Axis %d: %s %s %s\n",j,label,unit,format);
    	if (count==select && axis==j) {          /* get coordinate info */
    	    axlen = dimsizes[i];
    	    dprintf(0,"Also retrieving axis %d with length %d\n",axis,axlen);
            coord = (float *) allocate(axlen*sizeof(float));
            ret = DFSDgetdimscale(i+1,axlen,coord);    	    
    	}
    }
    ret = DFSDgetdatastrs(label,unit,format,coordsys);
    if (ret==0)
        dprintf(0,"Data: %s %s %s %s\n",label,unit,format,coordsys);
    fmax = fmin = 0.0;
  }

    /*
     * If output requested
     */

    if (hasvalue("out") && axis==0) {
        /* 
         * read all the SDS data; this must happen around here ...
         */
        ret = DFSDgetNT(&num_type);
        if (ret) error("Getting data type");
        if (num_type != 5)
            error("Can only read REAL32 floating point data now");
        buff = (float *) allocate(sizeof(float)*size);
        ret = DFSDgetdata(infile, rank, dimsizes, buff);
        if (ret) error("%d: Error reading data from SDS %s",ret,infile);

        /*
         *  revert axes for output FITS file
         */

    	outfile = getparam("out");
        dprintf(0,"%s: writing as FITS file\n",outfile);
        invert(rank,dimsizes);
        if (rank==2) {
     	    nz = 1;
	    ny = dimsizes[1];
	    nx = dimsizes[0];
        } else if (rank==3) {
	    nz = dimsizes[2];
            ny = dimsizes[1];
            nx = dimsizes[0];
        } else
            error("(%d) Cannot deal with rank != 2 or 3",rank);
	ff = fitopen(outfile,"new",rank, dimsizes);

        /*
         * and stuff it into fits format
         * but we need to revert the order of the axes
         * because of Fortran/C
         */

	fitwrhda(ff,"ORIGIN","NEMO: sdsfits");
	if (fmin == fmax) {
	    fmin = fmax = buff[0];
	    for (i=0, fp=buff; i<size; i++, fp++) {
	        fmin = MIN(fmin, *fp);
	        fmax = MAX(fmax, *fp);
	    }
	    dprintf(0,"Data min/max = %g %g\n",fmin,fmax);
	}
	if (fmin==fmax) warning("DATAMIN and DATAMAX are the same");
        fitwrhdr(ff,"DATAMIN", fmin);
        fitwrhdr(ff,"DATAMAX", fmax);
	fp = buff;
	for (k=0; k<nz; k++) {
            fitsetpl(ff, 1, &k);
            for (j=0; j<ny; j++) {
                fitwrite(ff, j, fp);
                fp += nx;
            }
	}
        fitclose(ff);
    } else if (hasvalue("out") && axis > 0) {

    	outfile = getparam("out");
        dprintf(0,"%s: writing axis %d as FITS file\n",outfile,axis);
        rank = 2;
	dimsizes[0] = axlen;
        dimsizes[1] = 1;
	ff = fitopen(outfile,"new",rank, dimsizes);

        /*
         * and stuff it into fits format
         * but we need to revert the order of the axes
         * because of Fortran/C
         */

	fitwrhda(ff,"ORIGIN","NEMO: sdsfits");
	if (fmin == fmax) {
	    fmin = fmax = coord[0];
	    for (i=0; i<axlen; i++) {
	        fmin = MIN(fmin, coord[i]);
	        fmax = MAX(fmax, coord[i]);
	    }
	    dprintf(0,"Data min/max = %g %g\n",fmin,fmax);
	}
	if (fmin==fmax) warning("DATAMIN and DATAMAX are the same");
        fitwrhdr(ff,"DATAMIN", fmin);
        fitwrhdr(ff,"DATAMAX", fmax);
	fitwrite(ff, 0, coord);         /* write one line */
        fitclose(ff);
    } 

}


/* invert all elements of an array idx[n] */

invert(int n, int *idx)
{
    int i, tmp, *a, *b;

    if (n<2) return;
    
    a = &idx[0];
    b = &idx[n-1];
    for (i=0; i<n/2; i++) {
        tmp = *a;
        *a = *b;
        *b = tmp;
        a++;
        b--;
    }
}

