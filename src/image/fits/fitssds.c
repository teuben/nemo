/*
 *  FITSSDS: convert FITS image/cube into an HDF SDS 
 *
 *      20-jan-95	V1.0  toy model, for Jim Stone		Peter Teuben
 */

#include <stdinc.h>
#include <getparam.h>
#include <fitsio.h>

/* #include <hdf.h> */

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input FITS file",
    "out=???\n			Output file (SDS HDF file)",
    "VERSION=1.0\n		20-jan-95 PJT",
    NULL,
};

string usage="convert FITS to SDS HDF";

#define MAXRANK    7

void nemo_main()
{
    string infile, outfile;
    char label[128], unit[128], format[128], coordsys[128];
    float *buff, *fp, fmin, fmax;
    int i, j, k, nsds, ret, size, rank, dimsizes[MAXRANK];
    int nx, ny, nz, num_type;
    FITS *ff;

    infile = getparam("in");
    ff = fitopen(infile,"old", MAXRANK, dimsizes);
    for (i=0, rank=0; i<MAXRANK; i++)
        if (dimsizes[i] > 1) rank++;
    dprintf(0,"%s: found rank %d: (",infile, rank);
    for (i=0, size=1; i<rank; i++) {
    	dprintf(0,"%d%s",dimsizes[i], i==rank-1 ? ")" : ",");
    	size *= dimsizes[i];
    }
    dprintf(0," - total of %d \"pixels\"\n",size);

    buff = (float *) allocate(size * sizeof(float));

    nx = dimsizes[0];
    ny = dimsizes[1];
    nz = dimsizes[2];

    fp = buff;
    for (k=0; k<nz; k++) {
        fitsetpl(ff, 1, &k);
        for (j=0; j<ny; j++) {
            fitread(ff, j, fp);
            fp += nx;
        }
    }
    fitclose(ff);

    fmin = fmax = buff[0];
    for (i=1; i<size; i++) {
        fmin = MIN(fmin, buff[i]);
        fmax = MAX(fmax, buff[i]);
    }
    dprintf(0,"Datamin/max: %g %g\n",fmin,fmax);

    outfile = getparam("out");
    dimsizes[0] = nz;
    dimsizes[1] = ny;
    dimsizes[2] = nx;
    ret = DFSDsetdims(rank, dimsizes);
    if (ret) error("Cannot setdims for new SDS file %s",outfile);
    ret = DFSDsetNT(5);
    if (ret) error("Cannot setNT for new SDS file %s",outfile);
    ret = DFSDputdata(outfile, rank, dimsizes, buff);
    if (ret) error("Writing SDS data to %s",outfile);
}






