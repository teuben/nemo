/*
 *  SDSFITS: convert HDF SDS into a FITS image/cube (2d/3d only)
 *
 *      19-jan-95	V1.0  toy model, for Jim Stone		Peter Teuben
 *      13-apr-95       V1.1  also made it work for 2D            pjt
 *	13-apr-96	V1.2  optionally dump out axis descriptor pjt 
 *      10-nov-02       V1.3  add HISTORY
 *      14-jan-03       V1.4  add flip= for Eve Ostriker          pjt
 *       2-dec-03       V1.5  add dummy=                          pjt
 *      11-dec-04       V1.6  add wcs=
 */

#include <stdinc.h>
#include <getparam.h>
#include <history.h>
#include <fitsio_nemo.h>

/* #include <hdf.h> */

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (SDS HDF file)",
    "out=\n			Optional output FITS file",
    "select=1\n			Select which SDS (1=first)",
    "axis=0\n                   Select this axis for 1D fits output",
    "reorder=f\n                Reorder the axes (true means memory intensive)",
    "dummy=t\n                  Write dummy axes also ?",
    "wcs=f\n                    Write a WCS of some sort",
    "VERSION=1.6\n		11-dec-03 PJT",
    NULL,
};

string usage="convert SDS HDF to FITS";

string cvsid="$Id$";

#define MAXRANK    3

void invert_array(int n, int *idx);
float *reorder_array(float *a, int size, int n, int *idx);

void nemo_main()
{
  string infile, outfile, *hitem;
  char label[128], unit[128], format[128], coordsys[128],  fitskey[80];
  float *buff, *fp, fmin, fmax, ratio, **coord;
  float crpix[MAXRANK], crval[MAXRANK], cdelt[MAXRANK];
  char ctype[3][80];
  int i, j, k, nsds, ret, size, rank, dimsizes[MAXRANK];
  int nx, ny, nz, num_type, select, axis, axlen, count = 0;
  FITS *ff;
  bool Qreorder = getbparam("reorder");
  bool Qdummy = getbparam("dummy");

  select = getiparam("select");
  axis = getiparam("axis");

  /*
   *  Open SDS and read and display it's header how many SDS# we have
   */

  infile = getparam("in");
  nsds = DFSDndatasets(infile);
  if (nsds<0) error("Error interpreting SDS file %s",infile);
  dprintf(0,"%s: has %d scientific data sets\n",infile,nsds);
  if (select < 1 || select > nsds)
    error("Illegal select=%d, #SDS=%d\n",select,nsds);

  for (count=1; count<=select; count++) {      /* loop over the SDS until we got the right one */
    
    ret = DFSDgetdims(infile, &rank, dimsizes, MAXRANK);
    if (ret)
      error("Error obtaining rank and dimensions (MAXRANK=%d)",MAXRANK);
    dprintf(0,"%s: found SDS rank %d: (",infile, rank);
    for (i=0, size=1; i<rank; i++) {
      dprintf(0,"%d%s",dimsizes[i], i==rank-1 ? ")" : ",");
      size *= dimsizes[i];
    } /* i */
    dprintf(0," - total of %d \"%s\"\n", 
	    size,
	    rank == 3 ? "voxels" : (rank == 2 ? "pixels" : "thingos"));
    
    for (i=rank-1, j=1; i>=0; i--, j++) {		/* reverse order !! */
      ret = DFSDgetdimstrs(i+1,label,unit,format);
      if (ret) {
	warning("Problem getdimstr for axis %s -- ignoring",j);
	continue;
      }
      dprintf(0,"Axis %d [%d] %s %s %s\n",
	      j,dimsizes[i],
	      label,unit,format);
    } /* i,j */
    if (count==select) { 
      if (axis > 0) {
	axlen = dimsizes[rank-axis];
	dprintf(0,"Also retrieving axis %d with length %d\n",axis,axlen);
      }
      coord = (float **) allocate(rank*sizeof(float *));
      for (i=0; i<rank; i++) {
	coord[i] = (float *) allocate(dimsizes[i]*sizeof(float));
	ret = DFSDgetdimscale(i+1,dimsizes[i],coord[i]);
	dprintf(0,"Axis %d [%d] WCS range: %g : %g \n",
	      rank-i,dimsizes[i],
	      coord[i][0], coord[i][dimsizes[i]-1]);
      }
    }
    ret = DFSDgetdatastrs(label,unit,format,coordsys);
    if (ret==0)
      dprintf(0,"Data: %s %s %s %s\n",label,unit,format,coordsys);
    fmax = fmin = 0.0;
  } /* count */
  
    /*
     * If some kind of output requested
     */
  
  if (hasvalue("out") && axis==0) {                 /* dump out the data for this SDS */
        /* 
         * read all the SDS data; this must happen around here ...
         */
        ret = DFSDgetNT(&num_type);
        if (ret) error("Getting data type");
        if (num_type != 5)
            error("Can only read REAL32 floating point data now");
        buff = (float *) allocate(sizeof(float)*size);
        ret = DFSDgetdata(infile, rank, dimsizes, buff);   /* get all data */
        if (ret) error("%d: Error reading data from SDS %s",ret,infile);
	for (i=0; i<rank; i++) {
	  j = rank-i-1;
	  crpix[j] = 1.0;
	  crval[j] = coord[i][0];
	  if (dimsizes[i] == 1) {
	    cdelt[j] = 1.0;
	  } else {
	    cdelt[j] = coord[i][1] - coord[i][0];
	    if (cdelt[j] == 0) {
	      warning("axis %d: resetting cdelt=0 -> 1",j+1);
	      cdelt[j] = 1;
	    }
	    if (dimsizes[i] > 2) {
	      ratio = (coord[i][2] - coord[i][1])/cdelt[j];
	      if (ratio != 1) warning("axis %d: cdelt_growth=%g\n",j+1,ratio);
	      dprintf(1,"axis %d: cdelt_growth=%g\n",j+1,ratio);
	    }
	  }
	}

        /*
         *  revert axes for output FITS file (HDF is native C, so by default
	 *  it reverts the axis; our keyword reorder=t restores the C order)
	 *
         */

    	outfile = getparam("out");
        dprintf(0,"%s: writing as FITS file\n",outfile);
	if (Qreorder) {                    /* Reordering can be a lot of work */
	  warning("if wcs=t and/or rank=3 be warned output may be flawed");
	  buff = reorder_array(buff, size, rank, dimsizes);
	  if (rank==2) {
     	    nz = 1;
	    ny = dimsizes[1];
	    nx = dimsizes[0];
	  } else if (rank==3) {
	    nx = dimsizes[2];
            ny = dimsizes[1];
            nx = dimsizes[0];
	  } else
            error("(%d) Cannot deal with rank != 2 or 3",rank);
	} else {                        /* this is actually the default */
	  invert_array(rank,dimsizes);
	  if (rank==2) {
     	    nz = 1;
	    ny = dimsizes[1];
	    nx = dimsizes[0];
	  } else if (rank==3) {
	    nz = dimsizes[2];
            ny = dimsizes[1];
            nx = dimsizes[0];
	    if (!Qdummy && nx==1) {
	      warning("trying to dummyfy NAXIS1");
	      rank=2;
	      nz = 1;
	      ny = dimsizes[2];
	      nx = dimsizes[1];
	      dimsizes[0] = nx;
	      dimsizes[1] = ny;
	    }
	  } else
            error("(%d) Cannot deal with rank != 2 or 3",rank);
	}
	ff = fitopen(outfile,"new",rank, dimsizes);

	strcpy(ctype[0],"X");   /* by lack of anything better....*/
	strcpy(ctype[1],"Y");
	strcpy(ctype[2],"Z");
  
	for (i=0; i<rank; i++) {
	  sprintf(fitskey,"CRPIX%d",i+1);
	  fitwrhdr(ff,fitskey,crpix[i]);
	  sprintf(fitskey,"CDELT%d",i+1);
	  fitwrhdr(ff,fitskey,cdelt[i]);
	  sprintf(fitskey,"CRVAL%d",i+1);
	  fitwrhdr(ff,fitskey,crval[i]);
	  sprintf(fitskey,"CTYPE%d",i+1);
	  fitwrhda(ff,fitskey,ctype[i]);
	}

        /*
         * and stuff it into fits format
         * but we need to revert the order of the axes
         * because of Fortran/C
         */

	fitwrhda(ff,"ORIGIN","NEMO: sdsfits w/history");

	/* only write 1 line */
	hitem = ask_history();
	fitwra(ff,"HISTORY",*hitem);

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
  } else if (hasvalue("out") && axis > 0) {     /* dump out this coordinate */

        if (Qreorder) error("Cannot deal with reorder=t here");
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

	fitwrhda(ff,"ORIGIN","NEMO: sdsfits w/ history");

	/* only write 1 line */
	hitem = ask_history();
	fitwra(ff,"HISTORY",*hitem);
	fitwra(ff,"COMMENT","Note the data are coordinate values, no WSC defined");

	if (fmin == fmax) {
	    fmin = fmax = coord[axis-1][0];
	    for (i=0; i<axlen; i++) {
	        fmin = MIN(fmin, coord[axis-1][i]);
	        fmax = MAX(fmax, coord[axis-1][i]);
	    }
	    dprintf(0,"Data min/max = %g %g\n",fmin,fmax);
	}
	if (fmin==fmax) warning("DATAMIN and DATAMAX are the same");
        fitwrhdr(ff,"DATAMIN", fmin);
        fitwrhdr(ff,"DATAMAX", fmax);
	fitwrite(ff, 0, coord[axis-1]);         /* write one line */
        fitclose(ff);
    } 

}


/* invert all elements of an array idx[n] */

void invert_array(int n, int *idx)
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

/* reorder the array/matrix: brute force approach (only 2D and 3D) */

float *reorder_array(float *a, int size, int n, int *idx)
{
  int x0, y0, z0, nx0, ny0, nz0;
  float *b = (float *) allocate(size*sizeof(float));

  if (n==2) {
    nx0 = idx[0];
    ny0 = idx[1];
    dprintf(0,"reorder: %d x %d\n",nx0,ny0);
    for (x0=0; x0<nx0; x0++)
      for (y0=0; y0<ny0; y0++)
	b[x0+nx0*y0] = a[y0+ny0*x0];
  } else {
    nx0 = idx[0];
    ny0 = idx[1];
    nz0 = idx[2];
    dprintf(0,"reorder: %d x %d x %d\n",nx0,ny0,nz0);
    for (x0=0; x0<nx0; x0++)
      for (y0=0; y0<ny0; y0++)
	for (z0=0; z0<nz0; z0++)
	b[x0+nx0*y0+nx0*ny0*z0] = a[z0+nz0*y0+nz0*ny0*x0];
  }		
  free(a);
  return b;
}

