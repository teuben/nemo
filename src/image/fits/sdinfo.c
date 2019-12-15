/*
 * sample program to use cfitsio on SDFITS files
 * benchmark some low level math
 *
 * 23-nov-2019   PJT       written
 * 11-dec-2019   PJT       mdarray reduction example
 *
 * Benchmark 6 N2347 files:  2.4"  (this is with mom=0 stats)
 * dims=5 for NGC5291:       31-35ms (depending in 1 or 3 levels)
 *
 * 
 */
#include <nemo.h>
#include <strlib.h>
#include <extstring.h>
#include <moment.h>
#include <mdarray.h>

#include <fitsio.h>  
#include <longnam.h>

string defv[] = {
    "in=???\n            Input SDFITS  fits file(s)",
    "stats=f\n           Doing simple stats",
    "mom=2\n             Highest moment for stats",
    "cols=\n             Column names to track - or names of the dimensions if given",
    "dims=\n             Dimensions to reduce [ex1: 2,11,2,2,4",
    "proc=\n             Reduction procedure (PS,FS,NOD,TP)",
    "bench=1\n           How many times to run benchmark",
    "VERSION=0.4\n       13-dec-2019 PJT",
    NULL,
};

string usage = "sdfits info and bench reduction procedure";


  // inline
float muladd(float x, float a, float b) {
    return x*a + b;
}


real dcmeantsys(int ndata, real *calon, real *caloff, real tcal)
{
  int nedge = ndata/10;
  int i;
  real tsys, sum1 = 0.0,  sum2 = 0.0;
  static int count = 0;

  for (i=nedge; i<ndata-nedge; i++) {
    sum1 += caloff[i];
    sum2 += (calon[i] - caloff[i]);
  }
  if (sum2 == 0.0)
    tsys = tcal;
  else
    tsys = tcal*sum1/sum2 + tcal/2.0;
  dprintf(1,"tsys=%g  %d %d %d\n",tsys,ndata,nedge,++count);
  return tsys;
}

void ta(int ndata, real *sig, real *ref, real tsys, real *ta)
{
  int i;

  for (i=0; i<ndata; i++)
    ta[i] = tsys * (sig[i]/ref[i] - 1.0);
}

void ta2(int ndata, real *sig1, real *sig2, real *ref1, real *ref2, real tsys, real *ta)
{
  int i;

  for (i=0; i<ndata; i++)
    ta[i] = tsys * ((sig1[i]+sig2[i])/(ref1[i]+ref2[i]) - 1.0);
  dprintf(1,"Ta(2) %g\n",ta[0]);
}

void average(int nvec, int ndata, real **data, real *ave)
{
  int i, j;

  for (i=0; i<ndata; ++i) ave[i] = 0;
  for (j=0; j<nvec; ++j)
    for (i=0; i<ndata; ++i)
      ave[i] += data[j][i];  // invalid read
  for (i=0; i<ndata; ++i) ave[i] /= nvec;  
}

void average2(int ndata, real *data1, real *data2, real *ave)
{
  int i, j;

  for (i=0; i<ndata; i++)
    ave[i] = 0.5*(data1[i]+data2[i]);
}

void nemo_main(void)
{
    fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    int status, fmode, ii, jj, i, j, k, colnum, data_col;
    long int nrows;
    int ncols, nchan, nfiles, found;
    string fname = getparam("in"), *fnames;
    string *colnames;
    int nsize, ndims, dims[MDMAXDIM];
    char template[64];
    char keyword[FLEN_KEYWORD], colname[FLEN_VALUE], data_fmt[FLEN_VALUE];
    bool Qstats = getbparam("stats");
    int mom = getiparam("mom");
    string *colcheck = burststring(getparam("cols"),", ");
    int ncolcheck = xstrlen(colcheck, sizeof(string))-1;
    int anynul = 0;
    int bench = getiparam("bench");
    real tsys;

    fnames = burststring(fname,", ");
    nfiles = xstrlen(fnames, sizeof(string)) -1;
    dprintf(0,"Found %d files\n",nfiles);
    dprintf(0,"Found %d cols to check\n",ncolcheck);

    if (hasvalue("dims")) {
      ndims = nemoinpi(getparam("dims"),dims,MDMAXDIM);
      for (i=0, nsize=1; i<ndims; i++)
	nsize *= dims[i];
      dprintf(0,"nsize=%d\n",nsize);
    } else 
      nsize = ndims = 0;

    for (j=0; j<nfiles; j++) {
      fname = fnames[j];
      status = 0;         /* initialize status before calling fitsio routines */
      fmode = READONLY;   /* from fitsio.h */
      fits_open_table(&fptr, fname, fmode, &status);
      if (status) {
	dprintf(0,"status=%d\n",status);
	fits_close_file(fptr, &status);     /* close the file, if that even makes sense */
	fits_report_error(stderr, status);  /* print out any error messages */	
	continue;
      }
      fits_get_num_rows(fptr, &nrows, &status);
      fits_get_num_cols(fptr, &ncols, &status);
      dprintf(1,"%s : Nrows: %d   Ncols: %d\n",fname,nrows,ncols);
      if (nsize>0 && nrows != nsize) {
	warning("nrows=%d nsize=%d",nrows,nsize);
	if (nsize > nrows) error("cannot continue");
      }
	
      colnames = (string *) allocate(ncols * sizeof(string));
      data_col = -1;
      for (i=0; i<ncols; i++) {    // loop over all columns to find the DATA column
	ii = i + 1;
	fits_make_keyn("TTYPE", ii, keyword, &status);
	fits_read_key(fptr, TSTRING, keyword, colname, NULL, &status);
	colnames[i] = scopy(colname);
	dprintf(1,"%s = %s\n",keyword,colname);
	if (streq(colname,"DATA")) {
	  data_col = ii;
	  fits_make_keyn("TFORM", ii, keyword, &status);
	  fits_read_key(fptr, TSTRING, keyword, data_fmt, NULL, &status);	
	  nchan = atoi(data_fmt);
	  dprintf(1,"DATA in column %d  nchan=%d\n",data_col,nchan);
	}
      } //for(i)
      dprintf(0,"%s : Nrows: %d   Ncols: %d  Nchan: %d\n",fname,nrows,ncols,nchan);

      if (ndims == 0) {

	for (k=0; k<ncolcheck; ++k) {
	  for (i=0, found=0; i<ncols; ++i) {
	    if (streq(colnames[i],colcheck[k])) {
	      found = 1;
	      fits_make_keyn("TFORM", i+1, keyword, &status);
	      fits_read_key(fptr, TSTRING, keyword, data_fmt, NULL, &status);	
	      dprintf(0,"Found %s in column %d   fmt=%s\n",colnames[i],i+1,data_fmt);
	    }
	  }
	  if (ncolcheck>0 && found==0) warning("Did not find %s",colcheck[k]);
	}
	//#define ONEDIM
#ifdef ONEDIM
	float *data1 = (float *) allocate(nchan*nrows*sizeof(float));
	float nulval = 0.0;
	fits_read_col(fptr, TFLOAT, data_col, 1, 1, nchan, &nulval, data1, &anynul, &status);
	dprintf(0,"DATA1 %g %g %g\n",data1[0],data1[1],data1[nchan]);
#else
	mdarray2 data2 = allocate_mdarray2(nrows,nchan);
	double nulval = 0.0;
	fits_read_col(fptr, TDOUBLE, data_col, 1, 1, nchan, &nulval, &data2[0][0], &anynul, &status);
	dprintf(0,"DATA2 %g %g %g\n",data2[0][0], data2[0][1], data2[1][0]);
#endif

	if (Qstats) {
	  float sum = 0.0;
	  int nmax = nchan * nrows;
	  Moment m;

	  ini_moment(&m,mom,0);
#ifdef ONEDIM
	  for (i=0; i<nmax; i++) {                  // empty loop: 1.1"
	    //sum += muladd(data[i], 1.0, 0.0);     // 2.1"
	    //sum += data2[0][i];
	    //sum += data1[i];                          // 1.6"
	    accum_moment(&m, data1[i], 1.0);         // 7.9" (6.6 if mom=1)
	  }//for(i)
#else
	  for (ii=0; ii<nrows; ii++)                 // empty loop: 1.1"
	    for (jj=0; jj<nchan; jj++)
	      if (mom==0)
		sum += data2[ii][jj];                         // 3.1"
	      else
		accum_moment(&m, data2[ii][jj], 1.0);       // 8.8" (7.3 if mom=1)
#endif
	  dprintf(0,"mean: %g\n", sum / nmax);
	  if (mom > 1)
	    printf("Mean: %g   Dispersion: %g Min: %g Max: %g   N: %d\n",
		   mean_moment(&m), sigma_moment(&m),
		   min_moment(&m), max_moment(&m), n_moment(&m));
	  else if (mom > 0)
	    printf("Mean: %g   Dispersion: %g Min: %g Max: %g   N: %d\n",
		   mean_moment(&m), -1.0,
		   min_moment(&m), max_moment(&m), n_moment(&m));
	}//if(Qstats)


	
      } else {

	if (ndims == 5) {
	  /* test case for NGC5291 with dims=2,11,2,2,4                                   */
	  warning("PS Reduction procedure in %d dimensions; data_col = %d",ndims,data_col);
	  /*                                 scan    sig     pol     int     cal     chan */

	  
	  mdarray6 data6 = allocate_mdarray6(dims[4],dims[3],dims[2],dims[1],dims[0],nchan);
	  mdarray4 data4 = allocate_mdarray4(dims[4],        dims[2],dims[1],        nchan);
	  mdarray3 data3 = allocate_mdarray3(dims[4],        dims[2],                nchan);
	  mdarray2 data2 = allocate_mdarray2(                dims[2],                nchan);
	  mdarray1 data1 = allocate_mdarray1(                                        nchan);	  
	  double nulval = 0.0;
	  dprintf(0,"DIMSIZE: %d\n",dims[4]*dims[3]*dims[2]*dims[1]*dims[0]);
	  fits_read_col(fptr, TDOUBLE, data_col, 1, 1, nchan*nsize, &nulval, &data6[0][0][0][0][0][0], &anynul, &status);
	  dprintf(0,"DATA6 %g %g %g\n",
		  data6[0][0][0][0][0][0],
		  data6[0][0][0][0][0][nchan-1],
		  data6[0][0][0][0][1][0]);


	  int i0,i1,i2,i3,i4;
	  real *s0,*s1,*s2,*s3, *s4;

	  while (bench--) {	  
	  
	  // calibration
	  for (i4=0; i4<dims[4]; ++i4) {  //scan
	    for (i2=0; i2<dims[2]; ++i2) { // pol
	      for (i1=0; i1<dims[1]; ++i1) { // int
		s0 = data6[i4][0][i2][i1][0];        // sig_calon
		s1 = data6[i4][0][i2][i1][1];        // sig_caloff
		s2 = data6[i4][1][i2][i1][0];        // ref_calon
		s3 = data6[i4][1][i2][i1][1];        // ref_caloff
		tsys = dcmeantsys(nchan, s2, s3, 1.4);
		s4 = data4[i4][i2][i1];
		ta2(nchan, s0, s1, s2, s3, tsys, s4);
	      }
	    }
	  }
	  
	  if (FALSE) {   // one shot averaging over time,scan,pol

	    // averaging over int,pol,scan
	    real **s5 = (real **) allocate(dims[1]*dims[2]*dims[4]*sizeof(real *));
	    i = 0;
	    for (i4=0; i4<dims[4]; ++i4) {  //scan
	      for (i2=0; i2<dims[2]; ++i2) { // pol
		for (i1=0; i1<dims[1]; ++i1) { // int
		  s5[i++] = data4[i4][i2][i1];
		}
	      }
	    }
	    average(dims[1]*dims[2]*dims[4],nchan,s5,data1);
	    
	  } else {  // averaging on 3 levels
	  
	    // time averaging
	    real **s5 = (real **) allocate(dims[1]*sizeof(real *));
	    for (i4=0; i4<dims[4]; ++i4) {  //scan
	      for (i2=0; i2<dims[2]; ++i2) { // pol
		for (i1=0; i1<dims[1]; ++i1) { // int
		  s5[i1] = data4[i4][i2][i1];            // invalid write
		}
		average(dims[1],nchan,s5,data3[i4][i2]);  // invalid read
	      }
	    }

	    // scan averaging
	    real **s6 = (real **) allocate(dims[4]*sizeof(real *));
	    for (i2=0; i2<dims[2]; ++i2) { // pol
	      for (i4=0; i4<dims[4]; ++i4) // scan
		s6[i4] = data3[i4][i2];
	      average(dims[4],nchan,s6,data2[i2]);
	    }
	  
	    // pol averaging
	    average2(nchan,data2[0],data2[1],data1);
	  }

	  } // bench


	  // and voila, we have a spectrum
	  real sum = 0.0;
	  for (i=0; i<nchan; i++) {
	    sum += data1[i];
	    //printf("%d %g\n",i,data1[i]);
          }
  	  printf("sum=%g\n",sum);


	  // free up
	  free_mdarray6(data6,dims[4],dims[3],dims[2],dims[1],dims[0],nchan);
	  free_mdarray4(data4,dims[4],        dims[2],dims[1],        nchan);
	  free_mdarray3(data3,dims[4],        dims[2],                nchan);
	  free_mdarray2(data2,                dims[2],                nchan);
	  free_mdarray1(data1,                                        nchan);	  
	  
	} else // ndims==5
	  warning("nothing to do for ndims=%d",ndims);
      } // ndims > 0
	  
      fits_close_file(fptr, &status);            /* close the file */
      fits_report_error(stderr, status);     /* print out any error messages */
	  
    }//for(j) loop over files
}

