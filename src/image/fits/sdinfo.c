/*
 * sample program to use cfitsio on SDFITS files
 *
 * 23-nov-2019   PJT       written
 * 11-dec-2019
 *
 * Benchmark 6 N2347 files:  2.4"
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
    "stats=t\n           Doing simple stats",
    "mom=2\n             Highest moment for stats",
    "cols=\n             Column names to track",
    "VERSION=0.3\n       11-dec-2019 PJT",
    NULL,
};

string usage = "sdfits info and bench";


  // inline
float muladd(float x, float a, float b) {
    return x*a + b;
}

void nemo_main(void)
{
    fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    int status, fmode, ii, jj, i, j, k, colnum,data_col;
    long int nrows;
    int ncols, nchan, nfiles, found;
    string fname = getparam("in"), *fnames;
    string *colnames;
    char template[64];
    char keyword[FLEN_KEYWORD], colname[FLEN_VALUE], data_fmt[FLEN_VALUE];
    bool Qstats = getbparam("stats");
    int mom = getiparam("mom");
    string *colcheck = burststring(getparam("cols"),", ");
    int ncolcheck = xstrlen(colcheck, sizeof(string))-1;

    fnames = burststring(fname,", ");
    nfiles = xstrlen(fnames, sizeof(string)) -1;
    dprintf(0,"Found %d files\n",nfiles);
    dprintf(0,"Found %d cols to check\n",ncolcheck);

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
      
      int anynul = 0;

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
      fits_close_file(fptr, &status);            /* close the file */
      fits_report_error(stderr, status);     /* print out any error messages */

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
	  !accum_moment(&m, data1[i], 1.0);         // 7.9" (6.6 if mom=1)
	}//for(i)
#else	
	for (ii=0; ii<nrows; ii++) {                  // empty loop: 1.1"
	  for (jj=0; jj<nchan; jj++) {
	    //sum += data2[ii][jj];                         // 3.1"
	    accum_moment(&m, data2[ii][jj], 1.0);       // 8.8" (7.3 if mom=1)
	  }
	}
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
    }//for(j) loop over files
}
