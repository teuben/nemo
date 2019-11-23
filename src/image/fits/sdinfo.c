/*
 * sample program to use cfitsio on SDFITS files
 *
 * 23-nov-2019   PJT       written
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

string defv[] = {
    "in=???\n            Input SDFITS  fits file(s)",
    "stats=t\n           Doing simple stats",
    "VERSION=0.2\n       23-nov-2019 PJT",
    NULL,
};

string usage = "sdfits info";

#if defined(HAVE_LIBCFITSIO)
# include <fitsio.h>  
# include <longnam.h>
#else
#error CFITSIO not available or properly configured
#endif

float muladd(float x, float a, float b) {
    return x*a + b;
}

void nemo_main(void)
{
    fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    int status, fmode, ii, jj, i, j, colnum,data_col;
    long int nrows;
    int ncols, nchan, nfiles;
    string fname = getparam("in"), *fnames;
    string *colnames;
    char template[64];
    char keyword[FLEN_KEYWORD], colname[FLEN_VALUE], data_fmt[FLEN_VALUE];
    bool Qstats = getbparam("stats");

    fnames = burststring(fname,", ");
    nfiles = xstrlen(fnames, sizeof(string)) -1;
    dprintf(0,"Found %d files\n",nfiles);

    for (j=0; j<nfiles; j++) {
      fname = fnames[j];
      status = 0;         /* initialize status before calling fitsio routines */
      fmode = READONLY;   /* from fitsio.h */
      fits_open_table(&fptr, fname, fmode, &status);
      fits_get_num_rows(fptr, &nrows, &status);
      fits_get_num_cols(fptr, &ncols, &status);
      dprintf(1,"%s : Nrows: %d   Ncols: %d\n",fname,nrows,ncols);
      colnames = (string *) allocate(ncols * sizeof(string));
      data_col = -1;
      for (i=0; i<ncols; i++) {
	ii = i + 1;
	fits_make_keyn("TTYPE", ii, keyword, &status);
	fits_read_key(fptr, TSTRING, keyword, colname, NULL, &status);
	//colnames[i] = scopy(colname);
	dprintf(1,"%s = %s\n",keyword,colname);
	if (streq(colname,"DATA")) {
	  data_col = ii;
	  fits_make_keyn("TFORM", ii, keyword, &status);
	  fits_read_key(fptr, TSTRING, keyword, data_fmt, NULL, &status);	
	  nchan = atoi(data_fmt);
	  dprintf(1,"DATA in column %d  nchan=%d\n",data_col,nchan);
	}
      }
      dprintf(0,"%s : Nrows: %d   Ncols: %d  Nchan: %d\n",fname,nrows,ncols,nchan);      
      int anynul = 0;
#if 1   
      float *data = (float *) allocate(nchan*nrows*sizeof(float));
      float nulval = 0.0;
      fits_read_col(fptr, TFLOAT, data_col, 1, 1, nchan, &nulval, data, &anynul, &status);
#else
      mdarray2 data = allocate_mdarray2(nrows,nchan);
      double nulval = 0.0;
      fits_read_col(fptr, TDOUBLE, data_col, 1, 1, nchan, &nulval, data, &anynul, &status);
    
#endif
      fits_close_file(fptr, &status);            /* close the file */
      fits_report_error(stderr, status);  /* print out any error messages */

      if (Qstats) {
	float sum = 0.0;
	int nmax = nchan * nrows;
	for (i=0; i<nmax; i++) {
	  sum += muladd(data[i], 1.0, 0.0);
	  //sum += data[0][i];	
	}
	dprintf(0,"mean: %g\n", sum / nmax);
      }
    }
}
