/*
 * sample program to use cfitsio on SDFITS files
 * also serves as benchmark and does some low level math
 *
 * 23-nov-2019   PJT       written
 * 11-dec-2019   PJT       mdarray reduction example
 * 29-feb-2020   PJT       verbose, raw I/O
 * 28-sep-2021   PJT       better cfitsio usage, report more SDFITS properties
 *    apr-2023   PJT       some mods for the DYSH work
 *
 * Benchmark 6 N2347 files:  2.4"  (this is with mom=0 stats)
 * dims=5 for NGC5291:       31-35ms (depending in 1 or 3 levels)
 *
 * 
 * 1. Open the file, assign HDUs and convert tables to Pandas minus the DATA column  ("Load")
 * 2. For each row, create a Spectrum1D with the default GWCS, i.e. I do not create a astropy.WCS object. ("Create obsblocks")
 * 3. Remove a baseline order 1,2,3 from every row.  I don't do ON-OFF/OFF or anything like that. ("Baseline_N")
 */
#include <nemo.h>
#include <strlib.h>
#include <extstring.h>
#include <moment.h>
#include <mdarray.h>
#include <lsq.h>

#include <fitsio.h>  
#include <longnam.h>

string defv[] = {
    "in=???\n            Input SDFITS fits file(s)",
    "mom=-1\n            Adding stats to this highest moment",
    "cols=\n             Column names to track - or names of the dimensions if given",
    "dims=\n             Dimensions to reduce [ex1: 2,11,2,2,4]",
    "proc=\n             Reduction procedure (PS,FS,NOD,TP)",
    "nrows=\n            Override the nrows derived from the data",
    "nchan=\n            Override the nchan derived from the data",
    "row=-1\n            Show ascii spectrum for this row (0=first)",
    "raw=f\n             Do only raw I/O ?",
    "hdu=2\n             Which HDU (BINTABLE SINGLE DISH) to process (1=first)",
    "blfit=-1\n          Do a baseline fit of this order (-1 skips)",
    "bench=1\n           How many times to run benchmark",
    "mode=-1\n           Mode how much to process the SDFITS files (-1 means all)",
    "datadim=2\n         1: DATA[nrows*nchan]   2: DATA[nrows][nchan]",
    "VERSION=0.9g\n      19-apr-2023 PJT",
    NULL,
};

string usage = "sdfits info and bench reduction procedures";

typedef struct {
  int nchan;
  float *data;
} Spectrum1D, *Spectrum1Dptr;


real dcmeantsys(int ndata, real *calon, real *caloff, real tcal)
{
  int nedge = ndata/10;
  int i;
  real tsys, sum1 = 0.0,  sum2 = 0.0;
  static int count = 0;

  //#pragma omp for reduction ( + : sum1 )
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

  //#pragma omp parallel shared(ndata, tsys, sig, ref, ta) private(i)
  //#pragma omp for
  for (i=0; i<ndata; i++)
    ta[i] = tsys * (sig[i]/ref[i] - 1.0);
}

void ta2(int ndata, real *sig1, real *sig2, real *ref1, real *ref2, real tsys, real *ta)
{
  int i;

  //#pragma omp parallel shared(ndata, tsys, sig1, sig2, ref1, ref2, ta) private(i)
  //#pragma omp for
  for (i=0; i<ndata; i++)
    ta[i] = tsys * ((sig1[i]+sig2[i])/(ref1[i]+ref2[i]) - 1.0);
  dprintf(1,"Ta(2) %g\n",ta[0]);
}

void average(int nvec, int ndata, real **data, real *ave)
{
  int i, j;

  //#pragma omp parallel shared(ave) private(i)
  //#pragma omp for
  for (i=0; i<ndata; ++i) ave[i] = 0;
  
  //#pragma omp parallel shared(ave, data) private(i,j)
  //#pragma omp for
  for (j=0; j<nvec; ++j)
    for (i=0; i<ndata; ++i)
      ave[i] += data[j][i];  // invalid read
  
  //#pragma omp parallel shared(ave, ndata, nvec) private(i)
  //#pragma omp for
  for (i=0; i<ndata; ++i) ave[i] /= nvec;  
}

void average2(int ndata, real *data1, real *data2, real *ave)
{
  int i;

  //#pragma omp parallel shared(ave, data1, data2) private(i)
  //#pragma omp for
  for (i=0; i<ndata; i++)
    ave[i] = 0.5*(data1[i]+data2[i]);
}

/*
 *   perform a baseline subtraction
 */


void baseline(int npt, real *x, real *y, int npoly, real *coeffs, real *errors)
{
  int i,j;
  real sum;
  static real *mat = NULL;
  static real *vec = NULL;
  static real *sol = NULL;
  static real *a   = NULL;
  static real *x0  = NULL;
  static int  n0 = 0;

  if (n0 == 0) {
    dprintf(0,"baseline fit:  npt=%d\n", npt);
  }
  if (mat == NULL) {
    mat = (real *) allocate(sizeof(real) * (npoly+1)*(npoly+1));
    vec = (real *) allocate(sizeof(real) * (npoly+1));
    sol = (real *) allocate(sizeof(real) * (npoly+1));
    a   = (real *) allocate(sizeof(real) * (npoly+2));
    n0  = npt;
    x0  = (real *) allocate(sizeof(real) * n0);
    for (i=0; i<n0; i++) x0[i] = i;
  }
  if (n0 > npt) error("First npt=%d, a new npt=%d too large", n0, npt);
  if (x == NULL)   x = x0;
  
  lsq_zero(npoly+1, mat, vec);
  for (i=0; i<npt; i++) {
    a[0] = 1.0;
    for (j=0; j<npoly; j++)
	a[j+1] = a[j] * x[i];
    a[npoly+1] = y[i];
    lsq_accum(npoly+1,mat,vec,a,1.0);
  }
  //if (npoly==0) printf("TEST = %g %g\n",mat[0], vec[0]);
  lsq_solve(npoly+1,mat,vec,sol);
  //for (j=0; j<=npoly; j++) printf("%g ",sol[j]);
  //printf("\n");

  for(i=0; i<npt; i++) {
    sum=sol[npoly];
    for (j=npoly-1; j>=0; j--)
      sum = (x[i] * sum + sol[j]);
    y[i] = y[i] - sum;
  }
}

/* find a string in an array of strings; return -1 if not found */

int keyindex(int ncols, string *colnames, string keyword)
{
  for (int i=0; i<ncols; i++)
    if (streq(colnames[i],keyword)) return i;
  return -1;
}

void minmaxi(int n, int *data, int *data_min, int *data_max, string label)
{
  if (data == NULL) {
    warning("%s: no alloc n=%d", label, n);
    *data_min = *data_max = -1;
    return;
  }
  *data_min = *data_max = data[0];
  for (int i=1; i<n; i++) {
    if (data[i] < *data_min) *data_min = data[i];
    if (data[i] > *data_max) *data_max = data[i];
  }
}


//

double *get_column_dbl(fitsfile *fptr, char *colname, int nrows, int ncols, char **colnames)
{
  int col = keyindex(ncols, colnames, colname) + 1;
  if (col < 1) return NULL;
  double *data = (double *) allocate(nrows * sizeof(double));
  double nulval = 0.0;
  int anynul;
  int fstatus = 0;
  fits_read_col(fptr, TDOUBLE, col, 1, 1, nrows, &nulval, data, &anynul, &fstatus);
  if (anynul) printf("get_column %s -> %d null's\n", colname, anynul);  
  return data;
}

int *get_column_int(fitsfile *fptr, char *colname, int nrows, int ncols, char **colnames)
{
  int col = keyindex(ncols, colnames, colname) + 1;
  if (col < 1) return NULL;
  int *data = (int *) allocate(nrows * sizeof(int));
  int nulval = 0;
  int anynul;
  int fstatus = 0;
  fits_read_col(fptr, TINT, col, 1, 1, nrows, &nulval, data, &anynul, &fstatus);
  if (anynul) printf("get_column %s -> %d null's\n", colname, anynul);
  return data;
}

char **get_column_str(fitsfile *fptr, char *colname, int nrows, int ncols, char **colnames)
{
  char keyword[FLEN_KEYWORD+1], data_fmt[FLEN_VALUE+1];
  int fstatus = 0;
  int col = keyindex(ncols, colnames, colname) + 1;
  if (col < 1) return NULL;
  char **data = (char **) allocate(nrows * sizeof(char *));
  // find the string length of this keyword
  fits_make_keyn("TFORM", col, keyword, &fstatus);
  fits_read_key(fptr, TSTRING, keyword, data_fmt, NULL, &fstatus);	
  int slen = atoi(data_fmt);  // 1 extra for terminating 0
  printf("DATA in column %d  slen=%d\n",col,slen);
  // allocate the full block of chars for (slen+1)*nrows
  char *vals = (char *) allocate(nrows*(slen+1)*sizeof(char));
  for (int i=0; i<nrows; i++)
    data[i] = &vals[(slen+1)*i];
  char *nulval = "\0";
  int anynul;
  fits_read_col_str(fptr, col, 1, 1, nrows, nulval, data, &anynul, &fstatus);  
  if (anynul) printf("get_column %s -> %d null's\n", colname, anynul);  
  return data;
}

void iferror(fitsfile *fptr, int fstatus)
{
  if (fstatus) {
    printf("status=%d\n",fstatus);
    fits_close_file(fptr, &fstatus);
    fits_report_error(stderr, fstatus);
    exit(fstatus);
  }
}

#define MAXPOLYFIT 10

void nemo_main(void)
{
    fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    int status, fmode, ii, jj, i, j, k, data_col;
    long int nrows;
    long int nrows0;
    int ncols, nchan, nfiles, found, row, nhdus;
    string fname = getparam("in"), *fnames;
    string *colnames;
    int nsize, ndims, dims[MDMAXDIM];
    char keyword[FLEN_KEYWORD+1], colname[FLEN_VALUE+1], data_fmt[FLEN_VALUE+1], tdim[FLEN_VALUE+1];
    bool Qraw = getbparam("raw");
    int mom = getiparam("mom");
    int mode = getiparam("mode");
    int datadim = getiparam("datadim");
    string *colcheck = burststring(getparam("cols"),", ");
    int ncolcheck = xstrlen(colcheck, sizeof(string))-1;
    int anynul = 0;
    int bench = getiparam("bench");
    int blfit = getiparam("blfit");
    real tsys;
    real coeffs[MAXPOLYFIT], errors[MAXPOLYFIT];

    if (blfit > MAXPOLYFIT) error("blfit=%d too large; MAXPOLYFIT=%d",blfit,MAXPOLYFIT);

    fnames = burststring(fname,", ");
    nfiles = xstrlen(fnames, sizeof(string)) -1;
    dprintf(0,"Found %d files\n",nfiles);
    dprintf(0,"Found %d cols to check\n",ncolcheck);

    row = getiparam("row");

    if (hasvalue("dims")) {
      ndims = nemoinpi(getparam("dims"),dims,MDMAXDIM);
      for (i=0, nsize=1; i<ndims; i++)
	nsize *= dims[i];
      dprintf(0,"nsize=%d\n",nsize);
    } else 
      nsize = ndims = 0;

    dprintf(0,"%s mode\n", datadim==1 ? "ONEDIM" : "TWODIM");
    
    for (j=0; j<nfiles; j++) {
      fname = fnames[j];

      if (Qraw) {
	stream fp;
	char fitsline[IOBUFLEN];    // IOBUFLEN = 2880 from fitsio.h
	int ndat, n=0;
	
	fp = stropen(fname,"r");
	while(1) {
	  ndat = fread(fitsline, IOBUFLEN, 1, fp);
	  if (ndat != 1)
	    break;
	  n++;
	}
	strclose(fp);
	dprintf(0,"RAW I/O: read %d  %d FITS records = %g GB\n",n, IOBUFLEN, (n/1024.0/1024.0/1024.0)*IOBUFLEN);
	continue;
      }

      dprintf(1,"MODE=0\n");      
      status = 0;         /* initialize status before calling fitsio routines */
      fmode = READONLY;   /* from fitsio.h */
      fits_open_table(&fptr, fname, fmode, &status);
      if (status) {
	dprintf(0,"status=%d\n",status);
	fits_close_file(fptr, &status);     /* close the file, if that even makes sense */
	fits_report_error(stderr, status);  /* print out any cfitsio error messages */	
	continue;                           /* go to next fitsfile to process */
      }
      fits_get_num_hdus(fptr, &nhdus, &status);
      dprintf(1,"%s : nHDUs: %d\n",fname,nhdus);
      fits_get_num_rows(fptr, &nrows, &status);
      fits_get_num_cols(fptr, &ncols, &status);
      dprintf(1,"%s : nRows: %d   nCols: %d\n",fname,nrows,ncols);
      if (nsize>0 && nrows != nsize)
	warning("nrows=%d nsize=%d",nrows,nsize);

      if (mode >= 0 && mode < 1) continue;
      dprintf(1,"MODE=1\n");
		
      colnames = (string *) allocate(ncols * sizeof(string));
      data_col = -1;               // need to find the DATA column
      for (i=0; i<ncols; i++) {    // loop over all columns to find the DATA column, awkward
	ii = i + 1;
	fits_make_keyn("TTYPE", ii, keyword, &status);
	fits_read_key(fptr, TSTRING, keyword, colname, NULL, &status);
	colnames[i] = scopy(colname);
	fits_make_keyn("TFORM", ii, keyword, &status);
	fits_read_key(fptr, TSTRING, keyword, data_fmt, NULL, &status);
	// DATA is special since we need to get the number of channels from TFORM
	if (streq(colname,"DATA") || streq(colname,"SPECTRUM")) {  // classic SDFITS or CLASS FITS
	  data_col = ii;
	  nchan = atoi(data_fmt);
	  dprintf(2,"DATA in column %d  nchan=%d\n",data_col,nchan);
	}
	dprintf(1,"%3d:  %-16s %s\n",ii,colname,data_fmt);
      } //for(i)
      real np = nrows * nchan * 4 / 1e6;
      dprintf(0,"%s : Nrows: %d   Ncols: %d  Nchan: %d  nP: %g Mp\n",fname,nrows,ncols,nchan,np);
      int col_tcal  = keyindex(ncols, colnames, "TCAL") + 1;
      int col_cal   = keyindex(ncols, colnames, "CAL") + 1;
      int col_sig   = keyindex(ncols, colnames, "SIG") + 1;
      int col_fdnum = keyindex(ncols, colnames, "FDNUM") + 1;
      int col_ifnum = keyindex(ncols, colnames, "IFNUM") + 1;
      int col_plnum = keyindex(ncols, colnames, "PLNUM") + 1;
      dprintf(0,"tcal: %d   cal: %d sig: %d fdnum: %d ifnum: %d plnum: %d\n",
	      col_tcal, col_cal, col_sig, col_fdnum, col_ifnum, col_plnum);


      // get optional columns
      int *fdnum_data = get_column_int(fptr, "FDNUM", nrows, ncols, colnames);
      int *ifnum_data = get_column_int(fptr, "IFNUM", nrows, ncols, colnames);
      int *plnum_data = get_column_int(fptr, "PLNUM", nrows, ncols, colnames);
      int *int_data   = get_column_int(fptr, "INT",   nrows, ncols, colnames);
      
      int fd_min, fd_max, if_min, if_max, pl_min, pl_max, int_min, int_max;
      minmaxi(nrows, fdnum_data, &fd_min, &fd_max, "fdnum");
      minmaxi(nrows, ifnum_data, &if_min, &if_max, "ifnum");
      minmaxi(nrows, plnum_data, &pl_min, &pl_max, "plnum");
      printf("FDNUM: %d %d\n", fd_min, fd_max);
      printf("IFNUM: %d %d\n", if_min, if_max);
      printf("PLNUM: %d %d\n", pl_min, pl_max);

      if (int_data) {
	minmaxi(nrows, int_data, &int_min, &int_max, "int");
	printf("INT:   %d %d\n", int_min, int_max);	
      }
      char **sig_data = get_column_str(fptr, "SIG", nrows, ncols, colnames);
      char **cal_data = get_column_str(fptr, "CAL", nrows, ncols, colnames);
      int nsig=0, ncal=0;
      if (sig_data && cal_data) {
	nsig = ncal = 1;
	for (i=1; i<nrows; i++)
	  if (sig_data[i][0] != sig_data[0][0]) {
	    printf("ODD SIG %d = %s %s\n",i,sig_data[0],sig_data[i]);
	    nsig++;
	    break;
	  }
	for (i=1; i<nrows; i++)
	  if (cal_data[i][0] != cal_data[0][0]) {
	    printf("ODD CAL %d = %s %s\n",i,cal_data[0],cal_data[i]);	    
	    ncal++;
	    break;
	  }
      }
      printf("SIG:   %d  = %c\n", nsig, nsig==1 ? sig_data[0][0] : ' ');
      printf("CAL:   %d  = %c\n", ncal, ncal==1 ? cal_data[0][0] : ' ');

      if (mode >= 0 && mode < 2) continue;            
      dprintf(1,"MODE=2\n");
	
      if (0) {
	for (i=0; i<nrows; i++)
	  printf("SIGCAL %d '%s' '%s'\n", i, sig_data[i], cal_data[i]);
      }

      // printf("INT: @ 0x%d\n", int_data);
      
      if (hasvalue("nchan")) {
	nchan = getiparam("nchan");
	warning("Overriding with nchan=%d",nchan);
      }
      
      char **tdim_data = (char **) allocate(nrows * sizeof(char **));
      fits_make_keyn("TDIM", data_col, keyword, &status);      
      int col_tdim = keyindex(ncols, colnames, keyword) + 1;
      char *nulvals = "\0";
      tdim_data[0] = tdim;
      fits_read_col_str(fptr, col_tdim, 1, 1, 1, nulvals, tdim_data, &anynul, &status);
      dprintf(0,"DATA keyword %s = %s\n",keyword, tdim);

      if (mode >= 0 && mode < 3) continue;            
      dprintf(1,"MODE=3\n");
	
      if (Qraw) ndims = -1;
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
       
	if (mode >= 0 && mode < 4) continue;
	dprintf(1,"MODE=4\n");

	/*
	 *    Now get at the big DATA[] - this is where virtually memory gets filled
	 *    at the data I/O rate
	 */

	float *data1 = NULL;     // data1[nrows*nchan]
	mdarray2 data2 = NULL;   // data2[nrows][nchan]
	
	if (datadim==1) {
	  // just one long 1dim array. Force single precision
	  dprintf(0,"ONEDIM DATA1: make Waterfall\n");
	  data1 = (float *) allocate(nchan*nrows*sizeof(float));
	  float nulval = 0.0;
	  fits_read_col(fptr, TFLOAT, data_col, 1, 1, nchan*nrows, &nulval, data1, &anynul, &status);
	  if (nrows > 1)
	    dprintf(0,"DATA1 %g %g %g\n",data1[0],data1[1],data1[nchan]);
	  else
	    dprintf(0,"DATA1 %g %g ... %g (only 1 row)\n",data1[0],data1[1],data1[nchan-1]);
	  if (blfit >= 0) warning("ONEDIM mode cannot subtract baselines");  // also due to float/real
	} else {
	  // 2dim array, suitable for waterfall plot. Note this could be double precision in NEMO
	  dprintf(0,"TWODIM DATA2: make Waterfall\n");
	  data2 = allocate_mdarray2(nrows,nchan);
	  double nulval = 0.0;
	  fits_read_col(fptr, TDOUBLE, data_col, 1, 1, nchan*nrows, &nulval, &data2[0][0], &anynul, &status);
	  if (nrows > 1)
	    dprintf(0,"DATA2 %g %g %g\n",data2[0][0], data2[0][1], data2[1][0]);
	  else
	    dprintf(0,"DATA2 %g %g ... %g (only 1 row)\n",data2[0][0], data2[0][1], data2[0][nchan-1]);
	  if (blfit >= 0) {
	    dprintf(0,"BASELINE %d\n",blfit);
	    for (ii=0; ii<nrows; ii++)
	      baseline(nchan, NULL, data2[ii], blfit, coeffs, errors);
	  }
	}

	if (mode >= 0 && mode < 5) continue;
	dprintf(1,"MODE=5\n");
	

	if (mom >= 0) {
	  dprintf(0,"Stats\n");
	  float sum = 0.0;
	  int nmax = nchan * nrows;
	  Moment m;

	  ini_moment(&m,mom,0);
	  if (datadim==1) {
	    for (i=0; i<nmax; i++) {                  // empty loop: 1.1"
	      if (mom==0)
		sum += data1[i];                          // 1.6"
	      else
		accum_moment(&m, data1[i], 1.0);         // 7.9" (6.6 if mom=1)
	    }//for(i)
	  } else {
	    for (ii=0; ii<nrows; ii++)                 // empty loop: 1.1"
	      for (jj=0; jj<nchan; jj++)
		if (mom==0)
		  sum += data2[ii][jj];                         // 3.1"
		else
		  accum_moment(&m, data2[ii][jj], 1.0);       // 8.8" (7.3 if mom=1)
	    if (row >= 0) {
	      warning("spectrum row");
	      for (jj=0; jj<nchan; jj++)
		printf("%d %g\n",jj,data2[row][jj]);
	    }
	  } 
	  dprintf(0,"mean: %g\n", sum / nmax);
	  if (mom > 2)
	    printf("MOM %d TBD\n",mom);
	  else if (mom > 1)
	    printf("Mean: %g   Dispersion: %g Min: %g Max: %g   N: %d\n",
		   mean_moment(&m), sigma_moment(&m),
		   min_moment(&m), max_moment(&m), n_moment(&m));
	  else if (mom > 0)
	    printf("Mean: %g   Dispersion: %g Min: %g Max: %g   N: %d\n",
		   mean_moment(&m), -1.0,
		   min_moment(&m), max_moment(&m), n_moment(&m));
	} // stats


	
      } else if (ndims > 0) {     // special processing

	if (hasvalue("nchan")) {
	  nchan = getiparam("nchan");
	  warning("Re-using the dims analysis with nchan=%d",nchan);
	}

	if (ndims == 5) {
	  /* test case for NGC5291 with dims=2,11,2,2,4                                   */
	  warning("PS Reduction procedure in %d dimensions; data_col = %d",ndims,data_col);
	  /*                                 scan    sig     pol     int     cal     chan */
	  mdarray6 data6 = allocate_mdarray6(dims[4],dims[3],dims[2],dims[1],dims[0],nchan);
	  mdarray4 data4 = allocate_mdarray4(dims[4],        dims[2],dims[1],        nchan);
	  mdarray3 data3 = allocate_mdarray3(dims[4],        dims[2],                nchan);
	  mdarray2 data2 = allocate_mdarray2(                dims[2],                nchan);
	  mdarray1 data1 = allocate_mdarray1(                                        nchan);
	  double  *tcal  = (double *) allocate(nrows*sizeof(double));
	  double nulval = 0.0;
	  dprintf(0,"DIMSIZE: %d\n",dims[4]*dims[3]*dims[2]*dims[1]*dims[0]);
	  // make sure nsize <= nrows
	  fits_read_col(fptr, TDOUBLE, data_col, 1, 1, nchan*nsize, &nulval, &data6[0][0][0][0][0][0], &anynul, &status);
	  if (col_tcal < 0)
	    for(i=0; i<nsize; i++) tcal[i] = 1.0;
	  fits_read_col(fptr, TDOUBLE, col_tcal, 1, 1,       nsize, &nulval, tcal, &anynul, &status);	  
	  // should also read the TCAL column
	  dprintf(0,"DATA6 %g %g %g     %g\n",
		  data6[0][0][0][0][0][0],
		  data6[0][0][0][0][0][nchan-1],
		  data6[0][0][0][0][1][0],
		  tcal[0]);


	  int i1,i2,i4;
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
		tsys = dcmeantsys(nchan, s2, s3, 1.42);
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
	} else if (ndims == 7) {	  // NOD example
	  warning("NOD masochism example here?");
	  /*                                   scan    sig     fd      if,     pol     int     cal     chan */
	  //mdarray8 data8 = allocate_mdarray8(dims[6],dims[5],dims[4],dims[3],dims[2],dims[1],dims[0],nchan);
	} else { // ndims==5
	  warning("nothing to do for ndims=%d",ndims);
	}
      } // ndims > 0
	  
      fits_close_file(fptr, &status);            /* close the file */
      fits_report_error(stderr, status);     /* print out any error messages */

	  
    } //for(j) loop over files
}

