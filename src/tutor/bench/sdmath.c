/*
 * benchmark some SD match
 *
 * The default nscan=1000 nchan=100000 iter=20 needs 14 Gop
 * d76:    4.6s -> 3 Gop/sec128
 * lma:    7.3s    2
 * jansky: 3.4s    4
 * 
 *     30-apr-2025    Created  // PJT
 */
 
#include <nemo.h>


string defv[] = {
  "nscan=1000\n        Number of scans",
  "nchan=100000\n      Number of channels",
  "iter=20\n           How many times to iterate (0 means do nothing)",
  "mode=1\n            Mode of math (not implemented)",
  "aver=f\n            Add time average over nscan",
  "in=\n               Read a file into memory",
  "bs=16\n             Blocksize in kB to read",
  "maxbuf=0\n          Max buffer for file data in MB",
  "VERSION=1.4\n       2-may-2025",
  NULL,
};

string usage = "SD math bench";


// n OPS
real mean1(int n, real *data)
{
  real sum = 0.0;
  #pragma omp parallel for reduction(+:sum)
  for (int i=0; i<n; i++) sum += data[i];
  return sum; 
}

// 2*n OPS
real mean2(int n, real *data1, real *data2)
{
  real sum = 0.0;
  #pragma omp parallel for reduction(+:sum)
  for (int i=0; i<n; i++) sum += data1[i]-data2[i];
  return sum; 
}




// this one went slower!!!
// 2*n OPS
real tsys2(int n, real *hot, real *cold)
{
  real a, b;
  // @todo these two can be done in parallel
#pragma omp parallel
  {
#pragma omp sections 
    {
#pragma omp section
      a = mean1(n, hot);
#pragma omp section
      b = mean2(n, hot, cold);
    }
  }
  return a/b;
}

real tsys(int n, real *hot, real *cold)
{
  // @todo these two can be done in parallel
  return mean1(n, hot) / mean2(n, hot, cold);
}

void file_read(string fname)
{
  size_t nread, n = nemo_file_size(fname);
  char *data;
  size_t i, nreset=0;
  off_t o;
  int bs = getiparam("bs");   // "bs" in kB
  size_t bufsize = bs*1024;   
  size_t maxbuf = 1024*1024*getiparam("maxbuf"); // "maxbuf" was in MB

  if (maxbuf == 0) {
    data = allocate(n);
    maxbuf = n;
  } else
    data = allocate(maxbuf);

  printf("file size: %ld bytes ~ %g MB   bufsize=%ld ~ %g MB\n",n, (real)n/1024/1024, bufsize, (real)bufsize/1024/1024);
  printf("maxbuf:    %ld bytes = %d MB\n",maxbuf, getiparam("maxbuf"));
  
  stream fp = stropen(fname,"r");
  for(i=0, o=0; ; i++) {
    //o=0;
    nread = fread(&data[o], bufsize, 1, fp);
    if (nread<=0) break;
    o += bufsize;
    if (o >= maxbuf) {
      nreset++;
      o = 0;
      dprintf(1,"buffer RESET at i=%ld\n",i+1);
    } else
      dprintf(2,"buffer at i=%ld at o=%ld maxbuf=%ld\n",i+1,o,maxbuf);
  }
  strclose(fp);
  printf("Last block# %ld of %ld,  #resets %ld\n", i, bufsize, nreset);
}



void nemo_main(void)
{
  size_t nchan = getiparam("nchan");
  size_t nscan = getiparam("nscan");
  size_t iter = getiparam("iter");
  size_t n = nchan * nscan * 4;
  int i, j, mode = getiparam("mode");
  real *data, *row0, *row1, *row2, *row3, *row4;
  real *aver, *spec, *onn, *off, t_onn, t_off;
  size_t ndata = n * sizeof(real);
  bool Qaver = getbparam("aver");
  int nphase = 4;

  if (hasvalue("in")) {
    file_read(getparam("in"));
    return;
  }

  data = (real *) allocate(n * sizeof(real));
  onn  = (real *) allocate(nchan * sizeof(real));
  off  = (real *) allocate(nchan * sizeof(real));
  spec = (real *) allocate(nchan * sizeof(real));
  aver = (real *) allocate(nchan * sizeof(real));
  

  printf("data size: %ld bytes ~ %g MB  nscan=%ld nchan=%ld iter=%ld\n",
	 ndata, (real)ndata/1024/1024, nscan, nchan, iter);
  printf("# ops = %ld  %g Gop\n",7*nchan*nscan*iter,  (7.0*nchan*nscan*iter)/1e9);
    

  // only random one row, xrandom() is expensive
  for (i=0; i<nphase*nchan; i++) data[i] = xrandom(0.0,1.0);
  // duplicate other rows, so we fill memory out of the caching
  row0 = data;
  row1 = data + nphase*nchan;
  for (j=1; j<nscan; j++) {
    for (i=0; i<nphase*nchan; i++) row1[i] = row0[i];
    row0 = row1;
    row1 = row0 + nphase*nchan;
  }

  
  while (iter--) {
    row0 = data;
    for (i=0; i<nscan; i++) {
      row1 = row0 + nchan;
      row2 = row1 + nchan;
      row3 = row2 + nchan;
      row4 = row3 + nchan;
      // Tsys =  Tc * <cold> / <hot-cold>
      t_onn = tsys(nchan, row1, row2);          // 3*nchan
      // t_off = tsys(nchan, row3, row4);       // 3*nchan
      // Ta = Tsys * (on/off-1)  - but on and off are their Calon&Caloff averaged

      #pragma omp parallel for
      for (j=0; j<nchan; j++) {
	spec[j] = t_onn*((row1[j]+row2[j])/(row3[j]+row4[j])-1);  // 4*nchan      3.90user 0.89system 0:04.79elapsed 99%CPU
	//spec[j] = ((row1[j]+row2[j])/(row3[j]+row4[j])-1);      // check-1      3.85user 0.90system 0:04.76elapsed 99%CPU
	//spec[j] = (row1[j]+row2[j])/(row3[j]+row4[j]);          // check-2      3.81user 0.88system 0:04.70elapsed 99%CPU 
	//spec[j] = row1[j]/row3[j];                              // check-3      3.77user 0.91system 0:04.69elapsed 99%CPU 
	//spec[j] = row1[j];                                      // check-4      3.37user 0.88system 0:04.26elapsed 99%CPU 
       	if (Qaver)
	  aver[j] = aver[j] + spec[j];
      } // j
    } // i
  } // iter

  dprintf(1,"data[first..last] = %g %g\n", data[0], data[n-1]);
  dprintf(1,"aver=%g\n",mean1(nchan,aver));
  
}
