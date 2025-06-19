/*
 * benchmark some SD : math and I/O
 *
 * The default nscan=1000 nchan=100000 iter=10
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
  "iter=10\n           How many times to iterate (0 means do nothing)",
  "mode=1\n            Mode to init the data (0=random 1=arange)",
  "aver=f\n            Add time average over nscan",
  "in=\n               Read a file into memory",
  "bs=16\n             Blocksize in kB to read",
  "maxbuf=0\n          Max buffer for file data in MB",
  "seed=0\n            Seed, if need be",
  "VERSION=1.5\n       5-may-2025",
  NULL,
};

string usage = "SD math bench";


//   assumed order:   on_cold, on_hot, off_cold, off_hot
//
//   Tsys =  Tc * <row1> / <row2-row1>      <cold> / <hot>-<cold>
//   on  = row1 + row2
//   off = row3 + row4
//   Ta = Tsys * (on/off - 1)
//   Spectrum = <Ta>


// n OPS
real sum1(int n, real *data)
{
  real sum = 0.0;
  #pragma omp parallel for reduction(+:sum)
  for (int i=0; i<n; i++) sum += data[i];
  return sum; 
}

// 2*n OPS
real sum2(int n, real *data1, real *data2)
{
  real sum = 0.0;
  #pragma omp parallel for reduction(+:sum)
  for (int i=0; i<n; i++) sum += data2[i]-data1[i];
  return sum; 
}




// this one went slower!!!
// 2*n OPS
real tsys2(int n, real *cold, real *hot)
{
  real a, b;
  // @todo these two can be done in parallel
#pragma omp parallel
  {
#pragma omp sections 
    {
#pragma omp section
      a = sum1(n, cold);
#pragma omp section
      b = sum2(n, cold, hot);
    }
  }
  return a/b;
}

real tsys(int n, real *cold, real *hot)
{
  static real tc = 10.0;
  // @todo these two can be done in parallel
  return tc * sum1(n, cold) / sum2(n, cold, hot);
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
  size_t n4 = nchan * nscan * 4;
  int i, j, mode = getiparam("mode");
  real *data, *row0, *row1, *row2, *row3, *row4;
  real *aver, *spec, t_onn;
  size_t ndata = n4 * sizeof(real);
  bool Qaver = getbparam("aver");
  int nphase = 4;

  if (hasvalue("in")) {
    file_read(getparam("in"));
    return;
  }

  data = (real *) allocate(n4 * sizeof(real));
  spec = (real *) allocate(nchan * sizeof(real));
  aver = (real *) allocate(nchan * sizeof(real));
  

  printf("data size: %ld bytes ~ %g MB  nscan=%ld nchan=%ld iter=%ld\n",
	 ndata, (real)ndata/1024/1024, nscan, nchan, iter);
  printf("# rops = %ld  ~ %g Grops\n",nchan*nscan*iter,  (nchan*nscan*iter)/1e9);
    

  // only random one scan, xrandom() is expensive
  if (mode==1) {
    int seed = init_xrandom(getparam("seed"));
    dprintf(1,"seed used = %d\n",seed);
    for (i=0; i<nphase*nchan; i++) data[i] = xrandom(0.0,1.0);
  } else {
    // simple 0,1,2 for reproducing python
    for (i=0; i<nphase*nchan; i++) data[i] = i;
  }
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
      row1 = row0;
      row2 = row1 + nchan;
      row3 = row2 + nchan;
      row4 = row3 + nchan;
      row0 = row4;  // for the next iter
      // Tsys =  Tc * <cold> / <hot-cold>
      t_onn = tsys(nchan, row1, row2);          // 3*nchan

      #pragma omp parallel for
      for (j=0; j<nchan; j++) {
	spec[j] = t_onn*((row1[j]+row2[j])/(row3[j]+row4[j])-1);  // 4*nchan      3.90user 0.89system 0:04.79elapsed 99%CPU
	//spec[j] = ((row1[j]+row2[j])/(row3[j]+row4[j])-1);      // check-1      3.85user 0.90system 0:04.76elapsed 99%CPU
	//spec[j] = (row1[j]+row2[j])/(row3[j]+row4[j]);          // check-2      3.81user 0.88system 0:04.70elapsed 99%CPU 
	//spec[j] = row1[j]/row3[j];                              // check-3      3.77user 0.91system 0:04.69elapsed 99%CPU 
	//spec[j] = row1[j];                                      // check-4      3.37user 0.88system 0:04.26elapsed 99%CPU 
       	if (Qaver)
	  aver[j] = aver[j] + spec[j];
      } // j-chan
    } // i-scan
  } // iter

  dprintf(1,"data[0,1..last] = %g %g %g\n", data[0], data[1], data[n4-1]);
  dprintf(1,"sum=%g\n",sum1(nchan,aver));
  
}
