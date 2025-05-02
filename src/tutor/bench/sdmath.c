/*
 * benchmark some SD match
 * 
 *     30-apr-2025    Created  // PJT
 */
 
#include <nemo.h>


string defv[] = {
  "nscan=100\n         Number of scans",
  "nchan=32768\n       Number of channels",
  "iter=1\n            How many times to iterate (0 means do nothing)",
  "mode=1\n            Mode of math (not implemented)",
  "aver=f\n            Add time average over nscan",
  "in=\n               Read a file into memory",
  "bs=16\n             Blocksize in kB to read",
  "VERSION=1.3\n       1-may-2025",
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
  printf("file size: %ld bytes ~ %g MB\n",n, (real)n/1024/1024);  
  char *data = allocate(n);
  size_t i;
  off_t o;
  int bs = getiparam("bs");
  int bufsize = bs*1024;   // bs in kB
  
  stream fp = stropen(fname,"r");
  for(i=0, o=0; ; i++) {
    nread = fread(&data[o], bufsize, 1, fp);
    if (nread<=0) break;
    o += bufsize;
  }
  strclose(fp);
  printf("Last block %ld of %d\n", i, bufsize);
}



void nemo_main(void)
{
  int i, j;
  size_t nchan = getiparam("nchan");
  size_t nscan = getiparam("nscan");
  int mode = getiparam("mode");
  int iter = getiparam("iter");
  size_t n = nchan * nscan * 4;
  real *data, *row0, *row1, *row2, *row3, *row4;
  real *aver, *spec, *onn, *off, t_onn, t_off;
  size_t ndata = n * sizeof(real);
  bool Qaver = getbparam("aver");

  if (hasvalue("in")) {
    file_read(getparam("in"));
    return;
  }

  data = (real *) allocate(n * sizeof(real));
  onn  = (real *) allocate(nchan * sizeof(real));
  off  = (real *) allocate(nchan * sizeof(real));
  spec = (real *) allocate(nchan * sizeof(real));
  aver = (real *) allocate(nchan * sizeof(real));
  

  printf("data size: %ld bytes ~ %g MB\n",ndata, (real)ndata/1024/1024);
    

  // only random one row, xrandom() is expensive
  for (i=0; i<nchan*4; i++) data[i] = xrandom(0.0,1.0);
  // duplicate other rows, so we fill memory out of the caching
  row0 = data;
  row1 = data + 4*nchan;
  for (j=1; j<nscan; j++) {
    for (i=0; i<nchan*4; i++) row1[i] = row0[i];
    row0 = row1;
    row1 = row0 + 4*nchan;
  }

  
  while (iter--) {
    row0 = data;
    for (i=0; i<nscan; i++) {
      row1 = row0 + nchan;
      row2 = row1 + nchan;
      row3 = row2 + nchan;
      row4 = row3 + nchan;
      // tsys = <hot> / <hot-cold>
      t_onn = tsys(nchan, row1, row2);          // 3*nchan
      // t_off = tsys(nchan, row3, row4);       // 3*nchan
      // Ta = tsys * (on/off-1)  - but on and off are their Calon&Caloff averaged

      #pragma omp parallel for
      for (j=0; j<nchan; j++) {
	spec[j] = t_onn*((row1[j]+row2[j])/(row3[j]+row4[j])-1);    // 4*nchan
        //spec[j] = t_onn*(row1[j]/row3[j]-1);    // 1*nchan
#if 1
       	if (Qaver)
	  aver[j] = aver[j] + spec[j];
#endif	
      } // j
    } // i
  } // iter
  
  dprintf(1,"aver=%g\n",mean1(nchan,aver));
}


// 1st version:   7n
// 2nd version:   7n
