/*
 * benchmark SD match
 * 
 *     30-apr-2025    Created  // PJT
 */
 
#include <nemo.h>


string defv[] = {
  "nscan=100\n         Number of scans",
  "nchan=32768\n       Number of channels",
  "mode=1\n            Mode of math",
  "iter=1\n            How many times to iterate (0 means do nothing)",
  "VERSION=1\n         30-apr-2013",
  NULL,
};

string usage = "SD bench";


real mean1(int n, real *data)
{
  real sum = 0.0;
  #pragma omp parallel for reduction(+:sum)
  for (int i=0; i<n; i++) sum += data[i];
  return sum;  //  no need to normalize, since they divide out later
}

real mean2(int n, real *data1, real *data2)
{
  real sum = 0.0;
  #pragma omp parallel for reduction(+:sum)
  for (int i=0; i<n; i++) sum += data1[i]-data2[i];
  return sum;  //  no need to normalize, since they divide out later
}


real tsys(int n, real *hot, real *cold)
{
  int i;
  return mean1(n, hot) / mean2(n, hot, cold);
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
  real *spec, *onn, *off, t_onn, t_off;
  size_t ndata = n * sizeof(real);

  data = (real *) allocate(n * sizeof(real));
  onn  = (real *) allocate(nchan * sizeof(real));
  off  = (real *) allocate(nchan * sizeof(real));
  spec = (real *) allocate(nchan * sizeof(real));  
  

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
      t_onn = tsys(nchan, row1, row2);
      t_off = tsys(nchan, row3, row4);
      // Ta = tsys * (on/off-1)
      #pragma omp parallel for
      for (j=0; j<nchan; j++)
	spec[j] = t_onn*(row1[j]/row3[j]-1);
    }
  }

}
